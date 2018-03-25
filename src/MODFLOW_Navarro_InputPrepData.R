## MODFLOW_Navarro_InputPrepData.R
#' This script is intended to prepare input files for MODFLOW
#' from various geospatial datasets, and save them as text files.
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source("src/paths+packages.R")

# X and Y resolution
DELR <- 250  # x spacing - column width
DELC <- 250  # y spacing - row height

# what variables to process?
ibound <- T
ztop <- F
riv <- T
wel <- T

## load common data
# domain boundary shapefile
shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent", stringsAsFactors=F)
shp.adj@data$HUC12 <- as.numeric(shp.adj@data$HUC12)
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

# river boundary shapefile
shp.riv.adj <- readOGR(dsn="data/NHDPlusV2/HYD", layer="NHDPlusV21_National_Seamless NHDFlowline_Network_WBDHU12_Navarro+Adjacent",
                       stringsAsFactors=F)
shp.riv.adj@data$StreamOrde <- as.numeric(shp.riv.adj@data$StreamOrde)

# Make boundary conditions (IBOUND) ---------------------------------------

if (ibound){
  
  # make empty raster
  r.empty <- raster(extent(shp.adj.UTM), crs = crs.MODFLOW, resolution=c(DELR, DELC))
  
  # rasterize
  r.ibound.in <- rasterize(shp.adj.UTM, r.empty, field="HUC12")
  r.ibound <- trim(r.ibound.in)
  
  # find coastline (StreamOrde = -9) for constant head
  shp.riv.adj.coast <- subset(shp.riv.adj, StreamOrde == -9)
  
  # reproject to UTM
  shp.riv.adj.coast.UTM <- spTransform(shp.riv.adj.coast, crs.MODFLOW)
  
  # rasterize
  r.coast <- rasterize(shp.riv.adj.coast.UTM, r.ibound, field="StreamOrde", fun=mean)
  
  # set constant head boundaries
  n.coast.temp <- -999  # placeholder value for coast
  r.ibound[r.coast==-9] <- n.coast.temp
  
  # extract as matrix
  m.ibound <- as.matrix(r.ibound)
  
  # in each row, convert all cells to the left of a coastline to inactive
  setCoast <- function(x, n.coast=n.coast.temp){
    if (sum(x==n.coast, na.rm=T) > 0){
      i.coast <- min(which(x==n.coast))
      x[1:(i.coast-1)] <- NA
      x[(i.coast+1):length(x)] <- abs(x[(i.coast+1):length(x)])  # convert any additional coast pixels to 1s
    }
    return(x)
  }
  m.ibound <- t(apply(m.ibound, 1, setCoast))
  
  # there are two little disconnected zones along the coastline (north and south ends) 
  # following trimming which should be eliminated to avoid convergence issues
  #  their HUC12 numbers are:
  #   180101080803 (north end - delete 0-1300 m)
  #   180101080905 (south end - delete 37300 - end m)
  i.top <- 1:round(1300/DELC)
  i.bot <- round(37300/DELC):dim(m.ibound)[1]
  num.replace <- function(x, to.replace, replace.with){
    x[x %in% to.replace] <- replace.with
    return(x)
  }
  
  m.ibound[i.top,] <- t(apply(m.ibound[i.top,], 1, num.replace, to.replace=c(n.coast.temp, abs(n.coast.temp), 180101080803), replace.with=NA))
  m.ibound[i.bot,] <- t(apply(m.ibound[i.bot,], 1, num.replace, to.replace=c(n.coast.temp, abs(n.coast.temp), 180101080905), replace.with=NA))
  
  ## search for "trouble areas": active cells surrounded by inactive cells which will
  ## cause instability when recharge is applied because it has nowhere to go
  # pad ibound with 0s to deal with edge effects
  m.ibound.pad <- rbind(NA,cbind(NA, m.ibound,NA), NA)
  active <- which(is.finite(m.ibound.pad), arr.ind=T)
  
  neighbor.sum <- numeric(0)
  for (i in 1:dim(active)[1]){
    r <- active[i,'row']
    c <- active[i,'col']
    
    # get neighbors
    neighbors <- c(m.ibound.pad[r-1,c],
                   m.ibound.pad[r+1,c],
                   m.ibound.pad[r,c-1],
                   m.ibound.pad[r,c+1])
    
    neighbor.sum <- c(neighbor.sum, sum(is.finite(neighbors)))
    
    # status update
    if (i %% 100 == 0){
      print(paste0(i, " complete"))
    }
  }
  
  # get rid of isolated points
  bad.neighbors <- active[neighbor.sum==0, ]
  if (sum(neighbor.sum==0)>1){
    for (j in 1:dim(bad.neighbors)[1]){
      m.ibound.pad[bad.neighbors[j,'row'], bad.neighbors[j,'col']] <- NA
    }
  } else if (length(bad.neighbors==1)){
    m.ibound.pad[bad.neighbors[1], bad.neighbors[2]] <- NA
  }
  
  # trim padded row/col
  m.ibound <- m.ibound.pad[2:(dim(m.ibound.pad)[1]-1), 2:(dim(m.ibound.pad)[2]-1)]
  
  # put back into raster
  r.ibound[] <- m.ibound[]
  
  # trim white-space
  r.ibound <- trim(r.ibound, values=NA)
  m.ibound <- as.matrix(r.ibound)
  
  # extract HUC12 boundaries
  r.ibound.HUC12 <- r.ibound
  m.ibound.HUC12 <- as.matrix(r.ibound.HUC12)
  
  # set to 1s and 0s
  m.ibound[is.finite(m.ibound) & m.ibound != n.coast.temp] <- 1
  m.ibound[is.finite(m.ibound) & m.ibound == n.coast.temp] <- -1
  m.ibound[is.na(m.ibound)] <- 0
  
  # put back into raster
  r.ibound[] <- m.ibound[]
  
  # number of active cells
  paste0(sum(abs(m.ibound)), " active cells")
  
  # number of cells in navarro
  paste0(sum(substr(as.character(r.ibound.in[is.finite(r.ibound.in[])]), 1, 10)==HUC), " cells in Navarro")
  
  # save as geotiff and text
  writeRaster(r.ibound, "modflow/input/ibound.tif", datatype="INT2S", overwrite=T)
  writeRaster(r.ibound.HUC12, "modflow/input/ibound_HUC12.tif", overwrite=T)
  
  write.table(m.ibound, "modflow/input/ibound.txt", sep=" ", row.names=F, col.names=F)
  write.table(m.ibound.HUC12, "modflow/input/ibound_HUC12.txt", sep=" ", row.names=F, col.names=F)
  
  # extract data for ggplot
  df <- data.frame(rasterToPoints(r.ibound))
  colnames(df) <- c("lon", "lat", "ibound")
}

# Elevation (ZTOP) --------------------------------------------------------

if (ztop){
  # load elevation raster
  r.dem <- raster(paste0(dir.dem.NED, "Navarro_NED10m.vrt"))
  
  # reproject
  r.dem.proj <- projectRaster(r.dem, r.ibound)
  
  # set nodata
  r.dem.proj[r.dem.proj < 0] <- 0
  
  # extract as matrix
  m.dem.proj <- as.matrix(r.dem.proj)
  m.dem.proj[m.ibound==-1] <- 0
  m.dem.proj[m.ibound==0] <- 0
  
  # put back into raster
  r.dem.proj[] <- m.dem.proj[]
  
  # save as text file
  write.table(m.dem.proj, "modflow/input/ztop.txt", sep=" ", row.names=F, col.names=F)
  
  # save as geotiff
  writeRaster(r.dem.proj, "modflow/input/ztop.tif", overwrite=T)
  
} else {
  # load already processed data
  r.dem.proj <- raster("modflow/input/ztop.tif")
  m.dem.proj <- as.matrix(r.dem.proj)
}

# extract data
df$dem.m <- r.dem.proj[]

# RIV and SFR2 input data ---------------------------------------------------------

if (riv){
  # python needs a dictionary of boundaries in the format:
  #  [lay, row, col, stage, cond, rbot]
  
  # get rid of coastline (StreamOrde = -9)
  shp.riv.adj.streams <- subset(shp.riv.adj, StreamOrde>0)
  
  # eliminate tiny streams - some useful options might be:
  #   StreamOrde = stream order
  #   TotDASqKM = total upstream drainage area [km2]
  # plot(shp.riv.adj.streams)
  # plot(subset(shp.riv.adj.streams, StreamOrde >= 2))
  # plot(subset(shp.riv.adj.streams, TotDASqKM >= 5))
  shp.riv.adj.streams <- subset(shp.riv.adj.streams, StreamOrde >= 2)
  
  # reproject to UTM
  shp.riv.adj.streams.UTM <- spTransform(shp.riv.adj.streams, crs.MODFLOW)
  
  # keep only useful columns
  #   slope units are [m/m]
  shp.streams <- subset(shp.riv.adj.streams.UTM, 
                        select=c("OBJECTID", "REACHCODE", "TerminalPa", 
                                 "TotDASqKM", "StreamOrde", "TerminalFl",
                                 "SLOPE", "FromNode", "ToNode"))
  
  shp.streams <- shp.streams[which(!duplicated(shp.streams@data$OBJECTID)), ]
  
  
  # make stream segment number column
  shp.streams@data$SegNum <- seq(1, dim(shp.streams@data)[1], 1)
  
  # rasterize
  r.riv <- rasterize(shp.streams, r.ibound, field="StreamOrde", fun='max')
  m.riv.order <- as.matrix(r.riv)
  m.riv.order[m.ibound==0] <- NA
  
  ## calculate length of river in a cell
  r.riv.id <- r.riv
  r.riv.id[is.finite(r.riv.id[])] <- 1:sum(is.finite(r.riv.id[]))
  shp.riv.id <- rasterToPolygons(r.riv.id)   # make polygons corresponding to each cell with a river in it
  riv.int <- intersect(shp.streams, shp.riv.id)  # extract info from SpatialLines for each polygon - the new 'layer' field is the ID corresponding to the number in r.riv.id
  riv.int$length_m <- gLength(riv.int, byid=TRUE)
  riv.int@data <- 
    data.frame(riv.id = r.riv.id[],
               elev_m = r.dem.proj[],
               ibound = r.ibound[]) %>% 
    left_join(x=riv.int@data, y=., by=c("layer"="riv.id"))
  riv.int <- subset(riv.int, ibound != 0)  # get rid of inactive cells
    
  # sum length for each cell
  x <- tapply(riv.int$length, riv.int$layer, sum)
  r.riv.length <- r.riv.id
  r.riv.length[is.finite(r.riv.length[])] <- x
  
  # extract as matrix
  m.riv <- as.matrix(r.riv.length)
  
  # get rid of those outside of domain
  m.riv[m.ibound==0] <- NA
  
  # set anything with a length < 100 m to NA
  plot(r.riv.length>100)
  m.riv[m.riv<100] <- NA
  
  # put back into raster
  r.riv.length[] <- m.riv[]
  
  # rasterize and extract HUC12 watershed number
  r.HUC <- rasterize(shp.adj.UTM, r.ibound, field='HUC12', fun='max')
  m.HUC <- as.matrix(r.HUC)
  m.HUC[m.ibound==0] <- NA
  
  # matrix indices
  i.riv <- data.frame(which(is.finite(m.riv), arr.ind=T))
  
  # add layer
  i.riv$lay <- 1
  
  # extract elevation, to use as stage
  i.riv$stage_m <- m.dem.proj[which(is.finite(m.riv))]
  
  # extract river length, which is part of conductance calculation
  i.riv$length_m <- m.riv[which(is.finite(m.riv))]
  
  # extract stream order
  i.riv$stream_order <- m.riv.order[which(is.finite(m.riv))]
  
  # add a reach number
  i.riv$reach <- seq(1,dim(i.riv)[1])
  
  # extract HUC number
  i.riv$HUC <- m.HUC[which(is.finite(m.riv))]
  
  # because lay/row/col are all 0-based in Python, subtract 1
  i.riv$row <- i.riv$row-1
  i.riv$col <- i.riv$col-1
  i.riv$lay <- i.riv$lay-1
  
  # save as text file
  write.table(i.riv, "modflow/input/iriv.txt", sep=" ", row.names=F, col.names=T)
  
  # save as geotiff
  writeRaster(r.riv.length, "modflow/input/iriv.tif", datatype="INT2S", overwrite=T)
  
  # save as shapefile
  writeOGR(shp.streams, "modflow/input", "iriv", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  # put in data frame
  df$iriv <- r.riv.length[]
  
  ### now: SFR2 data
  # get all unique terminal points
  term.pts <- unique(riv.int$TerminalPa)
  
  # scroll through each outlet to ocean
  seg.counter <- 0
  for (term in term.pts){
    # get all points with this outlet
    riv.int.term <- as.data.frame(subset(riv.int, TerminalPa==term))
    
    # summarize number of reaches per segment
    riv.int.term.summary <- 
      group_by(riv.int.term, SegNum, TotDASqKM) %>% 
      summarize(n.reach = sum(is.finite(elev_m)))
    
    # order by drainage area (small-->large = upstream-->downstream)
    riv.int.term.summary <- riv.int.term.summary[order(riv.int.term.summary$TotDASqKM), ]
    
    # make SFR segment number
    riv.int.term.summary$SFR_NSEG <- seg.counter + seq(1,dim(riv.int.term.summary)[1])
    
    # add summary data
    riv.int.term <- left_join(riv.int.term, riv.int.term.summary[,c("SegNum", "SFR_NSEG")], by=c("SegNum"))
    
    # order by SFR_NSEG and dem (high-->low = upstream-->downstream)
    riv.int.term <- riv.int.term[with(riv.int.term, order(SFR_NSEG, -elev_m)), ]
    
    # number by reach within each segment
    riv.int.term <- 
      riv.int.term %>% 
      group_by(SFR_NSEG) %>% 
      mutate(SFR_IREACH = row_number())
    
    # add to overall output data frame
    if (seg.counter==0){
      riv.seg.all <- riv.int.term
    } else {
      riv.seg.all <- rbind(riv.seg.all, riv.int.term)
    }
    
    # update seg.counter
    seg.counter <- max(riv.int.term.summary$SFR_NSEG)

  }
  
  # get row/col/layer of each cell
  df.riv.id <- data.frame(ncell = seq(1,length(r.riv.id[])))
  df.riv.id[,c("row", "col")] <- rowColFromCell(r.riv.id, df.riv.id$ncell)
  df.riv.id$id <- extract(r.riv.id, df.riv.id$ncell)
  df.riv.id$ibound <- extract(r.ibound, df.riv.id$ncell)
  df.riv.id <- df.riv.id[is.finite(df.riv.id$id) & df.riv.id$ibound != 0, ]
  
  # join data frame
  riv.seg.all <- left_join(df.riv.id[,c("id", "row", "col")], riv.seg.all, by=c("id"="layer"))
  
  ## testing: only reaches that are terminal
  riv.seg.all <- subset(riv.seg.all, SFR_NSEG %in% unique(riv.seg.all$SFR_NSEG[riv.seg.all$TerminalFl==1]))
  riv.seg.all <- riv.seg.all[order(riv.seg.all$SFR_NSEG), ]
  riv.seg.all$SFR_NSEG <- as.numeric(factor(riv.seg.all$SFR_NSEG))
  
  # for each segment, determine OUTSEG
  riv.seg.outseg <- unique(riv.seg.all[,c("SegNum", "FromNode", "ToNode", "SFR_NSEG", "TerminalFl", "TerminalPa")])  # get all unique segments

  # use ToNode/FromNode to map segment connections
  riv.seg.outseg$SFR_OUTSEG <- 
    riv.seg.outseg$SFR_NSEG[match(riv.seg.outseg$ToNode, riv.seg.outseg$FromNode)]
  
  # anything that has TerminalFl==1 (end of flowline, e.g. at ocean) set OUTSEG to 0
  riv.seg.outseg$SFR_OUTSEG[riv.seg.outseg$TerminalFl==1] <- 0
  
  # anything that has SFR_OUTSEG terminates at the edge of our model domain but not at the ocean
  # (the HUC12 basins surrounding Navarro)
  riv.seg.outseg$SFR_OUTSEG[is.na(riv.seg.outseg$SFR_OUTSEG)] <- 0
  
  # make sure all OUTSEG are in NSEG
  sum(!(riv.seg.outseg$SFR_OUTSEG[riv.seg.outseg$SFR_OUTSEG != 0] %in% riv.seg.outseg$SFR_NSEG))
  
  # python indexing is 0-based so subtract 1 from row/col
  riv.seg.all$row <- riv.seg.all$row-1
  riv.seg.all$col <- riv.seg.all$col-1
  
  # save as text file
  write.table(riv.seg.all, "modflow/input/isfr_ReachData.txt", sep=" ", row.names=F, col.names=T)
  write.table(riv.seg.outseg, "modflow/input/isfr_SegmentData.txt", sep=" ", row.names=F, col.names=T)
  
}

# Prep pumping well data --------------------------------------------------

if (wel){
  # m.riv has locations of river cells
  # m.HUC has the HUC watershed locations
  m.navarro <- matrix(as.numeric(substr(as.character(m.HUC), 1, 10)==as.character(HUC)), nrow=dim(m.HUC)[1], ncol=dim(m.HUC)[2])
  
  # disable cells with rivers
  m.navarro[is.finite(m.riv)] <- 0
  
  # create grid of wells
  wel.spacing <- 1000  # [m]
  wel.spacing.row <- round(wel.spacing/DELC)
  wel.spacing.col <- round(wel.spacing/DELR)
  i.wel.row <- seq(1, dim(m.navarro)[1], wel.spacing.row)
  i.wel.col <- seq(1, dim(m.navarro)[1], wel.spacing.col)
  i.wel <- expand.grid(i.wel.row, i.wel.col)
  colnames(i.wel) <- c("row", "col")
  i.wel$lay <- 0  # 0 because FloPy is 0-based
  
  # figure out which are within Navarro watershed and aren't rivers
  i.wel$navarro <- m.navarro[cbind(i.wel[,1], i.wel[,2])]
  i.wel <- subset(i.wel, navarro==1)
  
  # get lon/lat from row/col indices
  i.wel$lon <- xFromCol(r.ibound, col=i.wel$col)
  i.wel$lat <- yFromRow(r.ibound, row=i.wel$row)
  
  # make a well number column
  i.wel$WellNum <- seq(1, dim(i.wel)[1], 1)
  
  # grab elevation
  i.wel$ztop_m <- m.dem.proj[cbind(i.wel[,1], i.wel[,2])]
  
  # because lay/row/col are all 0-based in Python, subtract 1
  i.wel$row <- i.wel$row-1
  i.wel$col <- i.wel$col-1
  
  # get rid of Navarro column
  i.wel$navarro <- NULL

  # save as text file
  write.table(i.wel, "modflow/input/iwel.txt", sep=" ", row.names=F, col.names=T)
  
}

# Make plots --------------------------------------------------------------

## prep polygon boundaries
df.basin <- tidy(spTransform(readOGR(dsn="data/NHD/WBD", layer="WBDHU10_Navarro"), crs.MODFLOW))
df.riv <- tidy(shp.riv.adj.streams.UTM)

## boundary conditions: IBOUND and RIV
# prep data
df$BCs <- "Inactive"
df$BCs[df$ibound==1] <- "Active"
df$BCs[df$ibound==-1] <- "Constant Head"
df$BCs[is.finite(df$iriv)] <- "River"
df$BCs <- factor(df$BCs, levels=c("Active", "Constant Head", "River", "Inactive"), 
                 labels=c("Active", "Constant\nHead", "River", "Inactive"))

# plot
p.BCs <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=BCs)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_point(data=i.wel, aes(x=lon, y=lat), shape=46) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_manual(name="B.C.", breaks=c("Active", "Constant\nHead", "River"),
                    values=c("Inactive"="white", "Active"="gray65", "Constant\nHead"="cyan", "River"="blue")) +
  coord_equal() +
  theme_scz() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

## elevation
df$dem.m[df$ibound==0] <- NaN
p.dem <-
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=dem.m)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_viridis(name="Elevation [m]", na.value="white") +
  coord_equal() +
  theme_scz() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5))

## save plots
# align
p1 <- ggplotGrob(p.BCs+theme(legend.position="bottom"))
p2 <- ggplotGrob(p.dem+theme(legend.position="bottom"))
p <- cbind(p1, p2, size="first")
p$heights <- unit.pmax(p1$heights, p2$heights)

# save
ggsave("modflow/input/ibound+iriv+ztop.png",
       p, width=190, height=100, units="mm")
