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
DELR <- 100  # x spacing - column width
DELC <- 100  # y spacing - row height
DELV <- 20   # z spacing - layer height - only matters if nlay>1
nlay <- 5    # number of layers

# what variables to process?
ibound <- F
GLHYMPS <- F
ztop <- F
riv <- T
wel <- F

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
  
  # shapefile of active cells only
  r.ibound.out <- r.ibound
  r.ibound.out[r.ibound.out==0] <- NaN
  shp.ibound <- rasterToPolygons(r.ibound.out, na.rm=TRUE, digits=12, dissolve=FALSE)
  writeOGR(shp.ibound, "modflow/input", "ibound", driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  # save as geotiff and text
  writeRaster(r.ibound, file.path("modflow", "input", "ibound.tif"), datatype="INT2S", overwrite=T)
  writeRaster(r.ibound.HUC12, file.path("modflow", "input", "ibound_HUC12.tif"), overwrite=T)
  
  write.table(m.ibound, file.path("modflow", "input", "ibound.txt"), sep=" ", row.names=F, col.names=F)
  write.table(m.ibound.HUC12, file.path("modflow", "input", "ibound_HUC12.txt"), sep=" ", row.names=F, col.names=F)
  
} else {
  # load already processed data
  r.ibound <- raster(file.path("modflow", "input", "ibound.tif"))
  r.ibound.HUC12 <- raster(file.path("modflow", "input", "ibound_HUC12.tif"))
  
  m.ibound <- as.matrix(r.ibound)
  m.ibound.HUC12 <- as.matrix(r.ibound.HUC12)
}

# extract data for ggplot
df <- data.frame(rasterToPoints(r.ibound))
colnames(df) <- c("lon", "lat", "ibound")

# GLHYMPS data (k, n) -----------------------------------------------------

if (GLHYMPS){
  # read in rasters
  r.logk <- raster(file.path("results", "GIS", "Navarro_logK_GLHYMPS.tif"))
  r.porosity <- raster(file.path("results", "GIS", "Navarro_porosity_GLHYMPS.tif"))
  
  # first: just get some reasonable mean values to use for initial steady-state model
  logk.mean <- mean(r.logk[], na.rm=T)
  porosity.mean <- mean(r.porosity[], na.rm=T)
  
  # for Ksat, convert to m/day
  Ksat <- (10^logk.mean)*1e7*86400
  
  # for Sy, use 50% of porosity
  Sy <- porosity.mean*0.5
}


# Elevation (ZTOP) --------------------------------------------------------

if (ztop){
  # load elevation raster
  r.dem <- raster(file.path(dir.dem.NED, "Navarro_NED10m.vrt"))
  
  # reproject
  r.dem.proj <- projectRaster(r.dem, r.ibound)
  
  # set nodata
  r.dem.proj[r.dem.proj < 0] <- 0
  
  # extract as matrix
  m.dem.proj <- as.matrix(r.dem.proj)
  m.dem.proj[m.ibound==-1] <- 0
  m.dem.proj[m.ibound==0] <- 0
  
  # round off some digits
  m.dem.proj <- round(m.dem.proj, 2)
  
  # put back into raster
  r.dem.proj[] <- m.dem.proj[]
  
  # make additional layers
  if (nlay > 1){
    m.zbot.proj <- array(data=NA, dim=c(nrow(m.dem.proj), ncol(m.dem.proj), nlay))
    for (l in 1:nlay){
      if (l==1){
        m.zbot.proj[,,l] <- m.dem.proj - DELV
      } else {
        m.zbot.proj[,,l] <- m.zbot.proj[,,(l-1)] - DELV
      }
      write.table(m.zbot.proj[,,l], file.path("modflow", "input", paste0("zbot_", l, ".txt")), sep=" ", row.names=F, col.names=F)
    }
  }
  
  # save as text file
  write.table(m.dem.proj, file.path("modflow", "input", "ztop.txt"), sep=" ", row.names=F, col.names=F)
  
  # save as geotiff
  writeRaster(r.dem.proj, file.path("modflow", "input", "ztop.tif"), overwrite=T)
  
} else {
  # load already processed data
  r.dem.proj <- raster(file.path("modflow", "input", "ztop.tif"))
  m.dem.proj <- as.matrix(r.dem.proj)
}

# RIV and SFR2 input data ---------------------------------------------------------

if (riv){
  # # python needs a dictionary of boundaries in the format:
  # #  [lay, row, col, stage, cond, rbot]
  # 
  # ### only need to run this section after making changes to stream network (e.g. including different stream orders)
  # # get rid of coastline (StreamOrde = -9)
  # shp.riv.adj.streams <- subset(shp.riv.adj, StreamOrde>0)
  # 
  # # eliminate tiny streams - some useful options might be:
  # #   StreamOrde = stream order
  # #   TotDASqKM = total upstream drainage area [km2]
  # # plot(shp.riv.adj.streams)
  # # plot(subset(shp.riv.adj.streams, StreamOrde >= 2))
  # # plot(subset(shp.riv.adj.streams, TotDASqKM >= 5))
  # shp.riv.adj.streams <- subset(shp.riv.adj.streams, StreamOrde >= 2)
  # shp.riv.adj.streams@data$lineLength_m <- lengthLine(shp.riv.adj.streams)
  # 
  # # reproject to UTM
  # shp.riv.adj.streams.UTM <- spTransform(shp.riv.adj.streams, crs.MODFLOW)
  # 
  # # keep only useful columns
  # #   slope units are [m/m]
  # shp.streams <- subset(shp.riv.adj.streams.UTM,
  #                       select=c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m",
  #                                "TotDASqKM", "StreamOrde", "TerminalFl", "SLOPE", "FromNode", "ToNode"))
  # 
  # # extract duplicates
  # dups <- data.frame(i = c(which(duplicated(shp.streams@data$OBJECTID)),
  #                          which(duplicated(shp.streams@data$OBJECTID, fromLast=T)))) %>%
  #   cbind(., shp.streams@data[c(which(duplicated(shp.streams@data$OBJECTID)),
  #                               which(duplicated(shp.streams@data$OBJECTID, fromLast=T))),])
  # dups.keep <-
  #   dups %>%
  #   group_by(., OBJECTID) %>%
  #   summarize(lineLength_m = max(lineLength_m))
  # dups.null <- dups[!(paste0(dups$OBJECTID, "_", dups$lineLength_m) %in% paste0(dups.keep$OBJECTID, "_", dups.keep$lineLength_m)), ]
  # shp.streams <- shp.streams[-dups.null$i, ]
  # 
  # # make stream segment number column
  # shp.streams@data$SegNum <- seq(1, dim(shp.streams@data)[1], 1)
  # writeOGR(shp.streams, file.path("modflow", "input"), "iriv", driver="ESRI Shapefile", overwrite_layer=TRUE)
  shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")
  
  # shapefile shortens names; rename them
  names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                          "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")
  
  # rasterize
  r.riv <- rasterize(shp.streams, r.ibound, field="StreamOrde", fun='max')
  m.riv.order <- as.matrix(r.riv)
  m.riv.order[m.ibound==0] <- NA
  
  ## calculate length of river in a cell
  r.riv.id <- r.riv
  r.riv.id[is.finite(r.riv.id[])] <- 1:sum(is.finite(r.riv.id[]))
  df.riv.pts <- data.frame(rasterToPoints(r.riv.id)) %>% set_colnames(c("lon", "lat", "layer"))
  df.riv.pts$ncell <- cellFromXY(r.riv.id, xy=as.matrix(df.riv.pts[,c("lon", "lat")]))
  df.riv.pts[,c("row", "col")] <- rowColFromCell(r.riv.id, df.riv.pts$ncell)
  
  # ## this code block is VERY SLOW (hours to run)
  # ## only run it if you changed your grid resolution
  # shp.riv.id <- rasterToPolygons(r.riv.id)   # make polygons corresponding to each cell with a river in it
  # shp.riv.id.latlon <- spTransform(shp.riv.id, crs(r.dem))
  # df.riv.id.dem <- extract(r.dem, shp.riv.id.latlon, na.rm=T, df=T)
  # df.riv.id.dem.summary <-
  #   df.riv.id.dem %>% 
  #   group_by(ID) %>% 
  #   summarize(elev_m_min = min(Navarro_NED10m, na.rm=T),
  #             elev_m_med = median(Navarro_NED10m, na.rm=T),
  #             elev_m_mean = mean(Navarro_NED10m, na.rm=T),
  #             elev_m_max = max(Navarro_NED10m, na.rm=T))
  # df.riv.id.dem.summary$elev_m_min[df.riv.id.dem.summary$elev_m_min<0] <- 0
  # shp.riv.id@data <- left_join(shp.riv.id@data, df.riv.id.dem.summary, by=c("layer"="ID"))
  # writeOGR(shp.riv.id, "results/GIS", "MODFLOW_Navarro_InputPrepData_rivIDelev", driver="ESRI Shapefile", overwrite_layer=TRUE)
  shp.riv.id <- readOGR(file.path("results", "GIS"), "MODFLOW_Navarro_InputPrepData_rivIDelev")
  names(shp.riv.id) <- c("layer", "elev_m_min", "elev_m_med", "elev_m_mean", "elev_m_max")
  
  # get intersecting rivers
  riv.int <- intersect(shp.streams, shp.riv.id)  # extract info from SpatialLines for each polygon - the new 'layer' field is the ID corresponding to the number in r.riv.id
  riv.int$length_m <- gLength(riv.int, byid=TRUE)
  riv.int@data <- 
    data.frame(riv.id = r.riv.id[],
               elev_m = r.dem.proj[],
               ibound = r.ibound[]) %>% 
    left_join(x=riv.int@data, y=., by=c("layer"="riv.id")) %>% 
    left_join(x=., y=df.riv.pts, by=c("layer"))
  colnames(riv.int@data)[colnames(riv.int@data)=="layer"] <- "id"
  
  # summarize some metrics for each cell
  riv.int.cell <- 
    riv.int@data[,c("ncell", "id", "row", "col", "ibound", "elev_m_min", "elev_m", "TotDASqKM", "length_m", "SegNum")] %>% 
    group_by(ncell, id, row, col, ibound, elev_m_min, elev_m) %>% 
    summarize(totalLength_m = sum(length_m),
              drainageArea_km2 = max(TotDASqKM),
              n.segments = length(unique(SegNum))) %>% 
    subset(ibound != 0)

  # add n.segments and totalLength to riv.int
  riv.int@data <- left_join(riv.int@data, riv.int.cell)
  
  # calculate proportion of total length for a given cell that each segment has
  riv.int@data$seg_proportion <- riv.int@data$length_m/riv.int@data$totalLength_m
  
  # rasterize and extract HUC12 watershed number
  r.HUC <- rasterize(shp.adj.UTM, r.ibound, field='HUC12', fun='max')
  m.HUC <- as.matrix(r.HUC)  # need this for placing wells
  riv.int.cell$HUC <- raster::extract(r.HUC, riv.int.cell$ncell) 
  
  # matrix indices
  i.riv <- riv.int.cell[,c("row", "col", "elev_m_min", "elev_m", "drainageArea_km2", "totalLength_m")]
  
  # estimate river width based on drainage area
  i.riv$width_m <- WidthFromDA(DA=i.riv$drainageArea_km2, w.min=1, w.max=100)
  
  # add layer
  i.riv$lay <- 1

  # add a reach number
  i.riv$reach <- seq(1,dim(i.riv)[1])
  
  # because lay/row/col are all 0-based in Python, subtract 1
  i.riv$row <- i.riv$row-1
  i.riv$col <- i.riv$col-1
  i.riv$lay <- i.riv$lay-1
  
  ### now: SFR2 data
  # get all unique terminal points
  term.pts <- unique(riv.int$TerminalPa)
  
  # scroll through each outlet to ocean to define segment numbers
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
    
    # order by SFR_NSEG
    riv.int.term <- riv.int.term[with(riv.int.term, order(SFR_NSEG)), ]
    
    # add to overall output data frame
    if (seg.counter==0){
      riv.seg.all <- riv.int.term
    } else {
      riv.seg.all <- rbind(riv.seg.all, riv.int.term)
    }
    
    # update seg.counter
    seg.counter <- max(riv.int.term.summary$SFR_NSEG)
    
  }  # end of terminal paths loop
  
  ## for each river segment, figure out OUTSEG
  # get all unique segments
  riv.seg.info <- unique(riv.seg.all[,c("SegNum", "FromNode", "ToNode", "SFR_NSEG", "TotDASqKM", "TerminalFl", "TerminalPa")])
  
  # use ToNode/FromNode to map segment connections
  riv.seg.info$SFR_OUTSEG <- 
    riv.seg.info$SFR_NSEG[match(riv.seg.info$ToNode, riv.seg.info$FromNode)]
  
  # anything that has TerminalFl==1 (end of flowline, e.g. at ocean) or 
  # ends at the edge of our model domain, set OUTSEG to 0
  riv.seg.info$SFR_OUTSEG[riv.seg.info$TerminalFl==1] <- 0
  riv.seg.info$SFR_OUTSEG[is.na(riv.seg.info$SFR_OUTSEG)] <- 0
  
  ## now, figure out reach numbers...
  # for each segment, order reaches by upstream-->downstream
  nseg.all <- unique(riv.seg.all$SFR_NSEG)
  for (nseg in nseg.all){
    # subset data to only this segment
    riv.seg <- subset(riv.seg.all, SFR_NSEG==nseg)[,!(names(riv.seg.all) %in% c("totalLength_m", "n.segments", "seg_proportion"))]
    riv.seg[,c("lon", "lat", "ncell")] <- NULL
    
    # figure out which reach has water flowing out of it
    NSEG.outseg <- riv.seg.info$SFR_OUTSEG[riv.seg.info$SFR_NSEG==nseg]
    riv.seg.outseg <- subset(riv.seg.all, SFR_NSEG %in% NSEG.outseg)
    riv.seg.outseg <- riv.seg.outseg[which.max(riv.seg.outseg$elev_m_min),]
    if (dim(riv.seg.outseg)[1]>0){
      # figure out difference in row/col units (OK because delr=delc)
      tot.diff.outseg <- 
        sqrt(abs(outer(riv.seg$row, riv.seg.outseg$row, FUN="-"))^2 +
               abs(outer(riv.seg$col, riv.seg.outseg$col, FUN="-"))^2)
      
      # determine index in the segment of interest which is closest to the outseg
      outseg <- riv.seg$id[unique(which(tot.diff.outseg == min(tot.diff.outseg), arr.ind = TRUE)[,1])]
    } else {
      # if nothing flows out of this one, choose the cell with ibound==-1 (coastal cell)
      outseg <- riv.seg$id[which(riv.seg$ibound==-1)]
      
      # if there is no coastal cell, choose lowest elevation
      if (length(outseg)==0){
        outseg <- riv.seg$id[which.min(riv.seg$elev_m_min)]
      }
    }
    
    # if two are equally close, choose lower elevation
    if (length(outseg)>1){
      outseg <- riv.seg$id[riv.seg$elev_m_min==min(riv.seg$elev_m_min[riv.seg$id %in% outseg]) & riv.seg$id %in% outseg]
    }
    
    # if there are still two that have equal elevation, just choose 1
    if (length(outseg)>1){
      riv.seg <- subset(riv.seg, !(id %in% outseg[2:length(outseg)]))
      outseg <- outseg[1]
    }
    
    # get rid of any cells with an elevation lower than outseg
    riv.seg <- subset(riv.seg, elev_m_min >= riv.seg$elev_m_min[riv.seg$id==outseg])
    
    # figure out which reach has water flowing into it
    NSEG.inseg <- riv.seg.info$SFR_NSEG[riv.seg.info$SFR_OUTSEG==nseg]
    riv.seg.inseg <- subset(riv.seg.all, SFR_NSEG %in% NSEG.inseg)
    riv.seg.inseg <- riv.seg.inseg[which.min(riv.seg.inseg$elev_m_min),]
    if (dim(riv.seg.inseg)[1]>0){
      tot.diff.inseg <- 
        sqrt(abs(outer(riv.seg$row, riv.seg.inseg$row, FUN="-"))^2 +
               abs(outer(riv.seg$col, riv.seg.inseg$col, FUN="-"))^2)
      
      inseg <- riv.seg$id[unique(which(tot.diff.inseg == min(tot.diff.inseg), arr.ind = TRUE)[,1])]
    } else {
      # if there are no segments flowing into this one, choose the highest elevation as the start
      inseg <- riv.seg$id[which.max(riv.seg$elev_m_min)]
    }
    
    # if two are equally close, choose higher elevation for inseg
    if (length(inseg)>1){
      inseg <- riv.seg$id[riv.seg$elev_m_min==max(riv.seg$elev_m_min[riv.seg$id %in% inseg]) & riv.seg$id %in% inseg]
    }
    
    # if there are still two that have equal elevation, just choose 1
    if (length(inseg)>1){
      riv.seg <- subset(riv.seg, !(id %in% inseg[2:length(inseg)]))
      inseg <- inseg[1]
    }
    
    if (inseg==outseg){
      riv.seg <- subset(riv.seg, id==inseg)
    }
    
    # see how many reaches there are after eliminations based on inseg/outseg and only continue if more than 1
    if (dim(riv.seg)[1]>1){
      
      # find and eliminate local sinks (no neighboring cells with lower elevation)
      n.lower <- function(elev, row, col, elev_all, row_all, col_all){
        elev_ngb <- elev_all[(row_all %in% seq(row-1, row+1)) & 
                               (col_all %in% seq(col-1, col+1))]
        sum(elev_ngb<elev)
      }
      
      riv.seg$ngb.lower <- NaN
      for (i in 1:dim(riv.seg)[1]){
        riv.seg$ngb.lower[i] <- 
          n.lower(elev=riv.seg$elev_m_min[i], row=riv.seg$row[i], col=riv.seg$col[i], 
                  elev_all=riv.seg$elev_m_min, row_all=riv.seg$row, col_all=riv.seg$col)
      }
      
      riv.seg <- subset(riv.seg, (ngb.lower > 0) | (id == outseg) | (id==inseg))
      riv.seg$ngb.lower <- NULL
      
      ## if there are a bunch of reaches, move downgradient and make a network
      riv.seg$SFR_IREACH <- NaN
      path.complete <- F
      riv.seg$SFR_IREACH[riv.seg$id==inseg] <- 1
      while (!path.complete){
        # find maximum existing ireach (current terminus of path)
        term.ireach <- max(riv.seg$SFR_IREACH, na.rm=T)
        i.term <- which.max(riv.seg$SFR_IREACH)
        
        # get row/col of terminus
        term.row <- riv.seg$row[i.term]
        term.col <- riv.seg$col[i.term]
        term.elev <- riv.seg$elev_m_min[i.term]
        
        # find any neighboring (incl. diagonal) cells of lower elevation
        riv.seg.ngb <- subset(riv.seg, 
                              (row %in% seq((term.row-1), (term.row+1))) & 
                                (col %in% seq((term.col-1), (term.col+1))) &
                                elev_m_min < term.elev)
        
        if (dim(riv.seg.ngb)[1]>0){
          # if there is a neighboring cell, select highest elevation
          i.next <- which(riv.seg$id==riv.seg.ngb$id[which.max(riv.seg.ngb$elev_m_min)])
          
          # update terminal ireach
          riv.seg$SFR_IREACH[i.next] <- term.ireach+1
          
        } else {
          # if there is no neighboring cell, find closest cell with lower elevation
          dist.cells <- 
            sqrt(
              abs(term.row-riv.seg$row[riv.seg$elev_m_min < term.elev])^2 + 
                abs(term.col-riv.seg$col[riv.seg$elev_m_min < term.elev])^2
            )
          
          # get closest cell
          id.closest <- riv.seg$id[riv.seg$elev_m_min < term.elev][which(dist.cells==min(dist.cells))]
          
          # if >1 have same distance, get higher elevation
          id.closest <- id.closest[which.max(riv.seg$elev_m_min[riv.seg$id %in% id.closest])]
          i.closest <- which(riv.seg$id==id.closest)
          
          # get row/col of closest cell
          closest.row <- riv.seg$row[i.closest]
          closest.col <- riv.seg$col[i.closest]
          
          if (abs(closest.row - term.row)>1 & abs(closest.col - term.col)>1){
            # get indices of row/col in between
            path.row <- seq(term.row, closest.row)[-c(1, (length(seq(term.row, closest.row))))]
            path.col <- seq(term.col, closest.col)[-c(1, (length(seq(term.col, closest.col))))]
            
            # find number of diagonals
            n.diag <- min(c(length(path.row), length(path.col)))
            
            df.path <- data.frame(row=path.row[1:n.diag],
                                  col=path.col[1:n.diag],
                                  length_m = sqrt(DELR^2+DELC^2))
            
            # if there are still rows/col left
            if (length(path.row)>n.diag){
              df.path.strt <- data.frame(row = path.row[(n.diag+1):length(path.row)],
                                         col = closest.col,
                                         length_m = DELC)
              df.path <- rbind(df.path, df.path.strt)
              
            } else if (length(path.col)>n.diag){
              df.path.strt <- data.frame(row = closest.row,
                                         col = path.col[(n.diag+1):length(path.col)],
                                         length_m = DELR)
              df.path <- rbind(df.path, df.path.strt)
            }
            
          } else if (abs(closest.row - term.row)>1){
            path.row <- seq(term.row, closest.row)
            path.row <- path.row[!(path.row %in% c(term.row, closest.row))]
            df.path <- data.frame(row = path.row,
                                  col = term.col,
                                  length_m = DELC)
          } else if (abs(closest.col - term.col)>1){
            path.col <- seq(term.col, closest.col)
            path.col <- path.col[!(path.col %in% c(term.col, closest.col))]
            df.path <- data.frame(row = term.row,
                                  col = path.col,
                                  length_m = DELR)
          } else {
            stop("error in creating path between non-neighboring cells")
          }
          
          # fill in other necessary columns
          df.path$id <- -9999
          df.path$OBJECTID <- -9999
          df.path$SFR_IREACH <- seq(from=term.ireach+1, by=1, length.out=dim(df.path)[1])
          df.path$elev_m_min <- seq(from=riv.seg$elev_m_min[i.term], to=riv.seg$elev_m_min[i.closest],
                                    length.out=dim(df.path)[1]+2)[-c(1, (dim(df.path)[1]+2))]
          df.path$elev_m_mean <- seq(from=riv.seg$elev_m_mean[i.term], to=riv.seg$elev_m_mean[i.closest],
                                     length.out=dim(df.path)[1]+2)[-c(1, (dim(df.path)[1]+2))]
          df.path$elev_m_med <- seq(from=riv.seg$elev_m_med[i.term], to=riv.seg$elev_m_med[i.closest],
                                    length.out=dim(df.path)[1]+2)[-c(1, (dim(df.path)[1]+2))]
          df.path$elev_m_max <- seq(from=riv.seg$elev_m_max[i.term], to=riv.seg$elev_m_max[i.closest],
                                    length.out=dim(df.path)[1]+2)[-c(1, (dim(df.path)[1]+2))]
          df.path$elev_m <- seq(from=riv.seg$elev_m[i.term], to=riv.seg$elev_m[i.closest],
                                length.out=dim(df.path)[1]+2)[-c(1, (dim(df.path)[1]+2))]
          df.path[,c("REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", "TerminalFl", "SLOPE",
                     "FromNode", "ToNode", "SegNum", "ibound", "drainageArea_km2", "SFR_NSEG")] <-
            riv.seg[i.term, c("REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", "TerminalFl", "SLOPE",
                              "FromNode", "ToNode", "SegNum", "ibound", "drainageArea_km2", "SFR_NSEG")]
          
          # update IREACH for i.closest
          riv.seg$SFR_IREACH[i.closest] <- max(df.path$SFR_IREACH)+1
          
          # add df.path to riv.seg
          riv.seg <- rbind(riv.seg, df.path)
          rm(df.path)
          
        }
        
        # figure out if you've reached outseg yet
        if (riv.seg$id[which.max(riv.seg$SFR_IREACH)]==outseg){
          path.complete <- T
        }
        
        if (max(riv.seg$SFR_IREACH, na.rm=T)>dim(riv.seg)[1]){
          error("too many paths")
        }
        
        # status update
        print(paste0("nseg ", nseg, "; ireach ", max(riv.seg$SFR_IREACH, na.rm=T), " assigned, path.complete=", path.complete))
        
      }
      
      # retain only reaches that are connected via the stream network
      riv.seg <- subset(riv.seg, is.finite(SFR_IREACH))
      
    } else {
      # if there's only 1 reach in the segment, just set it to 1
      riv.seg$SFR_IREACH <- 1
    }
    
    # combine into output data frame
    if (nseg==nseg.all[1]){
      riv.seg.out <- riv.seg[order(riv.seg$SFR_IREACH), ]
    } else {
      riv.seg.out <- rbind(riv.seg.out, riv.seg[order(riv.seg$SFR_IREACH), ])
    }
    
  }
  
  # get rid of any elev < 0
  riv.seg.out$elev_m_min[riv.seg.out$elev_m_min<0] <- 0
  
  ## now: make sure each segment ends at or next to the start of the OUTSEG it drains into;
  ##      repeat this until no more elevation adjustment necessary
  n.toolow <- 1
  while (n.toolow > 0){
    # initialize n.toolow
    n.toolow <- 0
    
    # get start and end of each segment
    riv.seg.starts <- subset(riv.seg.out, SFR_IREACH==1)
    riv.seg.ends <- 
      riv.seg.out[,c("SFR_NSEG", "SFR_IREACH")] %>% 
      group_by(SFR_NSEG) %>% 
      summarize(SFR_IREACH = max(SFR_IREACH),
                end.of.seg = T) %>% 
      left_join(riv.seg.out, ., by=c("SFR_NSEG", "SFR_IREACH")) %>% 
      subset(end.of.seg)
    riv.seg.ends$end.of.seg <- NULL
    
    # scroll through all nseg that are not terminal segments
    nseg.check <- subset(riv.seg.info, SFR_OUTSEG != 0)$SFR_NSEG
    for (nseg in nseg.check){
      # get the row/col of the end of this segment
      row.this <- riv.seg.ends$row[riv.seg.ends$SFR_NSEG==nseg]
      col.this <- riv.seg.ends$col[riv.seg.ends$SFR_NSEG==nseg]
      elev.this <- riv.seg.ends$elev_m_min[riv.seg.ends$SFR_NSEG==nseg]
      
      # get OUTSEG for this segment
      outseg.this <- riv.seg.info$SFR_OUTSEG[riv.seg.info$SFR_NSEG==nseg]
      
      # get max reach # for this segment
      ireach.max.this <- riv.seg.ends$SFR_IREACH[riv.seg.ends$SFR_NSEG==nseg]
      
      # get the row/col for the start of the next segment
      row.next <- riv.seg.starts$row[riv.seg.starts$SFR_NSEG==outseg.this]
      col.next <- riv.seg.starts$col[riv.seg.starts$SFR_NSEG==outseg.this]
      elev.next <- riv.seg.starts$elev_m_min[riv.seg.starts$SFR_NSEG==outseg.this]
      
      # fix any cells that are too low
      if (elev.this < elev.next){
        i.toolow <- which(riv.seg.out$elev_m_min < elev.next & riv.seg.out$SFR_NSEG==nseg)
        riv.seg.out$elev_m_min[i.toolow] <- elev.next
        riv.seg.ends$elev_m_min[riv.seg.ends$SFR_NSEG==nseg] <- elev.next
        
        n.toolow <- n.toolow + length(i.toolow)
        print(paste0("Too low: ", length(i.toolow), " of ", sum(riv.seg.out$SFR_NSEG==nseg), " reaches on ", nseg))
      }
      
      # check if they are in neighboring cells (diagonal OK)
      if (abs(row.this-row.next)>1 | abs(col.this-col.next)>1){
        if (abs(row.this-row.next)>1 & abs(col.this-col.next)>1){
          
          # get indices of row/col in between
          path.row <- seq(row.this, row.next)[-c(1, (length(seq(row.this, row.next))))]
          path.col <- seq(col.this, col.next)[-c(1, (length(seq(col.this, col.next))))]
          
          # find number of diagonals
          n.diag <- min(c(length(path.row), length(path.col)))
          
          df.conn <- data.frame(row=path.row[1:n.diag],
                                col=path.col[1:n.diag],
                                length_m = sqrt(DELR^2+DELC^2),
                                SFR_NSEG = nseg)
          
          # if there are still rows/col left
          if (length(path.row)>n.diag){
            df.conn.strt <- data.frame(row = path.row[(n.diag+1):length(path.row)],
                                       col = col.next,
                                       length_m = DELC,
                                       SFR_NSEG = nseg)
            df.conn <- rbind(df.conn, df.conn.strt)
            
          } else if (length(path.col)>n.diag){
            df.conn.strt <- data.frame(row = row.next,
                                       col = path.col[(n.diag+1):length(path.col)],
                                       length_m = DELR,
                                       SFR_NSEG = nseg)
            df.conn <- rbind(df.conn, df.conn.strt)
          }
          
        } else if (abs(row.next - row.this)>1){
          path.row <- seq(row.this, row.next)
          path.row <- path.row[!(path.row %in% c(row.this, row.next))]
          df.conn <- data.frame(row = path.row,
                                col = col.this,
                                length_m = DELC,
                                SFR_NSEG = nseg)
        } else if (abs(col.next - col.this)>1){
          path.col <- seq(col.this, col.next)
          path.col <- path.col[!(path.col %in% c(col.this, col.next))]
          df.conn <- data.frame(row = row.this,
                                col = path.col,
                                length_m = DELR,
                                SFR_NSEG = nseg)
        }
        
        # fill in other necessary columns
        df.conn$id <- -9998
        df.conn$OBJECTID <- -9998
        df.conn$SFR_IREACH <- seq(from=ireach.max.this+1, by=1, length.out=dim(df.conn)[1])
        df.conn$elev_m_min <- seq(from=riv.seg.ends$elev_m_min[riv.seg.ends$SFR_NSEG==nseg], 
                                  to=riv.seg.starts$elev_m_min[riv.seg.starts$SFR_NSEG==outseg.this],
                                  length.out=dim(df.conn)[1]+2)[-c(1, (dim(df.conn)[1]+2))]
        df.conn$elev_m_mean <- seq(from=riv.seg.ends$elev_m_mean[riv.seg.ends$SFR_NSEG==nseg], 
                                   to=riv.seg.starts$elev_m_mean[riv.seg.starts$SFR_NSEG==outseg.this],
                                   length.out=dim(df.conn)[1]+2)[-c(1, (dim(df.conn)[1]+2))]
        df.conn$elev_m_med <- seq(from=riv.seg.ends$elev_m_med[riv.seg.ends$SFR_NSEG==nseg], 
                                  to=riv.seg.starts$elev_m_med[riv.seg.starts$SFR_NSEG==outseg.this],
                                  length.out=dim(df.conn)[1]+2)[-c(1, (dim(df.conn)[1]+2))]
        df.conn$elev_m_max <- seq(from=riv.seg.ends$elev_m_max[riv.seg.ends$SFR_NSEG==nseg], 
                                  to=riv.seg.starts$elev_m_max[riv.seg.starts$SFR_NSEG==outseg.this],
                                  length.out=dim(df.conn)[1]+2)[-c(1, (dim(df.conn)[1]+2))]
        df.conn$elev_m <- seq(from=riv.seg.ends$elev_m[riv.seg.ends$SFR_NSEG==nseg], 
                              to=riv.seg.starts$elev_m[riv.seg.starts$SFR_NSEG==outseg.this],
                              length.out=dim(df.conn)[1]+2)[-c(1, (dim(df.conn)[1]+2))]
        df.conn[,c("REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", "TerminalFl", "SLOPE",
                   "FromNode", "ToNode", "drainageArea_km2", "SegNum", "ibound")] <-
          riv.seg.starts[riv.seg.ends$SFR_NSEG==nseg, c("REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", "TerminalFl", "SLOPE",
                                                        "FromNode", "ToNode", "drainageArea_km2", "SegNum", "ibound")]
        
        # status update
        print(paste0("Connecting ", nseg, " to ", outseg.this, ": ", dim(df.conn)[1], " reach(es) created"))
        
        # add to overall data frame
        riv.seg.out <- rbind(riv.seg.out, df.conn)
        rm(df.conn)
        
      }
    }
    
  }  # end of while loop
  
  ## figure out direct drainage area for each segment
  riv.seg.info$DirectDASqKM <- NaN
  for (i in 1:length(riv.seg.info$SFR_NSEG)){
    nseg <- riv.seg.info$SFR_NSEG[i]
    
    # figure out which NSEG are draining into this
    nseg_upstream <- riv.seg.info$SFR_NSEG[riv.seg.info$SFR_OUTSEG==nseg]
    
    # calculate directDA
    riv.seg.info$DirectDASqKM[i] <- riv.seg.info$TotDASqKM[i] - 
      sum(riv.seg.info$TotDASqKM[riv.seg.info$SFR_NSEG %in% nseg_upstream])
    
    # some SegNums have a weird situation with two parallel channels draining into them and need special treatment
    if (riv.seg.info$SegNum[i] %in% c(71, 79, 134, 307)){
      riv.seg.info$DirectDASqKM[i] <- riv.seg.info$TotDASqKM[i] - 
        max(riv.seg.info$TotDASqKM[riv.seg.info$SFR_NSEG %in% nseg_upstream])
    }
    
    # other SegNums are weird for other reasons so just set direct drainage to 0
    if (riv.seg.info$SegNum[i] %in% c(72, 80, 136, 148, 354, 358, 368)){
      # 72, 80, 136 = parallel to other segment
      # 148, 354, 358, 368 = terminal segments outside Navarro
      riv.seg.info$DirectDASqKM[i] <- 0
    }
    
  }
  
  ## check slope
  # define minimum allowed slope [m/m]
  slope.min <- 0.001
  riv.seg.out$SLOPE[riv.seg.out$SLOPE<slope.min] <- slope.min
  
  # estimate river width based on drainage area
  riv.seg.out$width_m <- WidthFromDA(DA=riv.seg.out$TotDASqKM, w.min=1, w.max=100)
  riv.seg.info$width_m <- WidthFromDA(DA=riv.seg.info$TotDASqKM, w.min=1, w.max=100)
  
  ## extract coordinates and add to data frame
  riv.seg.out <-
    cellFromRowCol(r.riv.id, rownr=riv.seg.out$row, colnr=riv.seg.out$col) %>% 
    xyFromCell(r.riv.id, cell=.) %>% 
    data.frame(.) %>% 
    set_colnames(c("lon", "lat")) %>% 
    cbind(riv.seg.out, .)
  
  # python indexing is 0-based so subtract 1 from row/col
  riv.seg.out$row <- riv.seg.out$row-1
  riv.seg.out$col <- riv.seg.out$col-1
  
  # put in order
  riv.seg.info <- riv.seg.info[order(riv.seg.info$SFR_NSEG), ]
  riv.seg.out <- riv.seg.out[with(riv.seg.out, order(SFR_NSEG, SFR_IREACH)), ]
  
  ## subtract 1 from row/col numbers in riv.int@data to match other data frames
  riv.int@data$row <- riv.int@data$row-1
  riv.int@data$col <- riv.int@data$col-1
  riv.int@data <- subset(riv.int@data, ibound != 0)
  
  ## set up gaging stations
  # define gages by USGS gage numbers
  gage.stations <- station.outlet
  n.gage.stations <- length(gage.stations)
  
  # get coordinates
  df.gages <- 
    siteInfo(gage.stations)[,c("lng", "lat")] %>% 
    SpatialPointsDataFrame(coords=., 
                           data=data.frame(GageNum = seq(1,n.gage.stations)),
                           proj4string=crs(shp.riv.adj)) %>% 
    spTransform(., crs(crs.MODFLOW))
  
  # add coordinates to data frame
  df.gages@data <- 
    cbind(df.gages@data, 
          data.frame(lon=df.gages@coords[,"lng"], 
                     lat=df.gages@coords[,"lat"]))
  
  # find closest stream reach to gage (coordinates will not be exactly the same
  # because stream reach coordinates are the center of each MODFLOW grid cell)
  tot.diff.gages <- 
    sqrt(abs(outer(df.gages$lon, riv.seg.out$lon, FUN="-"))^2 +
           abs(outer(df.gages$lat, riv.seg.out$lat, FUN="-"))^2)
  
  df.gages@data <- 
    which(tot.diff.gages == min(tot.diff.gages), arr.ind = TRUE)[,"col"] %>% 
    riv.seg.out[., c("row", "col", "REACHCODE", "TotDASqKM", "StreamOrde", "SegNum", "SFR_NSEG", "SFR_IREACH")] %>% 
    cbind(df.gages@data, .)
  
  ## update elevation output
  # there can be multiple reaches with different min elevations per cell; set DEM equal to the highest one
  df.rowcol.elev <- 
    riv.seg.out[,c("row", "col", "lon", "lat", "elev_m_min")] %>% 
    rbind(., riv.int@data[,c("row", "col", "lon", "lat", "elev_m_min")]) %>% 
    group_by(row, col) %>% 
    summarize(elev_m_min_max = max(elev_m_min),
              lon = mean(lon),
              lat = mean(lat))
  m.dem.proj[as.matrix(df.rowcol.elev[,c("row", "col")])+1] <- 
    df.rowcol.elev$elev_m_min_max
  m.dem.proj <- round(m.dem.proj, 2)
  r.dem.proj[] <- m.dem.proj[]
  
  if (nlay > 1){
    m.zbot.proj <- array(data=NA, dim=c(nrow(m.dem.proj), ncol(m.dem.proj), nlay))
    for (l in 1:nlay){
      if (l==1){
        m.zbot.proj[,,l] <- m.dem.proj - DELV
      } else {
        m.zbot.proj[,,l] <- m.zbot.proj[,,(l-1)] - DELV
      }
      write.table(m.zbot.proj[,,l], file.path("modflow", "input", paste0("zbot_", l, ".txt")), sep=" ", row.names=F, col.names=F)
    }
  }
  
  write.table(m.dem.proj, file.path("modflow", "input", "ztop.txt"), sep=" ", row.names=F, col.names=F)
  writeRaster(r.dem.proj, file.path("modflow", "input", "ztop.tif"), overwrite=T)
  
  ## save output data
  # RIV
  write.table(i.riv, file.path("modflow", "input", "iriv.txt"), sep=" ", row.names=F, col.names=T)
  write.table(riv.int@data, file.path("modflow", "input", "iriv_ReachData.txt"), sep=" ", row.names=F, col.names=T)

  # SFR
  write.table(riv.seg.out, file.path("modflow", "input", "isfr_ReachData.txt"), sep=" ", row.names=F, col.names=T)
  write.table(riv.seg.info, file.path("modflow", "input", "isfr_SegmentData.txt"), sep=" ", row.names=F, col.names=T)
  
  # GAGE
  write.table(df.gages, file.path("modflow", "input", "gage_data.txt"), sep=" ", row.names=F, col.names=T)
  
  # add SFR data to output
  df <- left_join(df, df.rowcol.elev[,c("lon", "lat", "elev_m_min_max")], by=c("lon", "lat"))
}

# Prep pumping well data --------------------------------------------------

if (wel){
  # m.riv has locations of river cells
  # m.HUC has the HUC watershed locations
  m.navarro <- matrix(as.numeric(substr(as.character(m.HUC), 1, 10)==as.character(HUC)), nrow=dim(m.HUC)[1], ncol=dim(m.HUC)[2])
  
  # disable cells with either RIV or SFR
  m.navarro[as.matrix(df.rowcol.elev[,c("row", "col")])+1] <- 0
  m.navarro[as.matrix(riv.int@data[,c("row", "col")])+1] <- 0
  
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
  
  df <- left_join(df, i.wel[,c("lon", "lat", "WellNum")], by=c("lon", "lat"))
}

# Make plots --------------------------------------------------------------

# extract elevation data
df$dem.m <- r.dem.proj[]

## prep polygon boundaries
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)
df.riv <- tidy(shp.streams)

## boundary conditions: IBOUND and RIV
# prep data
df$BCs <- "Inactive"
df$BCs[df$ibound==1] <- "Active"
df$BCs[df$ibound==-1] <- "Constant Head"
df$BCs[is.finite(df$SFR_NSEG)] <- "SFR"
df$BCs <- factor(df$BCs, levels=c("Active", "Constant Head", "SFR", "Inactive"), 
                 labels=c("Active", "Constant\nHead", "SFR", "Inactive"))

# plot
p.BCs <- 
  ggplot() +
  geom_raster(data=df, aes(x=lon, y=lat, fill=BCs)) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_point(data=i.wel, aes(x=lon, y=lat), shape=46) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_manual(name="B.C.", breaks=c("Active", "Constant\nHead", "SFR"),
                    values=c("Inactive"="white", "Active"="gray65", "Constant\nHead"="cyan", "SFR"="blue")) +
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
