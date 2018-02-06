## Navarro_MODFLOW_InputPrepData.R
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
ztop <- T
riv <- T

## load common data
# domain boundary shapefile
shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent")
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
  #  their numbers are:
  #   15 (north end - delete 0-1300 m)
  #   19 (south end - delete 37300 - end m)
  i.15 <- 1:round(1300/DELC)
  i.19 <- round(37300/DELC):dim(m.ibound)[1]
  num.replace <- function(x, to.replace, replace.with){
    x[x %in% to.replace] <- replace.with
    return(x)
  }
  
  m.ibound[i.15,] <- t(apply(m.ibound[i.15,], 1, num.replace, to.replace=c(n.coast.temp, abs(n.coast.temp), 15), replace.with=NA))
  m.ibound[i.19,] <- t(apply(m.ibound[i.19,], 1, num.replace, to.replace=c(n.coast.temp, abs(n.coast.temp), 19), replace.with=NA))
  
  # set to 1s and 0s
  m.ibound[is.finite(m.ibound) & m.ibound != n.coast.temp] <- 1
  m.ibound[is.finite(m.ibound) & m.ibound == n.coast.temp] <- -1
  m.ibound[is.na(m.ibound)] <- 0
  
  ## search for "trouble areas": active cells surrounded by inactive cells which will
  ## cause instability when recharge is applied because it has nowhere to go
  
  # pad ibound with 0s to deal with edge effects
  m.ibound.pad <- rbind(0,cbind(0, m.ibound,0), 0)
  active <- which(abs(m.ibound.pad)==1, arr.ind=T)
  
  neighbor.sum <- numeric(0)
  for (i in 1:dim(active)[1]){
    r <- active[i,'row']
    c <- active[i,'col']
    
    # get neighbors
    neighbors <- c(m.ibound.pad[r-1,c],
                   m.ibound.pad[r+1,c],
                   m.ibound.pad[r,c-1],
                   m.ibound.pad[r,c+1])
    
    neighbor.sum <- c(neighbor.sum, sum(abs(neighbors)))
    
    # status update
    if (i %% 100 == 0){
      print(paste0(i, " complete"))
    }
  }
  
  # get rid of isolated points
  bad.neighbors <- active[neighbor.sum==0, ]
  if (sum(neighbor.sum==0)>1){
    for (j in 1:dim(bad.neighbors)[1]){
      m.ibound.pad[bad.neighbors[j,'row'], bad.neighbors[j,'col']] <- 0
    }
  } else if (length(bad.neighbors==1)){
    m.ibound.pad[bad.neighbors[1], bad.neighbors[2]] <- 0
  }
  
  # trim padded row/col
  m.ibound <- m.ibound.pad[2:(dim(m.ibound.pad)[1]-1), 2:(dim(m.ibound.pad)[2]-1)]
  
  # put back into raster
  r.ibound[] <- m.ibound[]
  
  # trim white-space
  r.ibound <- trim(r.ibound, values=0)
  
  # re-extract boundary conditions to matrix
  m.ibound <- as.matrix(r.ibound)
  
  # number of active cells
  paste0(sum(is.finite(m.ibound)), " active cells")
  
  # save as text file
  write.table(m.ibound, "modflow/input/ibound.txt", sep=" ", row.names=F, col.names=F)
  
  # save as geotiff
  writeRaster(r.ibound, "modflow/input/ibound.tif", datatype="INT2S", overwrite=T)
  
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
  
  # extract data
  df$dem.m <- r.dem.proj[]
  
}

# River locations ---------------------------------------------------------

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
  
  # rasterize
  r.riv <- rasterize(shp.riv.adj.streams.UTM, r.ibound, field="StreamOrde", fun='max')
  
  # extract as matrix
  m.riv <- as.matrix(r.riv)
  
  # get rid of those outside of domain
  m.riv[m.ibound==0] <- NA
  
  # put back into raster
  r.riv[] <- m.riv[]
  
  # matrix indices
  i.riv <- data.frame(which(is.finite(m.riv), arr.ind=T))
  
  # add layer
  i.riv$lay <- 1
  
  # extract elevation, to use as stage
  i.riv$stage <- m.dem.proj[which(is.finite(m.riv))]
  
  # because lay/row/col are all 0-based in Python, subtract 1
  i.riv$row <- i.riv$row-1
  i.riv$col <- i.riv$col-1
  i.riv$lay <- i.riv$lay-1
  
  # save as text file
  write.table(i.riv, "modflow/input/iriv.txt", sep=" ", row.names=F, col.names=T)
  
  # save as geotiff
  writeRaster(r.riv, "modflow/input/iriv.tif", datatype="INT2S", overwrite=T)
  
  # put in data frame
  df$iriv <- r.riv[]
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
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_manual(name="B.C.", breaks=c("Active", "Constant\nHead", "River"),
                    values=c("Inactive"="white", "Active"="gray65", "Constant\nHead"="black", "River"="blue")) +
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
