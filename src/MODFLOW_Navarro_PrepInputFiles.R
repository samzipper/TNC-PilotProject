## Navarro_MODFLOW_PrepInputFiles.R
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
  
  # save as text file
  write.table(m.dem.proj, "modflow/input/ztop.txt", sep=" ", row.names=F, col.names=F)
} else {
  # load output
  m.dem.proj <- as.matrix(read.table("modflow/input/ztop.txt", header=F))
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
}
