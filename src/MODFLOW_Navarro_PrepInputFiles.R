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
DELR <- 500  # x spacing - column width
DELC <- 500  # y spacing - row height

# what variables to process?
ibound <- T
ztop <- F
riv <- T

# Make boundary conditions (IBOUND) ---------------------------------------

if (ibound){
  # load shapefile
  shp.adj <- readOGR(dsn="data/NHD/WBD", layer="WBDHU12_Navarro+Adjacent")
  
  # reproject to UTM
  shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)
  
  # make empty raster
  r.empty <- raster(extent(shp.adj.UTM), crs = crs.MODFLOW, resolution=c(DELR, DELC))
  
  # rasterize
  r.ibound <- trim(rasterize(shp.adj.UTM, r.empty, field="HUC12"))
  
  # extract as matrix
  m.ibound <- as.matrix(r.ibound)
  
  # set to 1s and 0s
  m.ibound[is.finite(m.ibound)] <- 1
  m.ibound[is.na(m.ibound)] <- 0
  
  # number of active cells
  paste0(sum(m.ibound), " active cells")
  
  # save as text file
  write.table(m.ibound, "modflow/input/ibound.txt", sep=" ", row.names=F, col.names=F)
}

# Elevation (ZTOP) --------------------------------------------------------

if (ztop){
  # load elevation raster
  r.dem <- raster(paste0(dir.dem.NED, "Navarro_NED10m.vrt"))
  
  # crop
  r.dem.crop <- crop(r.dem, extent(spTransform(shp.adj, crs(r.dem))))
  
  # reproject
  r.dem.proj <- projectRaster(r.dem, r.ibound)
  
  # set nodata
  r.dem.proj[r.dem.proj < 0] <- 0
  
  # extract as matrix
  m.dem.proj <- as.matrix(r.dem.proj)
  m.dem.proj <- m.dem.proj*m.ibound
  
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
  
  # load shapefile
  shp.riv.adj <- readOGR(dsn="data/NHDPlusV2/HYD", layer="NHDPlusV21_National_Seamless NHDFlowline_Network_WBDHU12_Navarro+Adjacent",
                         stringsAsFactors=F)
  
  # convert stream order to numeric
  shp.riv.adj@data$StreamOrde <- as.numeric(shp.riv.adj@data$StreamOrde)
  
  # get rid of coastline (StreamOrde = -9)
  shp.riv.adj <- subset(shp.riv.adj, StreamOrde>0)
  
  # eliminate tiny streams - some useful options might be:
  #   StreamOrde = stream order
  #   TotDASqKM = total upstream drainage area [km2]
  # plot(shp.riv.adj)
  # plot(subset(shp.riv.adj, StreamOrde >= 2))
  # plot(subset(shp.riv.adj, TotDASqKM >= 5))
  shp.riv.adj <- subset(shp.riv.adj, StreamOrde >= 2)
  
  # reproject to UTM
  shp.riv.adj.UTM <- spTransform(shp.riv.adj, crs.MODFLOW)
  
  # rasterize
  r.riv <- rasterize(shp.riv.adj.UTM, r.ibound, field="StreamOrde", fun='max')
  
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
