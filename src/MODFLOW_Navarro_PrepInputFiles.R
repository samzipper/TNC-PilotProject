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
  shp.riv.adj <- readOGR(dsn="data/NHD/HYD", layer="NHDFlowline_HU12_Navarro+Adjacent")
  
  # reproject to UTM
  shp.riv.adj.UTM <- spTransform(shp.riv.adj, crs.MODFLOW)
  
  # rasterize
  r.riv <- rasterize(shp.riv.adj.UTM, r.ibound)
}
