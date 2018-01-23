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
DELR <- 1000  # x spacing - column width
DELC <- 1000  # y spacing - row height

# Make boundary conditions (IBOUND) ---------------------------------------

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

# save as text file
write.table(m.ibound, "modflow/input/ibound.txt", sep=" ", row.names=F, col.names=F)
