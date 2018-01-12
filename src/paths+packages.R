## paths+packages.R
# This script holds packages, paths, color schemes, etc which are shared among scripts.

## packages
require(ggplot2)
require(dplyr)
require(reshape2)
require(lubridate)
require(grid)
require(gridExtra)
require(hydrostats)
require(waterData)
require(rgdal)
require(raster)
require(rgeos)
require(maptools)
require(broom)

## metadata about watershed
station <- "11468000"  # USGS gage station number
station.name <- "NAVARRO R NR NAVARRO CA"
area.mi2 <- 303  # gage contributing area (from https://waterdata.usgs.gov/nwis/inventory/?site_no=11468000)
area.km2 <- area.mi2*1.609344*1.609344

## labels
labs.mo <- c("1"="Jan", "2"="Feb", "3"="Mar", "4"="Apr", "5"="May", "6"="Jun",
             "7"="Jul", "8"="Aug", "9"="Sep", "10"="Oct", "11"="Nov", "12"="Dec")

## path to save cropped geospatial datasets
dir.gis <- "results/GIS/"

## paths to external datasets (e.g. global geospatial datasets on GSAS server)
# dem, slope, river network from HydroSHEDS
path.dem <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/dem_topography/1original/HydroSHEDS/DEM_15s/na_dem_15s.tif"
path.slope <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/dem_topography/1original/HydroSHEDS/slope_15s/na_slope_15s.tif"
dir.riv <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/riv_river_network_lakes/1original/HydroSHEDS/na_riv_15s"

# soil and bedrock depth from SoilGrids1km
dir.SoilGrids <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/s_soils/"

# porosity and permeability from GLHYMPS
dir.GLHYMPS <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/k_permeability_porosity/2derived/GLHYMPS/v1_2014/geotiffs/EPSG4326/"

# water table depth from Fan et al. (2013)
path.wtd <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/wt_water_table_wells/1original/wt1_fan2103_water_table_depth/geotiffs/v2/N_America_model_wtd_v2.tif"

## ggplot theme
windowsFonts(Arial=windowsFont("TT Arial"))
theme_scz <- function(...){
  theme_bw(base_size=10, base_family="Arial") + 
    theme(
      text=element_text(color="black"),
      axis.title=element_text(face="bold", size=rel(1)),
      axis.text=element_text(size=rel(1)),
      strip.text=element_text(size=rel(1)),
      legend.title=element_text(face="bold", size=rel(1)),
      legend.text=element_text(size=rel(1)),
      panel.grid=element_blank())
}

## functions
# crop and mask a raster based on a shapefile
crop.mask <- function(r, shp){
  return(mask(crop(r, extent(shp)), shp))
}

# extract p-value from linear model fits
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
