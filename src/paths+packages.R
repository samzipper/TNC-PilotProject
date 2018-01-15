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
require(viridis)

## metadata about watershed
station.outlet <- "11476500"  # USGS gage station number for outlet gauge
station.outlet.name <- "SF EEL R NR MIRANDA CA"
HUC <- "18010106"
area.mi2 <- 303  # gage contributing area (from https://waterdata.usgs.gov/nwis/inventory/?site_no=11468000)
area.km2 <- area.mi2*1.609344*1.609344

## info about all USGS gauging stations in watershed
# list of stations from https://waterdata.usgs.gov/nwis/inventory?huc2_cd=18010106&format=station_list&sort_key=site_no&group_key=county&list_of_search_criteria=huc2_cd
# this was passed to stationInfo in package waterData to get station info, and supplemented with area info from USGS NWIS website
df.info <- structure(list(staid = c("11476500", "11476600", "11475560", "11475610", "11475800"), 
                          staname = c("SF EEL R NR MIRANDA CA", "BULL C NR WEOTT CA", "ELDER C NR BRANSCOMB CA", 
                                      "CAHTO C NR LAYTONVILLE CA", "SF EEL R A LEGGETT CA"), 
                          lat = c(40.18181004, 40.351528, 39.7295997, 39.6715454, 39.8745957), 
                          lng = c(-123.7761426, -124.0050423, -123.6439073, -123.4983497, -123.7205775)
                          ), 
                     .Names = c("staid", "staname", "lat", "lng"), 
                     row.names = c(NA, -5L), 
                     class = "data.frame")
df.info$area.km2 <- c(303, 28.1, 6.50, 5.09, 248)*1.609344*1.609344   # area from website, converted from mi2 to km2

## labels
labs.mo <- c("1"="Jan", "2"="Feb", "3"="Mar", "4"="Apr", "5"="May", "6"="Jun",
             "7"="Jul", "8"="Aug", "9"="Sep", "10"="Oct", "11"="Nov", "12"="Dec")

## CRS for WGS plots
crs.WGS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## path to save cropped geospatial datasets
dir.gis <- "results/GIS/"

## paths to external datasets (e.g. global geospatial datasets on GSAS server)
# dem from NED
path.dem.NED <- "Z:/2.active_projects/Zipper/1.Spatial_data/regional/SouthForkEel/dem_topography/1original/NED10mDEM/NED10mDEM.vrt"

# dem, slope, river network from HydroSHEDS
path.dem.hysheds <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/dem_topography/1original/HydroSHEDS/DEM_15s/na_dem_15s.tif"
path.slope.hysheds <- "Z:/2.active_projects/Zipper/1.Spatial_data/global/dem_topography/1original/HydroSHEDS/slope_15s/na_slope_15s.tif"
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
