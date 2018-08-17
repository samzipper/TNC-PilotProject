## paths+packages.R
# This script holds packages, paths, color schemes, etc which are shared among scripts.

## packages
require(ggtern)
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
require(stringr)
require(gtable)
require(EcoHydRology)
require(zoo)
require(gstat)
require(geosphere)
require(magrittr)
require(dismo)
require(hydroGOF)
require(RCurl)
require(SPEI)
require(tidyr)

## metadata about watershed
station.outlet <- "11468000"  # USGS gage station number for outlet gauge
station.outlet.name <- "NAVARRO R NR NAVARRO CA"
station.outlet.lat <- 39.17055556
station.outlet.lon <- -123.66694444
HUC <- "1801010804"  # navarro
outlet.TerminalPa <- 10013609
area.mi2 <- 303  # gage contributing area (from https://waterdata.usgs.gov/nwis/inventory/?site_no=11468000)
area.km2 <- area.mi2*1.609344*1.609344
baseflow.cfs <- 8.4  # baseflow requirement from Table 7 in final cannabis policy https://www.waterboards.ca.gov/board_decisions/adopted_orders/resolutions/2017/final_cannabis_policy_with_att_a.pdf
baseflow.mm_d <- baseflow.cfs*(0.3048^3)*(0.001^3)*(1/area.km2)*86400*1000*1000

# ## info about all USGS gauging stations in watershed - this is for south fork eel
# # list of stations from https://waterdata.usgs.gov/nwis/inventory?huc2_cd=18010106&format=station_list&sort_key=site_no&group_key=county&list_of_search_criteria=huc2_cd
# # this was passed to stationInfo in package waterData to get station info, and supplemented with area info from USGS NWIS website
# df.info <- structure(list(staid = c("11476500", "11476600", "11475560", "11475610", "11475800"), 
#                           staname = c("SF EEL R NR MIRANDA CA", "BULL C NR WEOTT CA", "ELDER C NR BRANSCOMB CA", 
#                                       "CAHTO C NR LAYTONVILLE CA", "SF EEL R A LEGGETT CA"), 
#                           lat = c(40.18181004, 40.351528, 39.7295997, 39.6715454, 39.8745957), 
#                           lng = c(-123.7761426, -124.0050423, -123.6439073, -123.4983497, -123.7205775)
#                           ), 
#                      .Names = c("staid", "staname", "lat", "lng"), 
#                      row.names = c(NA, -5L), 
#                      class = "data.frame")
# df.info$area.km2 <- c(303, 28.1, 6.50, 5.09, 248)*1.609344*1.609344   # area from website, converted from mi2 to km2
# df.info$baseflow.cfs <- c(54, 1.9, 1.3, 2.4, 25)   # baseflow requirement from Table 4 in final cannabis policy https://www.waterboards.ca.gov/board_decisions/adopted_orders/resolutions/2017/final_cannabis_policy_with_att_a.pdf
# df.info$baseflow.mm_d <- df.info$baseflow.cfs*(0.3048^3)*(0.001^3)*(1/df.info$area.km2)*86400*1000*1000

## labels
labs.mo <- c("1"="Jan", "2"="Feb", "3"="Mar", "4"="Apr", "5"="May", "6"="Jun",
             "7"="Jul", "8"="Aug", "9"="Sep", "10"="Oct", "11"="Nov", "12"="Dec")

labs.analytical <- c("glover"="Glover", "hunt"="Hunt")
labs.method <- c("Qf.InvDist"="Inverse\nDistance", "Qf.InvDistSq"="Inverse\nDistance\nSquared",
                 "Qf.Web"="Web", "Qf.WebSq"="Web\nSquared", "Qf.TPoly"="Thiessen\nPolygon")
labs.stream_BC <- c("RIV"="RIV", "SFR"="SFR")

## CRS for WGS plots
crs.WGS <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

## CRS for MODFLOW
crs.MODFLOW <- "+init=epsg:26910"   # NAD83 UTM Zone 10 N

## path to save cropped geospatial datasets
dir.gis <- "results/GIS/"

## paths to external datasets (e.g. global geospatial datasets on GSAS server)
# directory with CA-DWR water table depth data
dir.gw.dwr <- "Z:/2.active_projects/Zipper/1.Spatial_data/regional/NavarroRiver/wt_water_table_wells/1original/Statewide_GWL_Data_20170905/"

# directory with diversions data
dir.div <- "Z:/2.active_projects/Zipper/1.Spatial_data/regional/SouthForkEel/use_withdrawl_abstraction_use/1original/CaliforniaWaterRights"

# dem from NED
dir.dem.NED <- "Z:/2.active_projects/Zipper/1.Spatial_data/regional/NavarroRiver/dem_topography/1original/NED10mDEM/"

# NLCD: impervious, canopy, land cover
dir.nlcd <- "Z:/2.active_projects/Zipper/1.Spatial_data/regional/NavarroRiver/lulc_imagery/1original/"

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
      plot.title=element_text(face="bold", size=rel(1)),
      axis.title=element_text(face="bold", size=rel(1)),
      axis.text=element_text(size=rel(1)),
      strip.text=element_text(size=rel(1)),
      legend.title=element_text(face="bold", size=rel(1)),
      legend.text=element_text(size=rel(1)),
      panel.grid=element_blank(),
      plot.margin=unit(c(1,1,1,1), "mm"),
      strip.background=element_blank())
}

theme_set(theme_scz())

## axis breaks
map.breaks.x <- seq(440000, 480000, 20000)
map.breaks.y <- seq(4300000, 4340000, 20000)

## color palettes
# categorical color palette from https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
col.cat.grn <- "#3cb44b"   # green
col.cat.yel <- "#ffe119"   # yellow
col.cat.org <- "#f58231"   # orange
col.cat.red <- "#e6194b"   # red
col.cat.blu <- "#0082c8"   # blue
col.gray <- "gray65"       # gray for annotation lines, etc
pal.method <- c("f.TPoly"=col.cat.grn, "f.InvDist"=col.cat.yel, "f.InvDistSq"=col.cat.org, "f.Web"=col.cat.red, "f.WebSq"=col.cat.blu)
pal.method.Qf <- c("Qf.TPoly"=col.cat.grn, "Qf.InvDist"=col.cat.yel, "Qf.InvDistSq"=col.cat.org, "Qf.Web"=col.cat.red, "Qf.WebSq"=col.cat.blu, "Qf.NoApport"="black")
labels.method <- c("f.TPoly"="Thiessen", "f.InvDist"="Inverse", "f.InvDistSq"="Inverse\nSquared", "f.Web"="Web", "f.WebSq"="Web\nSquared")
labels.method.Qf <- c("Qf.TPoly"="Thiessen", "Qf.InvDist"="Inverse", "Qf.InvDistSq"="Inverse\nSquared", "Qf.Web"="Web", "Qf.WebSq"="Web\nSquared", "Qf.NoApport"="No\nApportionment")
labels.apportionment <- c("WholeDomain"="Whole Domain", "AdjacentOnly"="Adjacent", "LocalArea"="Local Area", "Dynamic"="Dynamic (10 yrs)")

# NLCD color palette
pal.NLCD <- c("11"="#5475A8", 
              "12"="#ffffff", 
              "21"="#E8D1D1", 
              "22"="#E29E8C", 
              "23"="#ff0000", 
              "24"="#B50000", 
              "31"="#D2CDC0", 
              "41"="#85C77E", 
              "42"="#38814E", 
              "43"="#D4E7B0", 
              "51"="#AF963C", 
              "52"="#DCCA8F", 
              "71"="#FDE9AA", 
              "72"="#D1D182", 
              "73"="#A3CC51", 
              "74"="#82BA9E", 
              "81"="#FBF65D", 
              "82"="#CA9146", 
              "90"="#C8E6F8", 
              "95"="#64B3D5")

labels.NLCD <- c("11"="Open Water", 
                 "12"="Perennial Ice", 
                 "21"="Developed-Open", 
                 "22"="Developed-Low", 
                 "23"="Developed-Med", 
                 "24"="Developed-High", 
                 "31"="Barren", 
                 "41"="Decid Forest", 
                 "42"="Everg Forest", 
                 "43"="Mixed Forest", 
                 "51"="Dwarf Scrub", 
                 "52"="Shrub/Scrub", 
                 "71"="Grassland", 
                 "72"="Sedge", 
                 "73"="Lichen", 
                 "74"="Moss", 
                 "81"="Pasture/Hay", 
                 "82"="Cropland", 
                 "90"="Woody Wetland", 
                 "95"="Herb. Wetland")

## functions
# estimate stream width based on drainage area
# (see Navarro_StreamWidthEstimates.xlsx file)
WidthFromDA <- function(DA, w.min, w.max){
  # DA = drainage area [km2]
  # w.min = minimum allowed width [m]
  # w.max = maximum allowed width [m]
  # w = estimated width [m]
  w <- 9.7133*exp(0.0023*DA)
  w[w>w.max] <- w.max
  w[w<w.min] <- w.min
  return(w)
}

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

# extract legend - https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

## fit functions
## functions from Gudmundsson et al. (2012) for modified version of KGE
# eq. 5 - units of output from these will be same as input units of sim and obs
#         the ideal value for each of these is 0.0
sd.p <- function(x){sqrt((length(x)-1)/length(x))*sd(x)}  # from https://stackoverflow.com/questions/44339070/calculating-population-standard-deviation-in-r
MSE.bias <- function(sim,obs){(mean(sim)-mean(obs))^2}
MSE.var <- function(sim,obs){(sd.p(sim)-sd.p(obs))^2}
MSE.cor <- function(sim,obs){
  if (sd(obs)==0 | sd(sim)==0){
    0 
  } else {
    2*sd.p(sim)*sd.p(obs)*(1-cor(sim,obs))
  }
}
MSE <- function(sim,obs){MSE.bias(sim,obs)+MSE.var(sim,obs)+MSE.cor(sim,obs)}   # this outputs slightly different results than mse() in the hydroGOF package

# eq. 6 - the ideal value for each is 0.0,but these are normalized and will always sum to 1.0
MSE.bias.norm <- function(sim,obs){
  MSE.bias(sim,obs)/MSE(sim,obs)
}

MSE.var.norm <- function(sim,obs){
  MSE.var(sim,obs)/MSE(sim,obs)
}

MSE.cor.norm <- function(sim,obs){
  MSE.cor(sim,obs)/MSE(sim,obs)
}

# R2
R2 <- function(sim, obs) {
  if (length(sim) != length(obs)) stop("vectors not the same size")
  return((sum((obs-mean(obs))*(sim-mean(sim)))/
            ((sum((obs-mean(obs))^2)^0.5)*(sum((sim-mean(sim))^2)^0.5)))^2)
}

## function for sourcing scripts directly from GitHub
source_github <- function(u) {
  # from: https://stackoverflow.com/questions/35720660/how-to-use-an-r-script-from-github
  # load package
  require(RCurl)
  
  # read script lines from website
  script <- getURL(u, ssl.verifypeer = FALSE)
  
  # parase lines and evaluate in the global environment
  eval(parse(text = script))
}
