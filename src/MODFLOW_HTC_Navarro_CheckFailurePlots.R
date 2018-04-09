## MODFLOW_HTC_Navarro_OutputPlots.R
#' This script is intended to plot output data from MODFLOW.
#' 
#' We are using EPSG:26910 as our projected CRS for MODFLOW, 
#' which has units of meters.
#' 
#' For the domain, we are using the Navarro River watershed
#' (HUC 1801010804) plus all adjacent HUC12 watersheds.

source("src/paths+packages.R")

## which stream BC are you using?
stream_BC = "RIV"  # RIV or SFR

## define which directory you are interested in
dir.runs <- file.path("modflow", "HTC", "Navarro", "SteadyState", stream_BC)

## load common data
# domain boundary shapefile
shp <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro")
shp.UTM <- spTransform(shp, crs.MODFLOW)

shp.adj <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

# streamlines
shp.streams <- readOGR(dsn="modflow/input", layer="iriv", stringsAsFactors=F)

## prep polygon boundaries for plots
df.basin <- tidy(shp.UTM)
df.basin.adj <- tidy(shp.adj.UTM)
df.riv <- tidy(shp.streams)

## load ibound raster
r.ibound <- raster("modflow/input/ibound.tif")
DELR <- res(r.ibound)[1]
DELC <- res(r.ibound)[2]

# make data frame
df <- data.frame(rasterToPoints(r.ibound))
colnames(df) <- c("lon", "lat", "ibound")

# Pumping well plots ------------------------------------------------------

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), header=T)

# Python indices are 0-based; add 1 for R
df.wel[,c("row", "col", "lay")] <- df.wel[,c("row", "col", "lay")]+1

# load model failure report and add to data frame
df.wel <- read.csv(file.path(dir.runs, "CheckFailure.csv"), stringsAsFactors=F) %>% 
  left_join(df.wel, ., by="WellNum")
df.wel$Success <- as.logical(df.wel$Success)

# well success
p.wel.succ <- 
  ggplot() +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_point(data=df.wel, aes(x=lon, y=lat, color=Success)) +
  labs(title=paste0(dir.runs, "; Success=", sum(df.wel$Success), ", Failure=", sum(df.wel$Success==F))) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_fill_viridis(name="Head [m]", na.value="white") +
  scale_color_manual(name="Converged?", values=c("FALSE"=col.cat.org, "TRUE"=col.cat.grn)) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position=c(0.01, 0.01),
        legend.justification=c(0,0))
ggsave(file.path(dir.runs, "CheckFailure.png"), p.wel.succ,
       width=150, height=150, units="mm")
