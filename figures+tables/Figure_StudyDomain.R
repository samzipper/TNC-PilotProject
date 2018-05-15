## Figure_StudyDomain.R
#' Map of study site with MODFLOW BCs and head at end of transient spin-up.

source(file.path("src", "paths+packages.R"))

# set up script
DELR <- 100  # should match script MODFLOW_Navarro_InputPrepData.R
DELC <- 100

## load data
# domain boundary shapefile
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv", stringsAsFactors=F)
df.streams <- tidy(shp.streams)

df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

# modflow inputs
r.ibound <- raster(file.path("modflow", "input", "ibound.tif"))

r.elev <- r.ibound
r.elev[] <- as.matrix(read.table(file.path("modflow", "input", "ztop.txt")))[]

df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), header=T)
df.riv <- read.table(file.path("modflow", "input", "iriv_ReachData.txt"), header=T)
df.sfr <- read.table(file.path("modflow", "input", "isfr_ReachData.txt"), header=T)

r.stream <- r.ibound
m.stream <- matrix(NaN, nrow=dim(r.stream)[1], ncol=dim(r.stream)[2])
m.stream[as.matrix(df.riv[,c("row", "col")])+1] <- 1
m.stream[as.matrix(df.sfr[,c("row", "col")])+1] <- 2
r.stream[] <- m.stream[]

r.wel <- r.ibound
m.wel <- matrix(NaN, nrow=dim(r.wel)[1], ncol=dim(r.wel)[2])
m.wel[as.matrix(df.wel[,c("row", "col")])+1] <- 1
r.wel[] <- m.wel[]

## data conversions
df <- 
  r.ibound %>% 
  rasterToPoints() %>% 
  data.frame() %>% 
  set_colnames(c("lon", "lat", "ibound"))
df$elev_m <- r.elev[]
df$stream <- r.stream[]
df$well <- r.wel[]

# clean up and make new columns
df$elev_m[df$ibound==0] <- NaN
df$ibound[df$ibound==0] <- NaN

df$BC <- NaN
df$BC[df$ibound==-1] <- "CHB"
df$BC[df$ibound==1] <- "Active"
df$BC[is.finite(df$stream)] <- "Stream"
df$BC <- factor(df$BC, levels=c("Stream", "CHB", "Active"))

## make plot
p.BC <- 
  ggplot() +
  geom_raster(data=subset(df, is.finite(BC)), aes(x=lon, y=lat, fill=BC), na.rm=T) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color=col.cat.red) +
  geom_point(data=df.wel, aes(x=lon, y=lat, color=is.finite(WellNum))) +
  scale_fill_manual(name="MODFLOW B.C.", 
                    values=c("Stream"=col.cat.blu, "CHB"=col.cat.org, "Active"=col.gray),
                    labels=c("Stream\n(RIV or SFR2)", "Constant\nHead", "Active")) +
  scale_color_manual(name=NULL, values=c("black"), labels=c("Synthetic\nWell")) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="bottom") + 
  guides(color=guide_legend(order=2),
         fill=guide_legend(order=1))

ggplot() +
  geom_raster(data=subset(df, is.finite(ibound)), aes(x=lon, y=lat, fill=elev_m), na.rm=T) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color=col.cat.red) +
  geom_path(data=shp.streams, aes(x=long, y=lat, group=group), color=col.cat.blu) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="bottom")
