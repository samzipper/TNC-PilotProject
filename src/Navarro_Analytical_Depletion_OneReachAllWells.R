## Navarro_Analytical_Depletion_OneReachAllWells.R

source(file.path("src", "paths+packages.R"))

## load depletion apportionment results
df.apportionment <- 
  file.path("results","Navarro_DepletionApportionment_AllMethods+Wells+Reaches.csv") %>% 
  read.csv()

## load stream data - created in MODFLOW_Navarro_InputPrepData.R
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# data frame for ggplots
df.streams <- tidy(shp.streams, id=SegNum)
df.streams$SegNum <- as.numeric(df.streams$id) + 1

## load well data
df.wel <- 
  file.path("modflow", "input", "iwel.txt") %>% 
  read.table(header=T)

## load basin outline
shp.basin <-
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW)
df.basin <- 
  shp.basin %>% 
  tidy(.)

## make raster of empty basin
r.empty <- 
  file.path("data", "NHD", "WBD") %>% 
  readOGR(dsn=., layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  extent() %>% 
  raster(., crs=crs.MODFLOW)
res(r.empty) <- c(100,100)

## figure out which reach to plot
reach <- 287
df.apportionment %>% 
  subset(SegNum==reach) %>% 
  ggplot(aes(x=f.WebSq)) +
  geom_histogram()

## build inverse distance model
df.depletion <- left_join(df.wel, df.apportionment, by="WellNum")
sp.wel <- SpatialPoints(df.depletion[,c("lon", "lat")], proj4string=crs(crs.MODFLOW))
sp.wel <- SpatialPointsDataFrame(sp.wel, df.depletion)

gs <- gstat(formula=as.formula("f.WebSq~1"), locations=sp.wel)

# interpolate to raster
r.interp <- interpolate(r.empty, gs)

# mask with shapefile
r.interp.mask <- mask(r.interp, shp.basin)

