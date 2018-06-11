## Navarro_Analytical_Depletion_OneWellAllReaches.R

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
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

## make a plot
wel <- 527
df.apportionment %>% 
  subset(WellNum==wel) %>% 
  left_join(., df.streams, by=c("SegNum")) %>% 
  ggplot() +
  geom_path(aes(x=long, y=lat, group=group, 
                color=cut(f.WebSq, 
                          breaks=c(0, 0.05, 0.1, 0.15, 0.2, 1),
                          include.lowest=T)),
            size=1.5) +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_point(data=subset(df.wel, WellNum==wel), aes(x=lon, y=lat), color="black",
             size=2, shape=21) +
  scale_color_manual(name="Depletion from\n100 cfs pumping",
                     labels=c("0-5 cfs", "5-10 cfs", "10-15 cfs", "15-20 cfs", ">20 cfs"),
                     values=c("#313695", "#00D9D9", "#18A718", "#D06F1D", "#D01D1D"), 
                     drop=F) +
  scale_x_continuous(name="Long") +
  scale_y_continuous(name="Lat") +
  coord_equal() +
  theme(legend.position="bottom") +
  ggsave("results/Navarro_Analytical_Depletion_OneWellAllReaches.png",
         width=6.25, height=6, units="in")
