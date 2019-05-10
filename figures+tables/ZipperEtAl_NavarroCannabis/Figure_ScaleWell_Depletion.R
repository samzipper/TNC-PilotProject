## Figure_ScaleWellt_Depletion.R
#' This script is intended to relate segment-level and well-level depletion
#' Very similar to Figure_Habitat+GrowLocations.R

source(file.path("src", "paths+packages.R"))

## load data from Navarro_Cannabis_HabitatIntrinsicPotential.R
# stream shapefile (from NHD)
sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  subset(TermnlP == outlet.TerminalPa) %>% 
  left_join(read.csv(file.path("results", "Navarro_Cannabis_HabitatIntrinsicPotential.csv")),
            by=c("SegNum"))

# grow locations from TNC (filtered and transformed in Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.R)
sf.grows <-
  sf::st_read(file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows.gpkg"), stringsAsFactors=F) %>%
  sf::st_transform(crs.MODFLOW)

# rasters
r.dem.30m <- raster(paste0(dir.gis, "Navarro_Cannabis_DEM_30m.tif"))
r.dtb <- raster(paste0(dir.gis, "Navarro_Cannabis_DTB_30m.tif"))
r.wte <- raster(paste0(dir.gis, "Navarro_Cannabis_WTE_30m.tif"))

## calculate saturated thickness of alluvial materials are saturated or not
r.wtd <- r.dem.30m - r.wte
r.alluvial.sat.thickness <- r.dtb - r.wtd
r.alluvial.sat.thickness[r.alluvial.sat.thickness < r.dtb] <- 0
r.alluvial.sat.thickness[r.alluvial.sat.thickness > r.dtb] <- r.dtb[r.alluvial.sat.thickness > r.dtb]

# domain boundary shapefile
df.basin <- 
  readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro") %>% 
  spTransform(., crs.MODFLOW) %>% 
  tidy(.)

sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp")) %>% 
  sf::st_transform(crs.MODFLOW)

sf.subbasin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU12_Navarro.shp")) %>% 
  sf::st_transform(crs.MODFLOW)
sf.subbasin$area_m2 <- as.numeric(sf::st_area(sf.subbasin))

#### load and process cannabis data
# depletion by segment associated with each well - output from Navarro_Cannabis-Grows_03_DepletionBySegment.R
df <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_03_DepletionBySegment.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  subset(GrowNum %in% sf.grows$GrowNum[sf.grows$Well.rf.pred=="Yes"])

# calculate year and month
df$year.dec <- df$time_days/365 + 1
df$year <- floor(df$year.dec)
df$month <- round((df$year.dec - df$year)*12)+1

# load well-stream pairs which has distance
df.pairs <- 
  file.path(dir.TNC, "DerivedData", "Navarro_Cannabis-Grows_01_CalculateWellStreamPairs.csv") %>% 
  read.csv(stringsAsFactors=F)

# calculate total water use for each grow
sf.grows$WaterUseTotal_m3y <- NaN
for (i in 1:dim(sf.grows)[1]){
  sf.grows$WaterUseTotal_m3y[i] <- 
    sum(
      as.numeric(
        st_drop_geometry(sf.grows)[i,c("JanWU_Estimate", "FebWU_Estimate", "MarWU_Estimate", "AprWU_Estimate",
                                       "MayWU_Estimate", "JunWU_Estimate", "JulWU_Estimate", "AugWU_Estimate",
                                       "SepWU_Estimate", "OctWU_Estimate", "NovWU_Estimate", "DecWU_Estimate")]
      )*days_in_month(seq(1,12))
    )
}

# what years to plot? options: 1  2  3  4  5 10 15 20 30 50
yrs.plot <- c(1, 10, 50)
mo.plot <- 9

# summarize total depletion by year, month, and grow
df.depletion.total.YearMonth <- 
  df %>% 
  subset(year %in% yrs.plot & month %in% mo.plot) %>% 
  group_by(year, month) %>% 
  summarize(depletion_m3d_total = sum(depletion_m3d))

df.depletion.total.YearMonthGrow <-
  df %>% 
  subset(year %in% yrs.plot & month %in% mo.plot) %>% 
  group_by(year, month, GrowNum) %>% 
  summarize(depletion_m3d_grow = sum(depletion_m3d)) %>% 
  left_join(df.depletion.total.YearMonth, by=c("year", "month")) %>% 
  transform(depletion_prc_grow = depletion_m3d_grow/depletion_m3d_total)

df.depletion.rank <- 
  df.depletion.total.YearMonthGrow[order(df.depletion.total.YearMonthGrow$year, df.depletion.total.YearMonthGrow$depletion_prc_grow), ] %>% 
  # rank from smallest contribution to depletion (1) to largest contribution to depletion (n.gw) within each year
  transform(depletion_rank = c(seq(sum(df.depletion.total.YearMonthGrow$year==yrs.plot[1]), 1, -1),
                               seq(sum(df.depletion.total.YearMonthGrow$year==yrs.plot[2]), 1, -1),
                               seq(sum(df.depletion.total.YearMonthGrow$year==yrs.plot[3]), 1, -1))) %>% 
  # calculate rank and depletion normalized 0-1
  transform(depletion_rank_prc = depletion_rank/length(unique(df.depletion.total.YearMonthGrow$GrowNum)),
            depletion_prc_cum = ave(depletion_prc_grow, year, FUN=cumsum))

# for each year, figure out what rank corresponds to 50% depletion
df.lines <- 
  df.depletion.rank %>% 
  group_by(year) %>% 
  filter(abs(depletion_prc_cum-0.5) == min(abs(depletion_prc_cum-0.5)))

## Figure: cumulative distribution function
ggplot(df.depletion.rank) +
  geom_line(aes(x=depletion_rank_prc, y=depletion_prc_cum, color=factor(year))) + 
  annotate("segment", x=-Inf, xend=df.lines$depletion_rank_prc, 
           y=0.5, yend=0.5,
           color=col.gray,
           linetype = "dashed") + 
  annotate("segment", x=df.lines$depletion_rank_prc, xend=df.lines$depletion_rank_prc, 
           y=-Inf, yend=0.5,
           color=c(col.cat.grn, col.cat.org, col.cat.red),
           linetype = "dashed") + 
  scale_x_continuous(name = "Percent of Wells", labels = scales::percent, 
                     limits=c(0,1), expand=c(0,0)) + 
  scale_y_continuous(name = "Percent of Total Streamflow Depletion", labels = scales::percent, 
                     limits=c(0,1), expand=c(0,0)) +
  scale_color_manual(name = "Years of\nPumping", 
                     values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  scale_fill_manual(name = "Years of\nPumping", 
                    values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background=element_blank(),
        plot.margin = unit(c(1.5,4.5,0,0), "mm")) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_PercentCDF.png"),
         width=95, height=85, units="mm")

#### calculate relative importance of different predictors
# calculate distance to closest stream for each grow
df.grow.closest <- 
  df.pairs %>% 
  # only use sites predicted to have wells
  subset(GrowNum %in% sf.grows$GrowNum[sf.grows$Well.rf.pred=="Yes"]) %>% 
  group_by(GrowNum) %>% 
  filter(dist_wellToStream_m == min(dist_wellToStream_m)) %>% 
  # some wells are equidistant from multiple stream segments; in this case, use highest Tr
  filter(Tr_bulk_m2d == max(Tr_bulk_m2d)) %>% 
  transform(SegNum_GrowNum = paste0(SegNum, "_", GrowNum))

# trim data frame
df.analysis <-
  df.depletion.total.YearMonthGrow %>% 
  subset(year %in% yrs.plot & month %in% mo.plot) %>% 
  left_join(df.grow.closest, by=c("GrowNum")) %>% 
  left_join(sf.grows[,c("GrowNum", "WaterUseTotal_m3y")], by="GrowNum")

ggplot(df.analysis, aes(y=depletion_m3d, x=dist_wellToStream_m, color=factor(year))) +
  geom_point(alpha=0.5) +
  stat_smooth(method="loess")

set.seed(1)
prc.sample <- 0.8
n.iter <- 100
for (yr in yrs.plot){
  # subset to only wells in this year
  df.analysis.yr <- subset(df.analysis, year==yr)
  
  # number to sample
  n.sample <- round(dim(df.analysis.yr)[1]*prc.sample)
  
  for (iter in 1:n.iter){
    
    ## depletion_m3d_grow = depletion caused in all segments by that site
    df.analysis.i <- dplyr::sample_n(df.analysis.yr, n.sample)
    
    # build statistical models and extract R2
    lm.full <- lm(depletion_m3d_grow ~ WaterUseTotal_m3y + dist_wellToStream_m + Tr_bulk_m2d + S_bulk + lmda_m2d + well_dtb_m, 
                  data=subset(df.analysis.i, year==yr & depletion_prc_grow > 0.001))
    
    lm.drop.water <- lm(depletion_m3d_grow ~ dist_wellToStream_m + Tr_bulk_m2d + S_bulk + lmda_m2d + well_dtb_m, 
                        data=subset(df.analysis.i, year==yr & depletion_prc_grow > 0.001))
    
    lm.drop.dist <- lm(depletion_m3d_grow ~ WaterUseTotal_m3y + Tr_bulk_m2d + S_bulk + lmda_m2d + well_dtb_m, 
                       data=subset(df.analysis.i, year==yr & depletion_prc_grow > 0.001))
    
    lm.drop.Tr <- lm(depletion_m3d_grow ~ WaterUseTotal_m3y + dist_wellToStream_m + S_bulk + lmda_m2d + well_dtb_m, 
                     data=subset(df.analysis.i, year==yr & depletion_prc_grow > 0.001))
    
    lm.drop.S <- lm(depletion_m3d_grow ~ WaterUseTotal_m3y + dist_wellToStream_m + Tr_bulk_m2d + lmda_m2d + well_dtb_m, 
                    data=subset(df.analysis.i, year==yr & depletion_prc_grow > 0.001))
    
    lm.drop.lmda <- lm(depletion_m3d_grow ~ WaterUseTotal_m3y + dist_wellToStream_m + Tr_bulk_m2d + S_bulk + well_dtb_m, 
                       data=subset(df.analysis.i, year==yr & depletion_prc_grow > 0.001))
    
    lm.drop.dtb <- lm(depletion_m3d_grow ~ WaterUseTotal_m3y + dist_wellToStream_m + Tr_bulk_m2d + S_bulk + lmda_m2d, 
                      data=subset(df.analysis.i, year==yr & depletion_prc_grow > 0.001))
    
    if (iter==1 & yr==yrs.plot[1]){
      df.fit <- data.frame(iter = iter, 
                           year = yr,
                           adj.R2.full = summary(lm.full)$adj.r.squared,
                           adj.R2.drop.dist = summary(lm.drop.dist)$adj.r.squared,
                           adj.R2.drop.water = summary(lm.drop.water)$adj.r.squared,
                           adj.R2.drop.Tr = summary(lm.drop.Tr)$adj.r.squared,
                           adj.R2.drop.S = summary(lm.drop.S)$adj.r.squared,
                           adj.R2.drop.lmda = summary(lm.drop.lmda)$adj.r.squared,
                           adj.R2.drop.dtb = summary(lm.drop.dtb)$adj.r.squared)
    } else {
      df.fit <- rbind(df.fit, 
                      data.frame(iter = iter, 
                                 year = yr,
                                 adj.R2.full = summary(lm.full)$adj.r.squared,
                                 adj.R2.drop.dist = summary(lm.drop.dist)$adj.r.squared,
                                 adj.R2.drop.water = summary(lm.drop.water)$adj.r.squared,
                                 adj.R2.drop.Tr = summary(lm.drop.Tr)$adj.r.squared,
                                 adj.R2.drop.S = summary(lm.drop.S)$adj.r.squared,
                                 adj.R2.drop.lmda = summary(lm.drop.lmda)$adj.r.squared,
                                 adj.R2.drop.dtb = summary(lm.drop.dtb)$adj.r.squared))
    }
    
    print(paste0("year ", yr, " iter ", iter, " complete"))
    
  }
}

# calculate change in R2 for each variable
df.fit$del.adj.R2.drop.dist <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.dist)/df.fit$adj.R2.full
df.fit$del.adj.R2.drop.water <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.water)/df.fit$adj.R2.full
df.fit$del.adj.R2.drop.Tr <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.Tr)/df.fit$adj.R2.full
df.fit$del.adj.R2.drop.S <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.S)/df.fit$adj.R2.full
df.fit$del.adj.R2.drop.lmda <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.lmda)/df.fit$adj.R2.full
df.fit$del.adj.R2.drop.dtb <- (df.fit$adj.R2.full - df.fit$adj.R2.drop.dtb)/df.fit$adj.R2.full

df.fit.summary <- 
  df.fit %>% 
  group_by(year) %>% 
  summarize(adj.R2.mean = mean(adj.R2.full, na.rm=T),
            adj.R2.std = sd(adj.R2.full, na.rm=T),
            change.dist.mean = mean(del.adj.R2.drop.dist, na.rm=T),
            change.dist.std = sd(del.adj.R2.drop.dist, na.rm=T),
            change.water.mean = mean(del.adj.R2.drop.water, na.rm=T),
            change.water.std = sd(del.adj.R2.drop.water, na.rm=T),
            change.Tr.mean = mean(del.adj.R2.drop.Tr, na.rm=T),
            change.Tr.std = sd(del.adj.R2.drop.Tr, na.rm=T),
            change.S.mean = mean(del.adj.R2.drop.S, na.rm=T),
            change.S.std = sd(del.adj.R2.drop.S, na.rm=T),
            change.lmda.mean = mean(del.adj.R2.drop.lmda, na.rm=T),
            change.lmda.std = sd(del.adj.R2.drop.lmda, na.rm=T),
            change.dtb.mean = mean(del.adj.R2.drop.dtb, na.rm=T),
            change.dtb.std = sd(del.adj.R2.drop.dtb, na.rm=T)) %>% 
  melt(id=c("year", "adj.R2.mean", "adj.R2.std")) %>% 
  transform(var = str_split_fixed(variable, pattern="[.]", n=3)[,2],
            metric = str_split_fixed(variable, pattern="[.]", n=3)[,3]) %>% 
  dplyr::select(year, value, var, metric) %>% 
  dcast(year + var ~ metric) %>% 
  transform(lab = factor(var, levels=c("water", "dist", "Tr", "S", "lmda", "dtb"), 
                         labels=c("Water Use", "Distance to\nClosest Stream", "Transmissivity", "Storativity", "Streambed\nConductance", "Depth to\nBedrock")))

df.fit.summary$error.min <- df.fit.summary$mean - df.fit.summary$std
df.fit.summary$error.max <- df.fit.summary$mean + df.fit.summary$std

df.fit.summary$error.min[df.fit.summary$error.min < 0] <- 0
df.fit.summary$error.max[df.fit.summary$error.max > 1] <- 1

## Figure: importance of different variables
ggplot(df.fit.summary, aes(x=lab, y=mean, fill=factor(year))) + 
  geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0, color=col.gray) +
  geom_errorbar(aes(ymin=error.min, ymax=error.max),
                width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(name="% Reduction in Adjusted R2", 
                     labels=scales::percent) +
  scale_x_discrete(name="Variable Dropped from Model") +
  scale_fill_manual(name = "Years of\nPumping", 
                    values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background=element_blank()) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_VariableImportance.png"),
         width=190, height=85, units="mm")

df.fit %>% 
  group_by(year) %>% 
  summarize(adj.R2.mean = mean(adj.R2.full, na.rm=T),
            adj.R2.std = sd(adj.R2.full, na.rm=T))
