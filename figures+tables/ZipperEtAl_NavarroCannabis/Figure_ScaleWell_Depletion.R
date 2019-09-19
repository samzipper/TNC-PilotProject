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


#### calculate relative importance of different predictors
# calculate distance to closest stream for each grow and bulk transmissivity
df.grow.closest <- 
  df.pairs %>% 
  # only use sites predicted to have wells
  subset(GrowNum %in% sf.grows$GrowNum[sf.grows$Well.rf.pred=="Yes"]) %>% 
  group_by(GrowNum) %>% 
  filter(dist_wellToStream_m == min(dist_wellToStream_m)) %>% 
  # some wells are equidistant from multiple stream segments; in this case, use highest Tr
  filter(Tr_bulk_m2d == max(Tr_bulk_m2d))

# trim data frame
df.analysis <-
  df.depletion.total.YearMonthGrow %>% 
  subset(year %in% yrs.plot & month %in% mo.plot) %>% 
  left_join(df.grow.closest, by=c("GrowNum")) %>% 
  left_join(sf.grows[,c("GrowNum", "WaterUseTotal_m3y")], by="GrowNum")

set.seed(1)
for (yr in yrs.plot){
  # subset to only wells in this year
  #df.analysis.yr <- subset(df.analysis, year==yr & depletion_prc_grow > 0.001)
  df.analysis.yr <- subset(df.analysis, year==yr)
  
  ## use relaimpo package
  ##  explanation: https://advstats.psychstat.org/book/mregression/importance.php
  # regression of all considered predictors
  lm.full <- lm(depletion_m3d_grow ~ WaterUseTotal_m3y + dist_wellToStream_m + Tr_bulk_m2d + lmda_m2d + well_dtb_m, 
                data=df.analysis.yr)
  
  # anova to select significant predictors
  anova.p <- anova(lm.full)$`Pr(>F)`
  anova.var <- row.names(anova(lm.full))
  var.keep <- anova.var[anova.p < 0.05 & is.finite(anova.p)]
  
  # build formula for boostrap relative importance
  vars.formula <- as.formula(paste0("depletion_m3d_grow ~ ", paste(var.keep, collapse ="+")))
  
  # bootstrap relative importance of retained variables
  rel.imp <- relaimpo::boot.relimp(formula = vars.formula, 
                                   data = df.analysis.yr,
                                   b = 1000, 
                                   type = "lmg")
  
  # extract important bits
  rel.imp.eval <- relaimpo::booteval.relimp(rel.imp)
  df.yr <- data.frame(var = var.keep,
                      R2.overall = rel.imp.eval$est["R2"],
                      R2.var = rel.imp.eval$est[paste0(var.keep, ".lmg")],
                      R2.var.upper =  rel.imp.eval$lmg.upper[1,],
                      R2.var.lower = rel.imp.eval$lmg.lower[1,],
                      year = yr)
  
  if (yr==yrs.plot[1]){
    df.fit <- df.yr
  } else {
    df.fit <- rbind(df.fit, df.yr)
  }
  
}

df.fit$var <- factor(df.fit$var, levels = c("WaterUseTotal_m3y", "dist_wellToStream_m", "Tr_bulk_m2d"),
                     labels =c("WaterUseTotal_m3y"="Water Use", 
                               "dist_wellToStream_m"="Distance to\nClosest Stream", 
                               "Tr_bulk_m2d"="Transmissivity"))

## Figure: cumulative distribution function
p.rank <-
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
        plot.margin = unit(c(1.5,4.5,0,1), "mm")) +
  #ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_PercentCDF.png"),
  #       width=95, height=85, units="mm") +
  NULL

## Figure: importance of different variables
p.vars <- 
  ggplot(df.fit, aes(x=var, y=R2.var, fill=factor(year))) + 
  geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0, color=col.gray) +
  geom_errorbar(aes(ymin=R2.var.lower, ymax=R2.var.upper),
                width=.2,
                position=position_dodge(.9)) +
  scale_y_continuous(name="Streamflow Depletion Variance Explained", 
                     labels=scales::percent_format(accuracy = 2)) +
  scale_x_discrete(name=NULL) +
  scale_fill_manual(name = "Years of\nPumping", 
                    values=c(col.cat.grn, col.cat.org, col.cat.red)) +
  theme(legend.position=c(1,1),
        legend.justification=c(1,1),
        legend.background=element_blank()) +
  #ggsave(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_VariableImportance.png"),
  #       width=190, height=85, units="mm") +
  NULL


## combined figures
plot_grid(p.rank, 
          p.vars, 
          nrow=1,
          align="tb",
          rel_widths=c(1,1),
          labels=c("(a)", "(b)"),
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = c(0.16, 0.14),
          label_y = 0.99) %>% 
  save_plot(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_PercentCDF+VariableImportance.png"),
            plot = .,
            nrow=1,
            base_width=190/25.4,
            base_height=95/25.4)

plot_grid(p.rank, 
          p.vars, 
          nrow=1,
          align="tb",
          rel_widths=c(1,1),
          labels=c("(a)", "(b)"),
          label_size = 10,
          label_fontfamily = "Arial",
          label_fontface = "plain",
          label_x = c(0.16, 0.14),
          label_y = 0.99) %>% 
  save_plot(file.path("figures+tables", "ZipperEtAl_NavarroCannabis", "Figure_ScaleWell_Depletion_PercentCDF+VariableImportance.pdf"),
            plot = .,
            nrow=1,
            base_width=190/25.4,
            base_height=95/25.4,
            device=cairo_pdf)
