## Figure_Error-Factors.R
#' This script is intended to compare estimated error to different landscape factors.

source(file.path("src", "paths+packages.R"))

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

## which MODFLOW to plot
modflow_v <- "mfnwt"
stream_BC_plot <- c("RIV")

## which conditions to plot
analytical_plot <- "hunt"
method_plot <- "Qf.WebSq"
domain_plot <- "Adjacent+Dynamic"

## various model parameters
# units: [m] and [d]
# flow parameters
hk <- 1e-12*1e7*86400  # horizontal K [m/d], convert k [m-2] to K [m/s] to K [m/d]
ss <- 1e-5             # specific storage
sy <- 0.10             # specific yield (using 50% of domain mean porosity for now)
vka <- 10              # anisotropy
vk <- hk/vka           # calculate vertical K [m/d] based on horizontal K and anisotropy
screen_length <- 50    # well screen length [m]

## streambed parameters
depth <- 5  # river depth?
riverbed_K <- hk/10
riverbed_thickness <- 1

## Prep spatial data

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario - this is used to define the screen interval
## so that screen interval is consistent between MODFLOW and analytical
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", "RIV", modflow_v, "wte.csv")), header=F)

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$wte_m <- m.wte[as.matrix(df.wel[,c("row", "col")])+1]

# make a spatial points data frame
xy <- df.wel[,c("lon", "lat")]
spdf.wel <- SpatialPointsDataFrame(coords = xy, data = df.wel,
                                   proj4string = CRS(crs.MODFLOW))

## load stream data - created in MODFLOW_Navarro_InputPrepData.R
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# domain boundary shapefile
shp <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro")
shp.UTM <- spTransform(shp, crs.MODFLOW)

shp.adj <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

## prep polygon boundaries for plots
df.basin <- tidy(shp.UTM)
df.basin.adj <- tidy(shp.adj.UTM)
df.riv <- tidy(shp.streams)

## process data
start.flag <- T
for (timeType in c("Transient", "Intermittent")) {
  
  ## open MODFLOW results
  df.MODFLOW.RIV <-
    file.path("modflow","HTC", "Navarro", timeType, "RIV", modflow_v, "Depletion_MODFLOW.csv") %>%
    read.csv(stringsAsFactors=F) %>%
    #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
    transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
              stream_BC = "RIV",
              stringsAsFactors=F)
  
  df.MODFLOW.SFR <-
    file.path("modflow","HTC", "Navarro", timeType, "SFR", modflow_v, "Depletion_MODFLOW.csv") %>%
    read.csv(stringsAsFactors=F) %>%
    #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
    transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
              stream_BC = "SFR",
              stringsAsFactors=F)
  
  df.MODFLOW <-
    rbind(df.MODFLOW.RIV, df.MODFLOW.SFR) %>%
    subset(depletion.prc.modflow > f.thres) %>% 
    subset(stream_BC %in% stream_BC_plot)
  df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
  
  ## calculate capture fraction
  df.MODFLOW.sum <- 
    df.MODFLOW %>%
    group_by(WellNum, Time, stream_BC) %>%
    summarize(Qf.total.modflow = sum(depletion.prc.modflow))
  
  ## find most affected segment
  df.MODFLOW.most <-
    df.MODFLOW %>%
    subset(depletion.prc.modflow > f.thres) %>%
    group_by(stream_BC, WellNum, Time) %>%
    filter(depletion.prc.modflow==max(depletion.prc.modflow))
  
  for (apportionment_name in domain_plot) {
    ## load analytical output
    df.analytical <-
      paste0("Depletion_Analytical_", timeType, "_", apportionment_name, "_AllMethods+Wells+Reaches.csv") %>%
      file.path("results", .) %>%
      read.csv(stringsAsFactors=F) %>%
      subset(analytical %in% analytical_plot) %>%
      dplyr::select(c("SegNum", "WellNum", "Time", "analytical", method_plot)) %>%
      melt(id=c("SegNum", "WellNum", "Time", "analytical"),
           value.name="depletion.prc", variable.name="method") %>%
      subset(depletion.prc > f.thres)
    
    ## Time has long decimals; round before merging to ensure time match
    df.analytical$Time <- round(df.analytical$Time, 1)
    
    ## capture fraction
    df.analytical.sum <- 
      df.analytical %>% 
      group_by(WellNum, Time, analytical, method) %>%
      summarize(Qf.total.analytical = sum(depletion.prc))
    
    ## most affected segment
    df.analytical.max <-
      df.analytical %>%
      group_by(analytical, WellNum, Time, method) %>%
      filter(depletion.prc==max(depletion.prc)) %>%
      dplyr::select(analytical, WellNum, Time, method, SegNum)
    
    for (BC in stream_BC_plot) {
      for (m in method_plot) {
        for (a in analytical_plot) {
          ## overall fit
          df <-
            full_join(df.MODFLOW, df.analytical, by=c("WellNum", "Time", "SegNum")) %>%
            replace_na(list("analytical"=a, "method"=m, "stream_BC"=BC, "depletion.prc"=0, "depletion.prc.modflow" = 0)) %>%
            transform(apportionment = apportionment_name,
                      pump = timeType,
                      stringsAsFactors=F)
          
          ## capture fraction
          df.sum <-
            df.MODFLOW.sum %>%
            subset(stream_BC == BC) %>%
            left_join(df.analytical.sum,
                      by=c("WellNum", "Time")) %>%
            replace_na(list("analytical"=a, "method"=m, "stream_BC"=BC, "Qf.total.analytical"=0, "Qf.total.modflow" = 0)) %>%
            transform(apportionment = apportionment_name,
                      pump = timeType,
                      stringsAsFactors=F)
          
          ## most affected segment
          df.max <-
            df.MODFLOW.most %>%
            subset(stream_BC == BC) %>%
            dplyr::select(SegNum, WellNum, Time, stream_BC, depletion.prc.modflow) %>%
            left_join(subset(df.analytical.max, method==m & analytical==a),
                      by=c("WellNum", "Time"), suffix=c(".modflow", ".analytical")) %>%
            # add depletion in the most affected MODFLOW segment
            left_join(subset(df.analytical, method==m & analytical==a),
                      by=c("WellNum", "Time", "SegNum.modflow"="SegNum", "analytical", "method")) %>%
            replace_na(list("analytical"=a, "method"=m, "SegNum.analytical"=9999, "depletion.prc" = 0)) %>%
            transform(apportionment = apportionment_name,
                      pump = timeType)
          
          
          if (start.flag) {
            df.all <- df
            df.sum.all <- df.sum
            df.max.all <- df.max
            start.flag <- F
          } else {
            df.all <- rbind(df.all, df)
            df.sum.all <- rbind(df.sum.all, df.sum)
            df.max.all <- rbind(df.max.all, df.max)
          }
          
          # status update
          print(paste(timeType, apportionment_name, BC, m, a, "complete"))
          
        }  # end of a loop
      }  # end of m loop
    }  # end of BC loop
  }  # end of apportionment_name loop
}  # end of timeType loop

# add column for domain analyzed
df.all$segments <- "All Segments"
df.max.all$segments <- "Most Affected Segment"

# which do you want to use for plots - most affected segments, or all segments?
df.plot <- df.all  # df.all or df.max.all
# what times do you want to plot?
times_all <- unique(df.plot$Time)
times_plot <- times_all[times_all >= (3650-365)]

# what pumping schedule do you want to plot?
pump_plot <- "Transient"

## calculate overall fit by well
df.fit.well <-
  df.plot %>% 
  subset(Time %in% times_plot &
           pump %in% pump_plot) %>% 
  group_by(pump, WellNum) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MAE.overall = mae(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            MAE.overall = mae(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"),
            depletion.prc.modflow.mean = mean(depletion.prc.modflow),
            depletion.prc.modflow.max = max(depletion.prc.modflow),
            depletion.prc.modflow.min = min(depletion.prc.modflow)) %>% 
  left_join(df.wel, by="WellNum")

## fit by distance to well
# load horizontal distance
df.dist <-
  file.path("results", "Navarro_DepletionApportionment_WholeDomain_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  left_join(df.wel[,c("WellNum", "ztop_m")], by=c("WellNum"))

# load streambed info to calculate vertical distance
df.stream.elev <- 
  read.table(file.path("modflow", "input", "isfr_ReachData.txt"), stringsAsFactors=F, header=T) %>% 
  group_by(SegNum) %>% 
  summarize(totalStreamLength_m = sum(length_m),
            streambed_elev_m = median(elev_m_min)-depth)

df.dist <- 
  left_join(df.dist, df.stream.elev, by=c("SegNum")) %>% 
  transform(dist.vert.m = (ztop_m - streambed_elev_m))

# join with all info
df.all.dist <-
  df.plot %>% 
  subset(pump=="Transient") %>% 
  left_join(df.dist[,c("SegNum", "WellNum", "distToWell.min.m", "dist.vert.m")], by=c("SegNum", "WellNum"))

dist.breaks <- quantile(df.all.dist$distToWell.min.m, seq(0,1,0.05))
df.all.dist$dist.cut <- cut(df.all.dist$distToWell.min.m, dist.breaks, include.lowest=T)
df.fit.dist <- 
  df.all.dist %>% 
  subset(Time %in% times_plot &
           pump %in% pump_plot) %>% 
  group_by(pump, dist.cut) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            MAE.overall = mae(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"),
            depletion.prc.modflow.mean = mean(depletion.prc.modflow),
            depletion.prc.modflow.max = max(depletion.prc.modflow),
            depletion.prc.modflow.min = min(depletion.prc.modflow),
            dist.m.median = median(distToWell.min.m),
            dist.m.min = min(distToWell.min.m),
            dist.m.max = max(distToWell.min.m))

dist.vert.breaks <- quantile(df.dist$dist.vert.m, seq(0,1,0.05))
df.all.dist$dist.vert.cut <- cut(df.all.dist$dist.vert.m, dist.vert.breaks, include.lowest=T)
df.fit.dist.vert <- 
  df.all.dist %>% 
  subset(Time %in% times_plot &
           pump %in% pump_plot) %>% 
  group_by(pump, dist.vert.cut) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            MAE.overall = mae(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"),
            depletion.prc.modflow.mean = mean(depletion.prc.modflow),
            depletion.prc.modflow.max = max(depletion.prc.modflow),
            depletion.prc.modflow.min = min(depletion.prc.modflow),
            dist.m.median = median(dist.vert.m),
            dist.m.min = min(dist.vert.m),
            dist.m.max = max(dist.vert.m))

## calculate overall fit by reach length
df.fit.length <-
  df.stream.elev %>% 
  transform(length.cut = cut(totalStreamLength_m, 
                             quantile(totalStreamLength_m, seq(0,1,0.05)), 
                             include.lowest=T)) %>%  # cut into 10 percentile chunks
  #transform(length.cut = cut(totalStreamLength_m, c(0, 0.1, 0.5, 1, 5, max(totalStreamLength_m))*1000)) %>%  # cut at same length intervals as Nanaimo paper
  left_join(df.plot, ., by=c("SegNum")) %>% 
  subset(Time %in% times_plot &
           pump %in% pump_plot) %>% 
  group_by(pump, length.cut) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            MAE.overall = mae(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"),
            depletion.prc.modflow.mean = mean(depletion.prc.modflow),
            depletion.prc.modflow.max = max(depletion.prc.modflow),
            depletion.prc.modflow.min = min(depletion.prc.modflow),
            length.m.median = median(totalStreamLength_m),
            length.m.min = min(totalStreamLength_m),
            length.m.max = max(totalStreamLength_m))

##### make plots
## map of fit by well
p.map.well <-
  df.fit.well %>% 
  ggplot() +
  geom_polygon(data=df.basin, aes(x=long, y=lat, group=group), fill=NA, color="red") +
  geom_path(data=df.riv, aes(x=long, y=lat, group=group), color="blue") +
  geom_point(aes(x=lon, y=lat, size=MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min)), alpha=0.5) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0), breaks=map.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0), breaks=map.breaks.y) +
  scale_size_area(name="Normalized\nMAE") +
  coord_equal() +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.title=element_text(hjust=0.5))

## scatter fit by WTE
p.fit.well <-
  df.fit.well %>% 
  ggplot(aes(x=wte_m, y=MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min))) +
  geom_point() +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Water Table Elevation at Well [m]") +
  scale_y_continuous(name="Normalized MAE") +
  theme(legend.position="bottom")

lm(MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min) ~ wte_m, data=df.fit.well) %>% 
  summary()

save_plot(file.path("figures+tables", "Figure_Error-Factors_ByWell.png"),
          plot_grid(p.map.well + theme(legend.position="top"),
                    p.fit.well,
                    ncol = 1, 
                    nrow = 2, 
                    align = 'v',
                    axis = 'l',
                    rel_heights = c(2, 1),
                    labels = c("(a)", "(b)"),
                    label_size = 10,
                    label_fontface = "plain",
                    label_colour = "grey30",
                    hjust = 0,
                    vjust = 1),
          base_width = 95/25.4, base_height=150/25.4)


## barplot fit by horizontal distance
p.bar.dist <-
  df.fit.dist %>% 
  ggplot(aes(x=dist.m.median, y=MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min))) +
  geom_hline(yintercept=0, color="gray65") +
  geom_point() +
  # geom_errorbarh(aes(xmin=dist.m.min, xmax=dist.m.max), height=0) +
  stat_smooth(method="lm") +
  scale_x_continuous(name="Lateral Well-Stream Distance [m]") +
  scale_y_continuous(name="Normalized MAE")

lm(MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min) ~ dist.m.median, data=df.fit.dist) %>% 
  summary()

df.fit.dist$MAE.overall/(df.fit.dist$depletion.prc.modflow.max-df.fit.dist$depletion.prc.modflow.min)
df.fit.dist$dist.cut

## barplot fit by vertical distance
p.bar.vert <- 
  df.fit.dist.vert %>% 
  ggplot(aes(x=dist.m.median, y=MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min))) +
  geom_hline(yintercept=0, color="gray65") +
  geom_point() +
  # geom_errorbarh(aes(xmin=dist.m.min, xmax=dist.m.max), height=0) +
  stat_smooth(method="loess") +
  scale_x_continuous(name="Vertical Well-Stream Distance [m]") +
  scale_y_continuous(name="Normalized MAE")

lm(MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min) ~ dist.m.median, data=df.fit.dist.vert) %>% 
  summary()

## barplot fit by reach length
p.bar.length <- 
  df.fit.length %>% 
  ggplot(aes(x=length.m.median, y=MAE.overall/(depletion.prc.modflow.max-depletion.prc.modflow.min))) +
  geom_hline(yintercept=0, color="gray65") +
  geom_point() +
  stat_smooth(method="loess") +
  #  geom_errorbarh(aes(xmin=length.m.min, xmax=length.m.max), height=0) +
  scale_x_continuous(name="Stream Segment Length [m]") +
  scale_y_continuous(name="Normalized MAE") +
  coord_cartesian(ylim=c(0, max(df.fit.length$MAE.overall/(df.fit.length$depletion.prc.modflow.max-df.fit.length$depletion.prc.modflow.min))))

save_plot(file.path("figures+tables", "Figure_Error-Factors_Geometry.png"),
          plot_grid(p.bar.dist + theme(plot.margin=unit(c(0.5, 0.5, 2, 0.5), "mm")),
                    p.bar.vert + theme(plot.margin=unit(c(1, 0.5, 1, 0.5), "mm")),
                    p.bar.length + theme(plot.margin=unit(c(2, 0.5, 0.5, 0.5), "mm")),
                    ncol = 1, 
                    nrow = 3, 
                    align = 'v',
                    axis = 'l',
                    labels = c("(a)", "(b)", "(c)"),
                    label_size = 10,
                    label_fontface = "plain",
                    label_colour = "grey30",
                    hjust = 0,
                    vjust = 1),
          base_width = 95/25.4, base_height=150/25.4)
