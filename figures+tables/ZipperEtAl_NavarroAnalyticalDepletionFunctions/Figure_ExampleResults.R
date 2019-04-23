## Figure_ExampleResults.R
#' Figure showing example of what streamflow depletion results look like.

source(file.path("src", "paths+packages.R"))

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## choose modflow version
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"
masked <- F           # use masked depletion apportionment? should only be T for SFR
timeType  <- "Transient" # "Transient" or "Intermittent"

# which methods to analyze?
method_plot <- c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly")

# what depletion apportionment output do you want?
domain_plot <- "Adjacent+Dynamic"

# what analytical model to plot?
analytical_plot <- "hunt"

# what MODFLOW stream BC?
stream_BC_plot <- c("RIV")

# what well to plot?
wel.plot <- 393  # 365, 393, 421

#### (0) Prep spatial data

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T) %>% 
  subset(WellNum %in% wel.plot)

sf.basin <-
  sf::st_read(file.path("data", "NHD", "WBD", "WBDHU10_Navarro.shp"))

#### (1) Load depletion apportionment results results and combine with MODFLOW
####     (Navarro_DepletionApportionment+Analytical_Transient.R)

# remember: this has been trimmed to only stream segments in Navarro,
# and any segments with depletion <= 0.0001 (0.01%) have been removed
df.analytical <- 
  paste0("Depletion_Analytical_", timeType, "_", domain_plot, "_AllMethods+Wells+Reaches.csv") %>% 
  file.path("results", .) %>% 
  read.csv(stringsAsFactors=F) %>% 
  dplyr::select(c("SegNum", "WellNum", "analytical", "Time", method_plot)) %>% 
  subset(Time == 3650 & analytical %in% analytical_plot & WellNum %in% wel.plot)

#### (3) Load MODFLOW results and figure out timesteps and wells for comparison 
####     (MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
df.MODFLOW <- 
  file.path("modflow","HTC", "Navarro", timeType, "RIV", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
  transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
            stream_BC = "RIV")

df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
df.analytical$Time <- round(df.analytical$Time, 1)

df.MODFLOW <- subset(df.MODFLOW, Time==3650 & WellNum %in% wel.plot)

# remove results with depletion <= f.thres (should be same as Navarro_Analytical_Transient.R)
f.thres <- 0.001  # =0.1%
df <- 
  full_join(df.analytical, 
            df.MODFLOW[,c("stream_BC", "SegNum", "WellNum", "Time", "depletion.prc.modflow")], 
            by=c("SegNum", "WellNum", "Time"))

sf.streams <- 
  sf::st_read(file.path("results", "GIS", "Navarro_Cannabis_StreamNetwork.shp"), stringsAsFactors=F) %>% 
  subset(TermnlP == outlet.TerminalPa) %>% 
  left_join(df, by="SegNum") %>% 
  dplyr::select(SegNum, Qf.Web, Qf.WebSq, Qf.InvDist, Qf.InvDistSq, Qf.TPoly, depletion.prc.modflow) %>% 
  melt(id=c("geometry", "SegNum"),
       value.name="depletion.prc", variable.name="method") %>%  
  replace_na(list(depletion.prc=0))

sf.streams$method <- factor(sf.streams$method, 
                            levels = c("depletion.prc.modflow", "Qf.TPoly", "Qf.InvDist", "Qf.InvDistSq", "Qf.Web", "Qf.WebSq"),
                            labels = c("(a) MODFLOW", "(b) Thiessen Polygon", "(c) Inverse Distance", "(d) Inverse Distance Squared", "(e) Web", "(f) Web Squared"))

zoom.breaks.x <- c(462000, 466000)
zoom.breaks.y <- c(4327000, 4332000)

if (wel.plot==365){
  shape.val <- 17
  shape.name <- "Proximate"
}
if (wel.plot==393){
  shape.val <- 15
  shape.name <- "Near"
}
if (wel.plot==421) {
  shape.val <- 18
  shape.name <- "Far"
}

ggplot(sf.streams) +
  geom_sf(aes(color=cut(depletion.prc, 
                        breaks=c(-0.01, 0.01, 0.05, 0.1, 0.25, 1.0),
                        labels = c("< 0.01", "0.01 - 0.05", "0.05 - 0.10", "0.10 - 0.25", "> 0.25"))),
          size = 2) +
  geom_point(data=df.wel, aes(x=lon, y=lat), color="black", shape=shape.val, size=2) +
  geom_text(data=df.wel, aes(x=lon, y=lat, label=shape.name), 
            color="black",  hjust = -0.2, vjust = 1.2) +
  geom_sf(data=sf.basin, color=col.gray, fill=NA) +
  facet_wrap(~method) + 
  scale_color_manual(name = "Depletion potential\nafter 10 years",
                     values = c(col.cat.blu, col.cat.grn, col.cat.yel, col.cat.org, col.cat.red),
                     drop = F) +
  scale_x_continuous(name="Easting [m]", expand=c(0,0),
                     breaks = zoom.breaks.x) +
  scale_y_continuous(name="Northing [m]", expand=c(0,0),
                     breaks = zoom.breaks.y) +
  coord_sf(crs=crs.MODFLOW, datum=crs.MODFLOW, 
           xlim = c(459500, 468000),    # limits from Figure_Error-IndividualWells.R
           ylim = c(4325000, 4335000)) +
  theme(panel.grid=element_line(color="transparent")) +
  theme(axis.text.y=element_text(angle=90, hjust=0.5),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.box.background=element_blank()) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroAnalyticalDepletionFunctions", paste0("Figure_ExampleResults_", shape.name, ".png")),
         width = 190, height=190, units="mm")
