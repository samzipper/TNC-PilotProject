## Figure_Error-Scatterplots.R
#' This script is intended to compare estimated transient streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
#' and geometric depletion apportionment methods (from script Navarro_Analytical_Transient.R)
#' for all wells and reaches.

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

#### make plots
f.plot <- 0.001
df.plot <- rbind(df.all[(df.all$depletion.prc > f.plot) | (df.all$depletion.prc.modflow > f.plot), 
                        c("Time", "depletion.prc.modflow", "depletion.prc", "segments", "pump", "WellNum")],
                 df.max.all[(df.max.all$depletion.prc > f.plot) | (df.max.all$depletion.prc.modflow > f.plot), 
                            c("Time", "depletion.prc.modflow", "depletion.prc", "segments", "pump", "WellNum")])
df.plot$pump <- factor(df.plot$pump, levels=c("Transient", "Intermittent"), labels=c("Continuous", "Intermittent"))
df.plot$segments <- factor(df.plot$segments, levels=c("Most Affected Segment", "All Segments"),
                           labels=c("Most Affected Seg.", "All Segments"))

# all points, scatter with time as color
ggplot(df.plot, aes(x=depletion.prc, y=depletion.prc.modflow, color=Time/365)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  facet_rep_grid(segments ~ pump) +
  scale_x_continuous(name="Analytical Depletion Potential", limits=c(0,1), expand=c(0,0), 
                     breaks=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", " ")) +
  scale_y_continuous(name="MODFLOW Depletion Potential", limits=c(0,1), expand=c(0,0), 
                     breaks=seq(0,1,0.2), labels=c("0.0", "0.2", "0.4", "0.6", "0.8", " ")) +
  scale_color_viridis(name="Time [years]", limits=c(0,10), breaks=seq(0,10,2)) +
  coord_equal() +
  theme(strip.text=element_text(face="bold.italic"),
        legend.position="bottom",
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "mm")) +
  ggsave(file.path("figures+tables", "Figure_Error-Scatterplots.png"),
         width=95, height=110, units="mm")

df.plot %>%
  group_by(pump, segments) %>% 
  summarize(MSE.all = MSE(depletion.prc, depletion.prc.modflow),
            RMSE.all = rmse(depletion.prc, depletion.prc.modflow),
            MAE.all = mae(depletion.prc, depletion.prc.modflow),
            KGE.all = KGE(depletion.prc, depletion.prc.modflow, method="2012"))

## other plots not saved...
# all points, binned plot
ggplot(df.plot, aes(x=depletion.prc, y=depletion.prc.modflow)) +
  geom_bin2d(binwidth=0.01) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  stat_smooth(method="lm") +
  facet_grid(segments ~ pump) +
  scale_x_continuous(name="Analytical Depletion Potential", limits=c(0,1), expand=c(0,0)) +
  scale_y_continuous(name="MODFLOW Depletion Potential", limits=c(0,1), expand=c(0,0)) +
  scale_fill_viridis(trans="log10")

# selected times, scatterplot
times_all <- unique(df.plot$Time)
times_all[times_all > 270 & times_all < 280]
times_plot <- c(times_all[times_all > 270 & times_all < 280],  # year 1 oct 1, 
                times_all[times_all > 480 & times_all < 490],  # year 2 may 1
                max(times_all))  # end of transient simulation
df.plot %>% 
  subset(pump=="Transient" & Time %in% times_plot) %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow)) +
  geom_point(shape=21) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  stat_smooth(method="lm") +
  facet_grid(Time ~ segments) +
  scale_x_continuous(name="Analytical Depletion Potential", limits=c(0,1), expand=c(0,0)) +
  scale_y_continuous(name="MODFLOW Depletion Potential", limits=c(0,1), expand=c(0,0)) +
  scale_fill_viridis(trans="log10")

#df.test <- 
  df %>% 
  subset(pump=="Transient") %>% 
  group_by(Time, WellNum, pump, stream_BC, analytical, method, apportionment) %>% 
  summarize(depletion.prc.sum = sum(depletion.prc),
            depletion.prc.modflow.sum = sum(depletion.prc.modflow)) %>% 
    qplot(depletion.prc.sum, data=.)
 # head(20)
