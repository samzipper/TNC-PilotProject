## Figure_CompareAll.R
#' This figure compares all analytical depletion functions.

source(file.path("src", "paths+packages.R"))

#### (0) Prep various inputs, parameters, etc
## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

## which MODFLOW to plot
modflow_v <- "mfnwt"
stream_BC_plot <- "SFR"

## which conditions to plot
analytical_plot <- c("glover", "hunt")
method_plot <- c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly")
domain_plot <- c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic", "Adjacent+Dynamic")

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
            df.MODFLOW %>%
            subset(stream_BC == BC) %>%
            left_join(subset(df.analytical, method==m & analytical==a),
                      by=c("WellNum", "Time")) %>%
            replace_na(list("analytical"=a, "method"=m, "stream_BC"=BC, "depletion.prc"=0, "depletion.prc.modflow" = 0)) %>%
            transform(apportionment = apportionment_name,
                      pump = timeType,
                      stringsAsFactors=F) %>%
            group_by(stream_BC, pump, analytical, apportionment, method, Time) %>%
            summarize(MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
                      KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))
          
          ## overall fit, only > 5% reaches
          df.gt5 <-
            df.MODFLOW %>%
            subset(stream_BC == BC) %>%
            left_join(subset(df.analytical, method==m & analytical==a),
                      by=c("WellNum", "Time")) %>%
            replace_na(list("analytical"=a, "method"=m, "stream_BC"=BC, "depletion.prc"=0, "depletion.prc.modflow" = 0)) %>%
            subset(depletion.prc > 0.05 | depletion.prc.modflow > 0.05) %>% 
            transform(apportionment = apportionment_name,
                      pump = timeType,
                      stringsAsFactors=F) %>%
            group_by(stream_BC, pump, analytical, apportionment, method, Time) %>%
            summarize(MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
                      KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))
          
          ## capture fraction
          df.sum <-
            df.MODFLOW.sum %>%
            subset(stream_BC == BC) %>%
            left_join(subset(df.analytical.sum, method==m & analytical==a),
                      by=c("WellNum", "Time")) %>%
            replace_na(list("analytical"=a, "method"=m, "stream_BC"=BC, "Qf.total.analytical"=0, "Qf.total.modflow" = 0)) %>%
            transform(apportionment = apportionment_name,
                      pump = timeType,
                      stringsAsFactors=F) %>%
            group_by(stream_BC, pump, analytical, apportionment, method, Time) %>%
            summarize(MSE.sum = MSE(Qf.total.analytical, Qf.total.modflow),
                      RMSE.sum = rmse(Qf.total.analytical, Qf.total.modflow),
                      MAE.sum = mae(Qf.total.analytical, Qf.total.modflow),
                      Qf.total.MODFLOW.mean = mean(Qf.total.modflow),
                      Qf.total.MODFLOW.min = min(Qf.total.modflow),
                      Qf.total.MODFLOW.max = max(Qf.total.modflow),
                      KGE.sum = KGE(Qf.total.analytical, Qf.total.modflow, method="2012"))
          
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
          
          # calculate fit
          df.max.perf <- 
            df.max %>% 
            group_by(Time, stream_BC, pump, apportionment, analytical, method) %>% 
            summarize(n.reach = sum(is.finite(SegNum.modflow)),
                      n.match = sum(SegNum.modflow==SegNum.analytical),
                      n.noMatch.noAnalytical = sum(SegNum.analytical==9999),
                      MSE.match = MSE(depletion.prc, depletion.prc.modflow),
                      RMSE.match = rmse(depletion.prc, depletion.prc.modflow),
                      MAE.match = mae(depletion.prc, depletion.prc.modflow),
                      depletion.prc.modflow.mean = mean(depletion.prc.modflow),
                      depletion.prc.modflow.min = min(depletion.prc.modflow),
                      depletion.prc.modflow.max = max(depletion.prc.modflow),
                      KGE.match = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
            transform(prc.match = n.match/n.reach,
                      prc.noAnalytical = n.noMatch.noAnalytical/n.reach,
                      Streams = "Qd > 0.1%")
          
          df.max.perf.gt5 <-
            df.max %>% 
            subset(depletion.prc.modflow > 0.05) %>% 
            group_by(Time, stream_BC, pump, apportionment, analytical, method) %>% 
            summarize(n.reach = sum(is.finite(SegNum.modflow)),
                      n.match = sum(SegNum.modflow==SegNum.analytical),
                      n.noMatch.noAnalytical = sum(SegNum.analytical==9999),
                      MSE.match = MSE(depletion.prc, depletion.prc.modflow),
                      RMSE.match = rmse(depletion.prc, depletion.prc.modflow),
                      MAE.match = mae(depletion.prc, depletion.prc.modflow),
                      depletion.prc.modflow.mean = mean(depletion.prc.modflow),
                      depletion.prc.modflow.min = min(depletion.prc.modflow),
                      depletion.prc.modflow.max = max(depletion.prc.modflow),
                      KGE.match = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
            transform(prc.match = n.match/n.reach,
                      prc.noAnalytical = n.noMatch.noAnalytical/n.reach,
                      Streams = "Qd > 5%")
          
          if (start.flag) {
            df.fit.all <- df
            df.fit.gt5.all <- df.gt5
            df.fit.sum <- df.sum
            df.fit.match <- df.max.perf
            df.fit.match.gt5 <- df.max.perf.gt5
            start.flag <- F
          } else {
            df.fit.all <- rbind(df.fit.all, df)
            df.fit.gt5.all <- rbind(df.fit.gt5.all, df.gt5)
            df.fit.sum <- rbind(df.fit.sum, df.sum)
            df.fit.match <- rbind(df.fit.match, df.max.perf)
            df.fit.match.gt5 <- rbind(df.fit.match.gt5, df.max.perf.gt5)
          }
          
          # status update
          print(paste(timeType, apportionment_name, BC, m, a, "complete"))
          
        }  # end of a loop
      }  # end of m loop
    }  # end of BC loop
  }  # end of apportionment_name loop
}  # end of timeType loop

#df.match.plot <- rbind(df.fit.match, df.fit.match.gt5)
df.match.plot <- df.fit.match
df.match.plot$apportionment <- factor(df.match.plot$apportionment, 
                                      levels=c("WholeDomain", "LocalArea", "AdjacentOnly", "Dynamic", "Adjacent+Dynamic"))
df.match.plot$method <- factor(df.match.plot$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly"))
df.match.plot$pump <- factor(df.match.plot$pump, levels=c("Transient", "Intermittent"))
df.match.plot$analytical <- factor(df.match.plot$analytical, levels=c("glover", "hunt"))

## overall statistics
#df.overall.plot <- rbind(df.fit.all, df.fit.gt5.all)
df.overall.plot <- df.fit.all
df.overall.plot$apportionment <- factor(df.overall.plot$apportionment, levels=c("WholeDomain", "LocalArea", "AdjacentOnly", "Dynamic", "Adjacent+Dynamic"))
df.overall.plot$pump <- factor(df.overall.plot$pump, levels=c("Transient", "Intermittent"))
df.overall.plot$method <- factor(df.overall.plot$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly"))
df.overall.plot$analytical <- factor(df.overall.plot$analytical, levels=c("glover", "hunt"))

## capture fraction statistics
df.fit.sum$apportionment <- factor(df.fit.sum$apportionment, levels=c("WholeDomain", "LocalArea", "AdjacentOnly", "Dynamic", "Adjacent+Dynamic"))
df.fit.sum$pump <- factor(df.fit.sum$pump, levels=c("Transient", "Intermittent"))
df.fit.sum$method <- factor(df.fit.sum$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly"))
df.fit.sum$analytical <- factor(df.fit.sum$analytical, levels=c("glover", "hunt"))

## times for annotation
ts.pump.start <- sum(days_in_month(seq(1,4))) + 1 # for Transient continuous pumping

# for transient intermittent pumping
ts.pump.starts <- seq(from=(sum(days_in_month(seq(1,5))) + 1),
                      by=365, 
                      length.out=10)
ts.pump.stops <- seq(from=(sum(days_in_month(seq(1,10))) + 1),
                     by=365, 
                     length.out=10)

df.NoPump.times <- 
  rbind(data.frame(starts = c(0, ts.pump.stops),
                   stops = c(ts.pump.starts-1, 3650), 
                   pump = "Intermittent"),
        data.frame(starts = 0,
                   stops = ts.pump.start-1, 
                   pump = "Transient"))

## what is the winner for each variable to plot?
analytical_best <- "hunt"
method_best <- "Qf.WebSq"
domain_best <- "Adjacent+Dynamic"

## plot of match percentage
p.match.prc <-
  df.match.plot %>% 
  subset(stream_BC %in% stream_BC_plot & 
           apportionment %in% domain_plot &
           analytical %in% analytical_plot) %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_line(aes(x=Time/365, y=prc.match, group=interaction(apportionment, analytical, method)), color="gray75") +
  geom_line(data=subset(df.match.plot, analytical==analytical_best & apportionment==domain_best & method=="Qf.Web" & Streams=="Qd > 0.1%"), 
            aes(x=Time/365, y=prc.match), color=col.cat.red, size=2) +
  geom_line(data=subset(df.match.plot, analytical==analytical_best & apportionment==domain_best & method==method_best & Streams=="Qd > 0.1%"), 
            aes(x=Time/365, y=prc.match), color=col.cat.blu, size=2) +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="% of Wells where\nMost Affected Segment is\nCorrectly Identified", 
                     limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.25),
                     labels=scales::percent) +
  theme(legend.position="bottom",
        strip.text=element_blank())

## plot of fit, most affected reach
p.match.fit <- 
  df.match.plot %>% 
  subset(stream_BC %in% stream_BC_plot & 
           apportionment %in% domain_plot &
           analytical %in% analytical_plot &
           Streams=="Qd > 0.1%") %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=MAE.match/(depletion.prc.modflow.max-depletion.prc.modflow.min), group=interaction(apportionment, analytical, method)), color="gray75") +
  #geom_line(data=subset(df.match.plot, analytical==analytical_best & apportionment==domain_best & method==method_best & Streams=="Qd > 0.1%"), 
  #          aes(x=Time/365, y=depletion.prc.modflow.mean), color="black") +
  geom_line(data=subset(df.match.plot, analytical==analytical_best & apportionment==domain_best & method=="Qf.Web" & Streams=="Qd > 0.1%"), 
            aes(x=Time/365, y=MAE.match/(depletion.prc.modflow.max-depletion.prc.modflow.min)), color=col.cat.red, size=2) +
  geom_line(data=subset(df.match.plot, analytical==analytical_best & apportionment==domain_best & method==method_best & Streams=="Qd > 0.1%"), 
            aes(x=Time/365, y=MAE.match/(depletion.prc.modflow.max-depletion.prc.modflow.min)), color=col.cat.blu, size=2) +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_linetype_discrete(name="Segments\nEvaluated") +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_reverse(name="Normalized MAE,\nMost Affected Segment") +
  theme(strip.text=element_blank()) +
  NULL

## plot of KGE for all wells/reaches
p.overall <-
  df.overall.plot %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=KGE.overall, group=interaction(apportionment, analytical, method)), color="gray75") +
  geom_line(data=subset(df.overall.plot, analytical==analytical_best & apportionment==domain_best & method=="Qf.Web"),
            aes(x=Time/365, y=KGE.overall), color=col.cat.red, size=2) +
  geom_line(data=subset(df.overall.plot, analytical==analytical_best & apportionment==domain_best & method==method_best),
            aes(x=Time/365, y=KGE.overall), color=col.cat.blu, size=2) +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="(a) Continuous Pumping", "Intermittent"="(b) Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="KGE,\nAll Segments", limits=c(min(df.overall.plot$KGE.overall), 1)) +
  theme(strip.text=element_blank()) +
  NULL

## plot of MAE for cumulative capture fraction
p.sum <- 
  df.fit.sum %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=MAE.sum, group=interaction(apportionment, analytical, method)), color="gray75") +
  #geom_line(data=subset(df.fit.sum, analytical==analytical_best & apportionment==domain_best & method==method_best), 
  #          aes(x=Time/365, y=Qf.total.MODFLOW.mean), color="black") +
  geom_line(data=subset(df.fit.sum, analytical==analytical_best & apportionment==domain_best & method=="Qf.Web"),
            aes(x=Time/365, y=MAE.sum/(Qf.total.MODFLOW.max-Qf.total.MODFLOW.min)), color=col.cat.red, size=2) +
  geom_line(data=subset(df.fit.sum, analytical==analytical_best & apportionment==domain_best & method==method_best),
            aes(x=Time/365, y=MAE.sum/(Qf.total.MODFLOW.max-Qf.total.MODFLOW.min)), color=col.cat.blu, size=2) +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="(a) Continuous Pumping", "Intermittent"="(b) Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_reverse(name="Normalized MAE,\nCapture Fraction") +
  theme(strip.text=element_blank()) +
  NULL

## combine and save

# # all in one command
# # grab legend
# p.legend <- g_legend(p.match.prc + theme(legend.title=element_text(hjust=0.5)))
# ggsave(file.path("figures+tables", "Figure_CompareProximity.png"),
#        grid.arrange(
#          plot_grid(p.match.prc + theme(legend.position="none",
#                                        axis.title.x = element_blank(),
#                                        axis.text.x = element_blank()), 
#                    p.match.KGE + theme(legend.position="none",
#                                        axis.title.x = element_blank(),
#                                        axis.text.x = element_blank()),
#                    p.overall.KGE + theme(legend.position="none",
#                                          axis.title.x = element_blank(),
#                                          axis.text.x = element_blank()),
#                    p.sum.KGE + theme(legend.position="none"), 
#                    ncol=1, align="v", axis="l", labels="auto"),
#          p.legend, ncol=1, heights=c(4,0.25)),
#        width=190, height=200, units="mm")

# version without axis or legend which can be added with InkScape
save_plot(file.path("figures+tables", paste0("Figure_CompareAll_", stream_BC_plot, "_NoText.pdf")),
          plot_grid(p.match.prc + theme(legend.position="none",
                                        axis.title.x = element_blank(),
                                        axis.text.x = element_blank()),
                    p.match.fit + theme(legend.position="none",
                                        axis.title.x = element_blank(),
                                        axis.text.x = element_blank()),
                    p.overall + theme(legend.position="none",
                                      axis.title.x = element_blank(),
                                      axis.text.x = element_blank()),
                    p.sum + theme(legend.position="none",
                                  axis.title.x = element_blank(),
                                  axis.text.x = element_blank()),
                    ncol=1, align="v", axis="l"),
          ncol = 1, nrow = 4, base_width = 185/25.4, base_height=50/25.4, device=cairo_pdf)
