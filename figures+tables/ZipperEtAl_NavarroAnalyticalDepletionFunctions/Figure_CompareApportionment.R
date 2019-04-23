## Figure_CompareApportionment.R
#' This figure compares the different depletion apportionment equations.
#' It requires output from Navarro_DecideBestMethod_01_CompareDepletionApportionment.R

source(file.path("src", "paths+packages.R"))

#### (0) Prep various inputs, parameters, etc
## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

## which MODFLOW to plot
modflow_v <- "mfnwt"
stream_BC_plot <- c("RIV")

## which conditions to plot
analytical_plot <- "hunt"
method_plot <- c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly")
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
    
    ## load no apportionment output
    df.NoApport <- 
      file.path("results", paste0("Depletion_Analytical_", timeType, "_NoApportionment_AllMethods+Wells+Reaches.csv")) %>% 
      read.csv(stringsAsFactors=F) %>% 
      dplyr::select(SegNum, WellNum, Time, analytical, Qf.NoApport) %>% 
      set_colnames(c("SegNum", "WellNum", "Time", "analytical", "depletion.prc")) %>% 
      transform(method = "Qf.NoApport")
    
    ## add no apportionment to analytical
    df.analytical <-
      rbind(df.analytical, df.NoApport)
    
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
      for (m in c(method_plot, "Qf.NoApport")) {
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
          
          
          if (start.flag) {
            df.fit.all <- df
            df.fit.gt5.all <- df.gt5
            df.fit.sum <- df.sum
            df.max.all <- df.max
            start.flag <- F
          } else {
            df.fit.all <- rbind(df.fit.all, df)
            df.fit.gt5.all <- rbind(df.fit.gt5.all, df.gt5)
            df.fit.sum <- rbind(df.fit.sum, df.sum)
            df.max.all <- rbind(df.max.all, df.max)
          }
          
          # status update
          print(paste(timeType, apportionment_name, BC, m, a, "complete"))
          
        }  # end of a loop
      }  # end of m loop
    }  # end of BC loop
  }  # end of apportionment_name loop
}  # end of timeType loop

## calculate fit statistics
# calculate match percent
df.fit.match <-
  df.max.all %>% 
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

df.fit.match.gt5 <-
  df.max.all %>% 
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

df.match.plot <- rbind(df.fit.match, df.fit.match.gt5)
df.match.plot$apportionment <- factor(df.match.plot$apportionment, 
                                      levels=c("WholeDomain", "LocalArea", "AdjacentOnly", "Dynamic", "Adjacent+Dynamic"))
df.match.plot$method <- factor(df.match.plot$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly", "Qf.NoApport"))
df.match.plot$pump <- factor(df.match.plot$pump, levels=c("Transient", "Intermittent"))
df.match.plot$analytical <- factor(df.match.plot$analytical, levels=c("glover", "hunt"))

## overall statistics
#df.overall.plot <- rbind(df.fit.all, df.fit.gt5.all)
df.overall.plot <- df.fit.all
df.overall.plot$apportionment <- factor(df.overall.plot$apportionment, levels=c("WholeDomain", "LocalArea", "AdjacentOnly", "Dynamic", "Adjacent+Dynamic"))
df.overall.plot$pump <- factor(df.overall.plot$pump, levels=c("Transient", "Intermittent"))
df.overall.plot$method <- factor(df.overall.plot$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly", "Qf.NoApport"))
df.overall.plot$analytical <- factor(df.overall.plot$analytical, levels=c("glover", "hunt"))

## capture fraction statistics
df.fit.sum$apportionment <- factor(df.fit.sum$apportionment, levels=c("WholeDomain", "LocalArea", "AdjacentOnly", "Dynamic", "Adjacent+Dynamic"))
df.fit.sum$pump <- factor(df.fit.sum$pump, levels=c("Transient", "Intermittent"))
df.fit.sum$method <- factor(df.fit.sum$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly", "Qf.NoApport"))
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

## plot of match percentage
p.match.prc <-
  df.match.plot %>% 
  subset(stream_BC %in% stream_BC_plot & 
           apportionment %in% domain_plot &
           analytical %in% analytical_plot &
           Streams=="Qd > 5%" &
           method != "Qf.NoApport") %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_line(aes(x=Time/365, y=prc.match, color=method)) +
  geom_line(data = subset(df.match.plot, 
                          stream_BC %in% stream_BC_plot &
                            apportionment %in% domain_plot &
                            analytical %in% analytical_plot &
                            Streams=="Qd > 5%" &
                            method == "Qf.NoApport"),
            aes(x=Time/365, y=prc.match), linetype = "dashed", color = "black") +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_linetype_discrete(name="Segments\nEvaluated") +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="% of Wells where\nMost Affected Segment is\nCorrectly Identified", 
                     limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.25),
                     labels=scales::percent) +
  scale_color_manual(name="Depletion\nApportionment\nEquation", values=pal.method.Qf, labels=labels.method.Qf) +
  theme(legend.position="bottom",
        strip.text=element_blank()) +
  guides(colour = guide_legend(order=1, nrow=2), 
         linetype = guide_legend(order=2, nrow=2))

subset(df.match.plot, pump=="Transient" & method=="Qf.WebSq" & Streams=="Qd > 0.1%" & Time > 3300)[,c("Time", "prc.match")]
subset(df.match.plot, pump=="Intermittent" & method=="Qf.WebSq" & Streams=="Qd > 0.1%" & Time > 3300)[,c("Time", "prc.match")]
subset(df.match.plot, pump=="Transient" & method=="Qf.WebSq" & Streams=="Qd > 5%" & Time > 3300)[,c("Time", "prc.match")]
subset(df.match.plot, pump=="Intermittent" & method=="Qf.WebSq" & Streams=="Qd > 5%" & Time > 3300)[,c("Time", "prc.match")]

## plot of fit, most affected reach
p.match.fit <- 
  df.match.plot %>% 
  subset(stream_BC %in% stream_BC_plot & 
           apportionment %in% domain_plot &
           analytical %in% analytical_plot &
           Streams=="Qd > 0.1%" &
           method != "Qf.NoApport") %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=MAE.match/(depletion.prc.modflow.max-depletion.prc.modflow.min), color=method)) +
 # geom_line(data=subset(df.match.plot, analytical=="hunt" & apportionment=="Adjacent+Dynamic" & method=="Qf.WebSq" & Streams=="Qd > 0.1%"), 
#            aes(x=Time/365, y=depletion.prc.modflow.mean), color="black") +
  geom_line(data = subset(df.match.plot, 
                          stream_BC %in% stream_BC_plot &
                            apportionment %in% domain_plot &
                            analytical %in% analytical_plot &
                            Streams=="Qd > 0.1%" &
                            method == "Qf.NoApport"),
            aes(x=Time/365, y=MAE.match/(depletion.prc.modflow.max-depletion.prc.modflow.min)), linetype = "dashed", color = "black") +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_linetype_discrete(name="Segments\nEvaluated") +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_reverse(name="Normalized MAE,\nMost Affected Segment") +
  scale_color_manual(name="Depletion\nApportionment\nEquation", values=pal.method.Qf, labels=labels.method.Qf) +
  theme(strip.text=element_blank()) +
  NULL

subset(df.match.plot, pump=="Transient" & method=="Qf.WebSq" & Streams=="Qd > 0.1%" & Time > 3200)[,c("Time", "MAE.match", "depletion.prc.modflow.mean")]
subset(df.match.plot, pump=="Intermittent" & method=="Qf.WebSq" & Streams=="Qd > 0.1%" & Time > 3300)[,c("Time", "MAE.match", "depletion.prc.modflow.mean")]

## plot of KGE for all wells/reaches
p.overall <-
  df.overall.plot %>% 
  subset(method != "Qf.NoApport") %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=KGE.overall, color=method)) +
  geom_line(data = subset(df.overall.plot, 
                          stream_BC %in% stream_BC_plot &
                            apportionment %in% domain_plot &
                            analytical %in% analytical_plot &
                            method == "Qf.NoApport"),
            aes(x=Time/365, y=KGE.overall), linetype = "dashed", color = "black") +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="(a) Continuous Pumping", "Intermittent"="(b) Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="KGE,\nAll Segments") +
  coord_cartesian(ylim=c(-2.5,1)) +
  scale_color_manual(name="Depletion\nApportionment\nEquation", values=pal.method.Qf, labels=labels.method.Qf) +
  theme(strip.text=element_blank()) +
  NULL

## plot of fit for cumulative capture fraction
p.sum <- 
  df.fit.sum %>% 
  subset(stream_BC %in% stream_BC_plot &
           method != "Qf.NoApport") %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=MAE.sum/(Qf.total.MODFLOW.max-Qf.total.MODFLOW.min), color=method)) +
  #geom_line(data=subset(df.fit.sum, analytical=="hunt" & apportionment=="Adjacent+Dynamic" & method=="Qf.WebSq"), 
  #          aes(x=Time/365, y=Qf.total.MODFLOW.mean), color="black") +
  geom_line(data = subset(df.fit.sum, 
                          stream_BC %in% stream_BC_plot &
                            apportionment %in% domain_plot &
                            analytical %in% analytical_plot &
                            method == "Qf.NoApport"),
            aes(x=Time/365, y=MAE.sum/(Qf.total.MODFLOW.max-Qf.total.MODFLOW.min)), linetype = "dashed", color = "black") +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="(a) Continuous Pumping", "Intermittent"="(b) Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_reverse(name="Normalized MAE,\nCapture Fraction") +
  scale_color_manual(name="Depletion\nApportionment\nEquation", values=pal.method.Qf, labels=labels.method.Qf) +
  theme(strip.text=element_blank()) +
  NULL

subset(df.fit.sum, pump=="Transient" & method=="Qf.WebSq" & Time > 3200)[,c("Time", "MAE.sum", "Qf.total.MODFLOW.mean")]
subset(df.fit.sum, pump=="Intermittent" & method=="Qf.WebSq" & Time > 3300)[,c("Time", "MAE.sum", "Qf.total.MODFLOW.mean")]

qplot(MAE.sum, data=subset(df.fit.sum, pump=="Transient" & method=="Qf.WebSq" & Time > 3200))
qplot(MAE.sum/Qf.total.MODFLOW.mean, data=subset(df.fit.sum, pump=="Transient" & method=="Qf.WebSq" & Time > 3200))

qplot(MAE.sum, data=subset(df.fit.sum, pump=="Intermittent" & method=="Qf.WebSq" & Time > 3200))
qplot(MAE.sum/Qf.total.MODFLOW.mean, data=subset(df.fit.sum, pump=="Intermittent" & method=="Qf.WebSq" & Time > 3200))

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
save_plot(file.path("figures+tables", "ZipperEtAl_NavarroAnalyticalDepletionFunctions", "Figure_CompareApportionment_NoText.pdf"),
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


#### (2) Figure for SI: Comparison among steady-state results
## read in fit statistics
df.fit.SS <- read.csv(file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-SS.csv"),
                      stringsAsFactors=F)
df.fit.SS$apportionment <- factor(df.fit.SS$apportionment, levels=c("AdjacentOnly", "LocalArea", "WholeDomain"))
df.fit.SS$method <- factor(df.fit.SS$method, levels=c("f.Web", "f.WebSq", "f.InvDist", "f.InvDistSq", "f.TPoly"))

df.fit.SS %>% 
  subset(stream_BC %in% stream_BC_plot & 
           apportionment==domain_plot_SS) %>% 
  dplyr::select(method, apportionment, KGE.overall) %>% 
  melt(id=c("method", "apportionment"),
       variable.name="FitStatistic") %>%
  ggplot(aes(x=method, y=value, fill=method)) +
  geom_bar(stat="identity", position="dodge") +
  geom_hline(yintercept=0, color=col.gray) +
  scale_x_discrete(name="Depletion Apportionment Equation", 
                   labels=labels.method) +
  scale_y_continuous(name="Steady-State KGE") +
  #facet_wrap(~FitStatistic, scales="free_y", ncol=1) +
  scale_fill_manual(values=pal.method, guide=F) +
  ggsave(file.path("figures+tables", "Figure_SI_CompareApportionment_SteadyState.png"),
         width=95, height=60, units="mm") +
  NULL
