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

## which analytical to plot
analytical_plot <- "hunt"
domain_plot <- "Adjacent+Dynamic"
domain_plot_SS <- "WholeDomain"  # Adjacent+Dynamic doesn't exist for SS; plot WholeDomain instead

#### (1) Match % through time for Adjacent+Dynamic, Hunt
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

  df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)
  df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)

  ## for each stream_BC, WellNum, and Time find the most affected segment
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
      dplyr::select(SegNum, WellNum, Time, analytical, Qf.InvDist, Qf.InvDistSq, Qf.Web, Qf.WebSq, Qf.TPoly) %>%
      melt(id=c("SegNum", "WellNum", "Time", "analytical"),
           value.name="depletion.prc", variable.name="method") %>%
      subset(depletion.prc > f.thres)

    ## Time has long decimals; round before merging to ensure time match
    df.analytical$Time <- round(df.analytical$Time, 1)

    df.analytical.max <-
      df.analytical %>%
      group_by(analytical, WellNum, Time, method) %>%
      filter(depletion.prc==max(depletion.prc)) %>%
      dplyr::select(analytical, WellNum, Time, method, SegNum)

    for (BC in stream_BC_plot) {
      for (m in unique(df.analytical$method)) {
        for (a in analytical_plot) {
          # combine
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
            df.max.all <- df.max
            start.flag <- F
          } else {
            df.max.all <- rbind(df.max.all, df.max)
          }

          # status update
          print(paste(timeType, apportionment_name, BC, m, a, "complete"))

        }  # end of a loop
      }  # end of m loop
    }  # end of BC loop
  }  # end of apportionment_name loop
}  # end of timeType loop

## calculate match percent
df.fit.match <-
  df.max.all %>% 
  group_by(Time, stream_BC, pump, apportionment, analytical, method) %>% 
  summarize(n.reach = sum(is.finite(SegNum.modflow)),
            n.match = sum(SegNum.modflow==SegNum.analytical),
            n.noMatch.noAnalytical = sum(SegNum.analytical==9999),
            MSE.match = MSE(depletion.prc, depletion.prc.modflow),
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
            KGE.match = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
  transform(prc.match = n.match/n.reach,
            prc.noAnalytical = n.noMatch.noAnalytical/n.reach,
            Streams = "Qd > 5%")

# set factor levels
df.fit.match$apportionment <- factor(df.fit.match$apportionment, 
                                     levels=c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic", "Adjacent+Dynamic"))
df.fit.match$method <- factor(df.fit.match$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly"))
df.fit.match$pump <- factor(df.fit.match$pump, levels=c("Transient", "Intermittent"))


df.fit.match.gt5$apportionment <- factor(df.fit.match.gt5$apportionment, 
                                         levels=c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic", "Adjacent+Dynamic"))
df.fit.match.gt5$method <- factor(df.fit.match.gt5$method, levels=c("Qf.Web", "Qf.WebSq", "Qf.InvDist", "Qf.InvDistSq", "Qf.TPoly"))
df.fit.match.gt5$pump <- factor(df.fit.match.gt5$pump, levels=c("Transient", "Intermittent"))

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

## correct identification of most affected segment through time
df.plot <- rbind(df.fit.match, df.fit.match.gt5)

df.plot %>% 
  subset(stream_BC %in% stream_BC_plot & 
           apportionment==domain_plot &
           analytical==analytical_plot) %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts, xmax=stops, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.5) +
  geom_line(aes(x=Time, y=prc.match, color=method, linetype=Streams)) +
  facet_wrap(pump ~ ., ncol=1, scales="free", 
             labeller=as_labeller(c("Transient"="(a) Continuous Pumping", "Intermittent"="(b) Intermittent Pumping"))) +
  scale_linetype_discrete(name="Segments\nEvaluated") +
  scale_x_continuous(name="Time [days]", limits=c(0,3650), expand=c(0,0)) +
  scale_y_continuous(name="% of Wells where Most-Affected Reach is Correctly Identified", 
                     limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.25)) +
  scale_color_manual(name="Depletion\nApportionment\nEquation", values=pal.method.Qf, labels=labels.method.Qf) +
  theme(strip.text=element_text(hjust=0)) +
  guides(colour = guide_legend(order=1), 
         linetype = guide_legend(order=2)) +
  ggsave(file.path("figures+tables", "Figure_CompareApportionment_Match-prc.png"),
         width=190, height=120, units="mm") +
  NULL
  
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
