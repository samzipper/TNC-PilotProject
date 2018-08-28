## Figure_Error-Ternary.R
#' This script is intended to compare estimated transient streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
#' and geometric depletion apportionment methods (from script Navarro_Analytical_Transient.R)
#' for all wells and reaches and make a ternary diagram.

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

## combine data
f.plot <- 0.001
df.plot <- rbind(df.all[(df.all$depletion.prc > f.plot) | (df.all$depletion.prc.modflow > f.plot), 
                        c("Time", "depletion.prc.modflow", "depletion.prc", "segments", "pump", "WellNum")],
                 df.max.all[(df.max.all$depletion.prc > f.plot) | (df.max.all$depletion.prc.modflow > f.plot), 
                            c("Time", "depletion.prc.modflow", "depletion.prc", "segments", "pump", "WellNum")])
df.plot$pump <- factor(df.plot$pump, levels=c("Transient", "Intermittent"), labels=c("Continuous\nPumping", "Intermittent\nPumping"))
df.plot$segments <- factor(df.plot$segments, levels=c("Most Affected Segment", "All Segments"),
                           labels=c("Most Affected Segment", "All Segments"))

# calculate fit
df.fit.plot <- 
  df.plot %>%
  group_by(Time, pump, segments) %>% 
  summarize(MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            RMSE.overall = rmse(depletion.prc, depletion.prc.modflow),
            MAE.overall = mae(depletion.prc, depletion.prc.modflow),
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
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))

#### make plots
p.tern <- 
  ggtern(subset(df.fit.plot, pump=="Continuous\nPumping" & segments=="Most Affected Segment"), 
         aes(x=MSE.bias.norm, y=MSE.var.norm, z=MSE.cor.norm, color=Time/365, shape=pump)) +
  geom_point(shape=21) +
  labs(x="% MSE due to Bias", y="% MSE due to Variability", z="% MSE due to Correlation") +
  scale_color_viridis(name="Time [years]", limits=c(0,10), breaks=seq(0,10,2)) +
  theme_rgbw(base_size=10, base_family="Arial") +
  theme(text=element_text(color="black"),
        axis.title=element_text(face="bold", size=rel(1)),
        axis.text=element_text(size=rel(1)),
        tern.axis.title.show=F,
        tern.panel.grid.major=element_line(linetype="solid", size=rel(0.5)),
        legend.position="bottom",
        legend.box="vertical",
        legend.box.margin=margin(0,0,0,0, "mm"),
        legend.background=element_blank(),
        legend.title=element_text(face="bold", size=rel(1)),
        legend.text=element_text(size=rel(1)),
        legend.spacing.y=unit(2, "mm"),
        legend.key=element_blank(),
        legend.text.align=0,
        strip.background=element_blank(),
        strip.text=element_text(face="bold.italic"),
        plot.margin=unit(c(-8,-8,-1,-8), "mm"),
        tern.axis.arrow.sep=0.085,
        tern.axis.arrow.start=0.15,
        tern.axis.arrow.finish=0.85,
        tern.axis.arrow.text=element_text(face="bold", size=rel(1)))
ggtern::ggsave(file.path("figures+tables", "Figure_Error-Ternary.png"),
               p.tern, width=95, height=95, units="mm")

p.tern.all <- 
  ggtern(df.fit.plot, 
       aes(x=MSE.bias.norm, y=MSE.var.norm, z=MSE.cor.norm, color=Time/365, shape=pump)) +
  geom_point(shape=21) +
  facet_grid(segments ~ pump) +
  labs(x="% MSE due to Bias", y="% MSE due to Variability", z="% MSE due to Correlation") +
  scale_color_viridis(name="Time [years]", limits=c(0,10), breaks=seq(0,10,2)) +
  theme_rgbw(base_size=10, base_family="Arial") +
  theme(text=element_text(color="black"),
        axis.title=element_text(face="bold", size=rel(1)),
        axis.text=element_text(size=rel(1)),
        tern.axis.title.show=F,
        tern.panel.grid.major=element_line(linetype="solid", size=rel(0.5)),
        legend.position="bottom",
        legend.box="vertical",
        legend.box.margin=margin(0,0,0,0, "mm"),
        legend.background=element_blank(),
        legend.title=element_text(face="bold", size=rel(1)),
        legend.text=element_text(size=rel(1)),
        legend.spacing.y=unit(2, "mm"),
        legend.key=element_blank(),
        legend.text.align=0,
        strip.background=element_blank(),
        strip.text=element_text(face="bold.italic"),
        tern.axis.arrow.sep=0.085,
        tern.axis.arrow.start=0.15,
        tern.axis.arrow.finish=0.85,
        tern.axis.arrow.text=element_text(face="bold", size=rel(1)))
ggtern::ggsave(file.path("figures+tables", "Figure_Error-Ternary_All.png"),
               p.tern.all, width=190, height=210, units="mm")
