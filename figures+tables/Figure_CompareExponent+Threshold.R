## Figure_CompareExponent+Threshold.R
#' This script is intended to compare different web exponents and expanding percent thresholds.
#' Requires output from Navarro_DecideBestMethod_05_CompareModels.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

#### (0) Prep various inputs, parameters, etc
## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

# which MODFLOW version? 
modflow_v <- "mfnwt"  # (no other option)

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%


## which exponents to plot?
web.exp_plot <- seq(1,3,0.5)

#### load fir statistics
## load data and plot
# fit based on most affected segment only
df.fit.match <- 
  read.csv(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match.csv"), 
           stringsAsFactors=F) %>% 
  subset(web.exp %in% web.exp_plot)
df.fit.match$pump <- factor(df.fit.match$pump, levels=c("Transient", "Intermittent"))

# fit based on all reaches
df.fit.all <- 
  read.csv(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-All.csv"), 
           stringsAsFactors=F) %>% 
  subset(web.exp %in% web.exp_plot)
df.fit.all$pump <- factor(df.fit.all$pump, levels=c("Transient", "Intermittent"))

#### plot
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

# color ramp for web.exp
pal.exp <- colorRampPalette(c("#e6194b","#0082c8","#d2f53c"))(6)

## based on Figure_CompareAll, for web.exp we should compare: 
#   -df.fit.match$MAE.match
#   -df.fit.all$bias.match
#   -df.fit.all$KGE.overall

p.match.fit <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_line(aes(x=Time/365, y=MAE.match/(depletion.prc.modflow.max-depletion.prc.modflow.mean), color=factor(web.exp))) +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="Normalized MAE,\nMost Affected Segment") +
  scale_color_manual(name="Web Exponent", values=pal.exp) +
  theme(legend.position="bottom") +
  NULL

p.match.bias <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=bias.match, color=factor(web.exp))) +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="% Bias of Depletion Potential,\nMost Affected Segment") +
  coord_cartesian(ylim=c(min(df.fit.match$bias.match), max(subset(df.fit.match, pump=="Intermittent")$bias.match))) +
  scale_color_manual(name="Web Exponent", values=pal.exp) +
  theme(legend.position="bottom") +
  NULL

p.overall.fit <-
  df.fit.all %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=KGE.overall, color=factor(web.exp))) +
  facet_wrap(pump ~ ., ncol=2, 
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="KGE, All Segments") +
  scale_color_manual(name="Web Exponent", values=pal.exp) +
  coord_cartesian(ylim=c(-1, 1)) +
  theme(legend.position="bottom") +
  NULL

# version without axis or legend which can be added with InkScape
save_plot(file.path("figures+tables", "Figure_CompareExponent+Threshold_Exponent_NoText.pdf"),
          plot_grid(p.match.fit + theme(legend.position="none",
                                        axis.title.x = element_blank(),
                                        axis.text.x = element_blank(),
                                        strip.text = element_blank()),
                    p.match.bias + theme(legend.position="none",
                                        axis.title.x = element_blank(),
                                        axis.text.x = element_blank(),
                                        strip.text = element_blank()),
                    p.overall.fit + theme(legend.position="none",
                                      axis.title.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      strip.text = element_blank()),
                    ncol=1, align="v", axis="l"),
          ncol = 1, nrow = 3, base_width = 185/25.4, base_height=60/25.4, device=cairo_pdf)


## based on Figure_CompareAll, for Qf.thres we should compare: 
#   -df.fit.all$KGE.overall

df.fit.all %>% 
  subset(web.exp == 2) %>% 
  ggplot() +
  geom_rect(data=df.NoPump.times, 
            aes(xmin=starts/365, xmax=stops/365, ymin=-Inf, ymax=Inf), 
            fill=col.gray, alpha=0.25) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line(aes(x=Time/365, y=KGE.overall, linetype=factor(Qf.thres))) +
  facet_wrap(pump ~ ., ncol=2, scales="free",
             labeller=as_labeller(c("Transient"="Continuous Pumping", "Intermittent"="Intermittent Pumping"))) +
  scale_x_continuous(name="Time [years]", expand=c(0,0), breaks=seq(0,10,2)) +
  scale_y_continuous(name="KGE, All Segments") +
  scale_linetype_manual(name="Adjacent + Expanding\nPercent Threshold", 
                        values=c("dotted", "dashed", "solid"),
                        labels=c("0.01%", "0.1%", "1%")) +
  coord_cartesian(ylim=c(-1, 1)) +
  theme(legend.position="bottom",
        legend.title=element_text(hjust=0.5),
        strip.text=element_text(face="bold"),
        plot.margin=unit(c(0.5, 2, 0.5, 0.5), "mm")) +
  ggsave(file.path("figures+tables", "Figure_CompareExponent+Threshold_Threshold.png"),
         width=190, height=95, units="mm") +
  NULL
