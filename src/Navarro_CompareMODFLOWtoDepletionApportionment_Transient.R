## Navarro_CompareMODFLOWtoDepletionApportionment_Transient.R
#' This script is intended to compare estimated transient streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
#' and geometric depletion apportionment methods (from script Navarro_Analytical_Transient.R).

source(file.path("src", "paths+packages.R"))

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## choose modflow version
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"
masked <- F           # use masked depletion apportionment? should only be T for SFR

# which methods to analyze?
methods.plot <- c("Qf.InvDistSq", "Qf.WebSq", "Qf.TPoly")

## load MODFLOW results
df.MODFLOW.RIV <- 
  file.path("modflow","HTC", "Navarro", "Transient", "RIV", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,
            stream_BC = "RIV")

df.MODFLOW.SFR <- 
  file.path("modflow","HTC", "Navarro", "Transient", "SFR", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,
            stream_BC = "SFR")

df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)

## load analytical results
# remember: this has been trimmed to only stream segments in Navarro,
# and any segments with depletion <= 0.0001 (0.01%) have been removed
df.analytical <- 
  file.path("results", "Depletion_Analytical_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  dplyr::select(c("SegNum", "WellNum", "analytical", "Time", methods.plot))

## Time has long decimals; round before merging to ensure time match
df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
df.analytical$Time <- round(df.analytical$Time, 1)

## combine MODFLOW with analytical
# remove analytical results with depletion <= f.thres (should be same as Navarro_Analytical_Transient.R)
f.thres <- 0.0001  # =0.01%
df <- 
  full_join(df.analytical, df.MODFLOW[,c("stream_BC", "SegNum", "WellNum", "Time", "depletion.prc.modflow")], by=c("SegNum", "WellNum", "Time")) %>% 
  melt(id=c("stream_BC", "SegNum", "WellNum", "Time", "analytical", "depletion.prc.modflow"),
       value.name="depletion.prc", variable.name="method") %>% 
  subset(is.na(depletion.prc) | depletion.prc > f.thres)

# missing values should be 0 (0s were filtered out in previous scripts)
df$depletion.prc.modflow[is.na(df$depletion.prc.modflow)] <- 0
df$analytical[is.na(df$depletion.prc)] <- "hunt"
tmp <- subset(df, is.na(depletion.prc))
tmp$analytical <- "glover"
df <- rbind(df, tmp)
df$depletion.prc[is.na(df$depletion.prc)] <- 0

tmp <- subset(df, is.na(stream_BC))
tmp$stream_BC <- "RIV"
df <- rbind(df, tmp)
df$stream_BC[is.na(df$stream_BC)] <- "SFR"

## calculate fit
df.fit <- 
  df %>% 
  group_by(stream_BC, analytical, method, Time) %>% 
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
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))

## calculate sum of depletion.prc for each well/timestep/method
df.sum <-
  df %>% 
  group_by(stream_BC, Time, WellNum, analytical, method) %>% 
  summarize(Qf.MODFLOW = sum(depletion.prc.modflow),
            Qf.analytical = sum(depletion.prc))

df.fit.sum <- 
  df.sum %>% 
  group_by(stream_BC, analytical, method) %>% 
  summarize(n.reach = sum(is.finite(Qf.analytical)),
            cor = cor(Qf.analytical, Qf.MODFLOW, method="pearson"),
            bias = pbias(Qf.analytical, Qf.MODFLOW),
            R2 = R2(Qf.analytical, Qf.MODFLOW),
            MSE.bias = MSE.bias(Qf.analytical, Qf.MODFLOW),
            MSE.var = MSE.var(Qf.analytical, Qf.MODFLOW),
            MSE.cor = MSE.cor(Qf.analytical, Qf.MODFLOW),
            MSE.bias.norm = MSE.bias.norm(Qf.analytical, Qf.MODFLOW),
            MSE.var.norm = MSE.var.norm(Qf.analytical, Qf.MODFLOW),
            MSE.cor.norm = MSE.cor.norm(Qf.analytical, Qf.MODFLOW),
            MSE.overall = MSE(Qf.analytical, Qf.MODFLOW),
            KGE.overall = KGE(Qf.analytical, Qf.MODFLOW, method="2012"))

df.fit.sum.WellNum <- 
  df.sum %>% 
  group_by(stream_BC, analytical, method, WellNum) %>% 
  summarize(n.reach = sum(is.finite(Qf.analytical)),
            cor = cor(Qf.analytical, Qf.MODFLOW, method="pearson"),
            bias = pbias(Qf.analytical, Qf.MODFLOW),
            R2 = R2(Qf.analytical, Qf.MODFLOW),
            MSE.bias = MSE.bias(Qf.analytical, Qf.MODFLOW),
            MSE.var = MSE.var(Qf.analytical, Qf.MODFLOW),
            MSE.cor = MSE.cor(Qf.analytical, Qf.MODFLOW),
            MSE.bias.norm = MSE.bias.norm(Qf.analytical, Qf.MODFLOW),
            MSE.var.norm = MSE.var.norm(Qf.analytical, Qf.MODFLOW),
            MSE.cor.norm = MSE.cor.norm(Qf.analytical, Qf.MODFLOW),
            MSE.overall = MSE(Qf.analytical, Qf.MODFLOW),
            KGE.overall = KGE(Qf.analytical, Qf.MODFLOW, method="2012"))
  

## make plots
# some lines/annotations
ts.pump.start <- sum(days_in_month(seq(1,4))) + 1
ts.all <- unique(df$Time)
ts.plot <- ts.all[c(6, 73, 365)]

p.ts.KGE <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=method)) +
  geom_hline(yintercept=0, color=col.gray) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Kling-Gupta Efficiency") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.ts.KGE.png")),
       p.ts.KGE, width=8, height=8, units="in")

p.ts.MSE <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=MSE.overall, color=method)) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Mean Squared Error", labels=scales::percent) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.ts.MSE.png")),
       p.ts.MSE, width=8, height=8, units="in")

p.ts.bias <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=bias, color=method)) +
  geom_hline(yintercept=0, color=col.gray) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Percent Bias") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.ts.bias.png")),
       p.ts.bias, width=8, height=8, units="in")

p.ts.cor <-
  df.fit %>% 
  subset(method %in% methods.plot) %>% 
  ggplot(aes(x=Time, y=cor, color=method)) +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_grid(stream_BC~analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Correlation") +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.ts.cor.png")),
       p.ts.cor, width=8, height=8, units="in")

## scatterplots at different timesteps
p.depletion.scatter <-
  df %>% 
  subset(Time %in% ts.plot & 
           method %in% methods.plot) %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(Time ~ stream_BC+analytical) +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.depletion.scatter.png")),
       p.depletion.scatter, width=8, height=8, units="in")

p.depletion.scatter.limits <-
  df %>% 
  subset(Time %in% ts.plot & 
           method %in% methods.plot) %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(Time ~ stream_BC+analytical) +
  scale_x_continuous(name="Analytical Depletion Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", limits=c(0,1)) +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.depletion.scatter.limits.png")),
       p.depletion.scatter.limits, width=8, height=8, units="in")

# capture fraction is independent of depletion apportionment equation so only plot one
p.sum.scatter <-
  df.sum %>% 
  subset(method %in% methods.plot[1]) %>% 
  ggplot(aes(x=Qf.analytical, y=Qf.MODFLOW, color=Time-ts.pump.start)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(stream_BC ~ analytical, scales="free_y", 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Analytical Capture Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Capture Fraction") +
  scale_color_continuous(name="Time Since Start\nOf Pumping [days]") +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.sum.scatter.png")),
       p.sum.scatter, width=6, height=6, units="in")

p.sum.scatter.limits <-
  df.sum %>% 
  subset(method %in% methods.plot[1]) %>% 
  ggplot(aes(x=Qf.analytical, y=Qf.MODFLOW, color=Time-ts.pump.start)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(stream_BC ~ analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Analytical Capture Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Capture Fraction", limits=c(0,1)) +
  scale_color_continuous(name="Time Since Start\nOf Pumping [days]") +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.sum.scatter.limits.png")),
       p.sum.scatter.limits, width=6, height=6, units="in")

### some RIV only plots for presentation
p.depletion.scatter.limits <-
  df %>% 
  subset(Time %in% ts.plot & 
           method %in% methods.plot & 
           stream_BC=="RIV") %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point(shape=21) +
  facet_grid(analytical ~ Time) +
  scale_x_continuous(name="Analytical Depletion Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", limits=c(0,1)) +
  scale_color_manual(name=NULL, values=pal.method.Qf, labels=labs.method) +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.depletion.scatter.limits_RIV.png")),
       p.depletion.scatter.limits, width=8, height=6, units="in")

p.sum.scatter.limits <-
  df.sum %>% 
  subset(method %in% methods.plot[1] & 
           stream_BC=="RIV") %>% 
  ggplot(aes(x=Qf.analytical, y=Qf.MODFLOW, color=Time-ts.pump.start)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(. ~ analytical, 
             labeller=as_labeller(c(labs.analytical, labs.stream_BC))) +
  scale_x_continuous(name="Analytical Capture Fraction", limits=c(0,1)) +
  scale_y_continuous(name="MODFLOW Capture Fraction", limits=c(0,1)) +
  scale_color_continuous(name="Time Since Start\nOf Pumping [days]") +
  theme(legend.position="bottom")
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_Transient_p.sum.scatter.limits_RIV.png")),
       p.sum.scatter.limits, width=6, height=4, units="in")

df %>% 
  subset(Time %in% ts.plot & 
           method %in% methods.plot & 
           stream_BC=="RIV" &
           analytical=="glover") %>% 
  group_by(stream_BC, analytical, Time, method) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            bias = pbias(depletion.prc, depletion.prc.modflow),
            R2 = R2(depletion.prc, depletion.prc.modflow),
            #MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            #MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            #MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))
