## Navarro_CompareMODFLOWtoDepletionApportionment.R
#' This script is intended to compare estimated streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_Riv-SummarizeLeakage.py)
#' and geometric depletion apportionment methods.

source("src/paths+packages.R")

## what is the pumping rate?
Qw <- -5000  # [m3/d]

## choose stream boundary condition
stream_BC <- "RIV"  # RIV or SFR (have not processed SFR data yet)

## load stream shapefile
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# figure out which segments are part of Navarro
segs.navarro <- shp.streams@data$SegNum[shp.streams@data$TerminalPa==outlet.TerminalPa]

## load geometric data (from script Navarro_DepletionApportionment.R)
df.apportion <- 
  read.csv(file.path("results","Navarro_DepletionApportionment_AllMethods+Wells+Reaches.csv"),
           stringsAsFactors=F) %>% 
  subset(SegNum %in% segs.navarro)

## load MODFLOW output
df.MODFLOW <- 
  file.path("modflow","HTC", "Navarro", "SteadyState", stream_BC, "RIV-SummarizeLeakage.csv") %>% 
  read.csv(stringsAsFactors=F)

## calculate MODFLOW depletion
# make separate column for no pumping depletion
df.MODFLOW<- 
  df.MODFLOW %>% 
  subset(WellNum==0, select=c("SegNum", "leakage")) %>% 
  set_colnames(c("SegNum", "leakage_NoPump")) %>% 
  left_join(subset(df.MODFLOW, WellNum != 0), ., by=c("SegNum")) %>% 
  arrange(WellNum, SegNum)

# calculate change in leakage
# convention: negative value for leakage means gaining stream;
# so if (leakage_NoPump - leakage) < 0, stream is depleted due to pumping
df.MODFLOW$depletion_m3.d <- df.MODFLOW$leakage_NoPump - df.MODFLOW$leakage

df.MODFLOW <- 
  df.MODFLOW %>% 
  group_by(WellNum) %>% 
  summarize(leakage.sum = sum(leakage),
            leakage_NoPump.sum = sum(leakage_NoPump),
            depletion.sum = sum(depletion_m3.d)) %>% 
  left_join(df.MODFLOW, ., by="WellNum") %>% 
  subset(SegNum %in% segs.navarro)

# calculate depletion as percentage of depletion.sum
df.MODFLOW$depletion.prc.modflow <- df.MODFLOW$depletion_m3.d/df.MODFLOW$depletion.sum

## combine MODFLOW with analytical
df <- 
  left_join(df.apportion, df.MODFLOW[,c("SegNum", "WellNum", "depletion.prc.modflow")], by=c("SegNum", "WellNum")) %>% 
  subset(is.finite(depletion.prc.modflow)) %>% 
  melt(id=c("SegNum", "WellNum", "distToWell.min.m", "depletion.prc.modflow"),
       value.name="depletion.prc", variable.name="method")

## make plots
p.depletion <-
  ggplot(df, aes(x=depletion.prc, y=depletion.prc.modflow)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point(shape=21) +
  stat_smooth(method="lm") +
  facet_wrap(~method, labeller=as_labeller(labels.method)) +
  labs(title="Comparison of Depletion Apportionment Equations to Steady-State MODFLOW",
       subtitle="Each dot = 1 stream segment response to 1 pumping well") +
  scale_x_continuous(name="Analytical Depletion Fraction", limits=c(0,1), expand=c(0,0)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", expand=c(0,0))
ggsave(file.path("results", "Navarro_CompareMODFLOWtoDepletionApportionment_p.depletion.png"),
       p.depletion, width=8, height=6, units="in")

# calculate fit statistics
df.fit <- 
  df %>% 
  group_by(method) %>% 
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

# print to screen
df.fit[,c("method", "KGE.overall")]
