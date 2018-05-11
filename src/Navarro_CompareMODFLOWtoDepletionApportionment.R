## Navarro_CompareMODFLOWtoDepletionApportionment.R
#' This script is intended to compare estimated streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_Riv-SummarizeLeakage.py)
#' and geometric depletion apportionment methods.

source(file.path("src", "paths+packages.R"))

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## choose stream boundary condition and modflow version
stream_BC <- "SFR"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"
masked <- T           # use masked depletion apportionment? should only be T for SFR
if (stream_BC != "SFR" & masked){
  stop("Error: Using masked depletion apportionment results with RIV MODFLOW output")
}

## load stream shapefile
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# figure out which segments are part of Navarro
segs.navarro <- shp.streams@data$SegNum[shp.streams@data$TerminalPa==outlet.TerminalPa]

## load geometric data (from script Navarro_DepletionApportionment.R)
if (masked){
  df.apportion <- 
    read.csv(file.path("results","Navarro_DepletionApportionment_MaskDryStreams_AllMethods+Wells+Reaches.csv"),
             stringsAsFactors=F) %>% 
    subset(SegNum %in% segs.navarro)
  
} else {
  df.apportion <- 
    read.csv(file.path("results","Navarro_DepletionApportionment_AllMethods+Wells+Reaches.csv"),
             stringsAsFactors=F) %>% 
    subset(SegNum %in% segs.navarro)
}

## load MODFLOW output
# leakage
df.MODFLOW <- 
  file.path("modflow","HTC", "Navarro", "SteadyState", stream_BC, modflow_v, paste0(stream_BC, "-SummarizeLeakage.csv")) %>% 
  read.csv(stringsAsFactors=F)

## if it's SFR: need to get SegNum
if (stream_BC=="SFR"){
  df.MODFLOW <- 
    read.table(file.path("modflow", "input", "isfr_SegmentData.txt"), 
               stringsAsFactors=F, header=T) %>% 
    left_join(df.MODFLOW, ., by="SFR_NSEG")
  
}

## calculate MODFLOW depletion
# make separate column for no pumping depletion
df.MODFLOW<- 
  df.MODFLOW %>% 
  subset(WellNum==0, select=c("SegNum", "leakage", "MNW_net")) %>% 
  set_colnames(c("SegNum", "leakage_NoPump", "MNW_NoPump")) %>% 
  left_join(subset(df.MODFLOW, WellNum != 0), ., by=c("SegNum")) %>% 
  arrange(WellNum, SegNum) %>% 
  subset(SegNum %in% segs.navarro)

# calculate change in leakage and MNW
# convention: negative value means flow out of aquifer (gaining stream, pumping well)
# so if (leakage_NoPump - leakage) < 0, stream is depleted due to pumping
df.MODFLOW$Qw_m3.d <- df.MODFLOW$MNW_net - df.MODFLOW$MNW_NoPump
df.MODFLOW$depletion_m3.d <- df.MODFLOW$leakage_NoPump - df.MODFLOW$leakage

## two approaches to calculating modflow depletion:
# (1) calculate depletion.prc as percentage of pumping rate (this is the definition of capture fraction from Leake et al., 2010)
df.MODFLOW$depletion.prc.modflow <- df.MODFLOW$depletion_m3.d/df.MODFLOW$Qw_m3.d

# # (2) calculate depletion.prc as percentage of all changes in leakage (this forces range 0-1)
# df.MODFLOW$depletion.prc.modflow <- df.MODFLOW$depletion_m3.d/df.MODFLOW$depletion.sum
# 
# ggplot(df.MODFLOW, aes(x=depletion.prc.modflow.1, y=depletion.prc.modflow.2)) + 
#   geom_abline(intercept=0, slope=1, color="gray65") + 
#   geom_point() + 
#   stat_smooth(method="lm")

## combine MODFLOW with analytical
df <- 
  left_join(df.apportion, df.MODFLOW[,c("SegNum", "WellNum", "depletion.prc.modflow")], by=c("SegNum", "WellNum"), all.x=T, all.y=T) %>% 
  subset(is.finite(depletion.prc.modflow)) %>% 
  melt(id=c("SegNum", "WellNum", "distToWell.min.m", "depletion.prc.modflow"),
       value.name="depletion.prc", variable.name="method")

## make plots
# figure out limits
min.depletion <- min(c(min(df$depletion.prc.modflow), min(df$depletion.prc), 0))
max.depletion <- max(c(max(df$depletion.prc.modflow), max(df$depletion.prc), 1))

p.depletion <-
  ggplot(df, aes(x=depletion.prc, y=depletion.prc.modflow)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point(shape=21) +
  stat_smooth(method="lm") +
  facet_wrap(~method, labeller=as_labeller(labels.method)) +
  labs(title="Comparison of Depletion Apportionment Equations to Steady-State MODFLOW",
       subtitle=paste0("MODFLOW: ", modflow_v, "; Stream BC: ", stream_BC, "; Mask: ", masked, "; 1 dot = 1 stream segment response to 1 pumping well")) +
  scale_x_continuous(name="Analytical Depletion Fraction", 
                     limits=c(0,1), breaks=seq(0,1,0.25), expand=c(0,0)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", 
                     limits=c(min.depletion,max.depletion), expand=c(0,0))
#p.depletion
ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_p.depletion_", stream_BC, "_", modflow_v, "_mask", masked, ".png")),
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

df.fit[,c("method", "MSE.bias.norm", "MSE.var.norm", "MSE.cor.norm")]
