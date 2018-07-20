## Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState.R
#' This script is intended to compare estimated steady-state streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_RIV-SummarizeLeakage.py)
#' and geometric depletion apportionment methods.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

#### First: process steady-state MODFLOW data
## MODFLOW - steady state results
stream_BC <- "RIV"
modflow_v <- "mfnwt"

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## threshold for analysis
f.thres <- 0.001

## load stream shapefile
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# figure out which segments are part of Navarro
segs.navarro <- shp.streams@data$SegNum[shp.streams@data$TerminalPa==outlet.TerminalPa]

# leakage
df_MODFLOW <- 
  file.path("modflow","HTC", "Navarro", "SteadyState", stream_BC, modflow_v, paste0(stream_BC, "-SummarizeLeakage.csv")) %>% 
  read.csv(stringsAsFactors=F)

## if it's SFR: need to get SegNum
if (stream_BC=="SFR"){
  df_MODFLOW <- 
    read.table(file.path("modflow", "input", "isfr_SegmentData.txt"), 
               stringsAsFactors=F, header=T) %>% 
    left_join(df_MODFLOW, ., by="SFR_NSEG")
  
}

## calculate MODFLOW depletion
# make separate column for no pumping depletion
df_MODFLOW<- 
  df_MODFLOW %>% 
  subset(WellNum==0, select=c("SegNum", "leakage", "MNW_net")) %>% 
  set_colnames(c("SegNum", "leakage_NoPump", "MNW_NoPump")) %>% 
  left_join(subset(df_MODFLOW, WellNum != 0), ., by=c("SegNum")) %>% 
  arrange(WellNum, SegNum) %>% 
  subset(SegNum %in% segs.navarro)

# calculate change in leakage and MNW
# convention: negative value means flow out of aquifer (gaining stream, pumping well)
# so if (leakage_NoPump - leakage) < 0, stream is depleted due to pumping
df_MODFLOW$Qw_m3.d <- df_MODFLOW$MNW_net - df_MODFLOW$MNW_NoPump
df_MODFLOW$depletion_m3.d <- df_MODFLOW$leakage_NoPump - df_MODFLOW$leakage

## calculate MODFLOW depletion
#df_MODFLOW$depletion.prc.modflow <- df_MODFLOW$depletion_m3.d/df_MODFLOW$Qw_m3.d
df_MODFLOW$depletion.prc.modflow <- df_MODFLOW$depletion_m3.d/Qw

df_MODFLOW <- subset(df_MODFLOW, depletion.prc.modflow > f.thres)

#### load depletion apportionment output
# AdjacentOnly
df_adj <- 
  file.path("results", "Navarro_DepletionApportionment_AdjacentOnly_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  dplyr::select(WellNum, SegNum, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "depletion.prc.modflow")],
            by=c("SegNum", "WellNum")) %>% 
  transform(apportionment = "AdjacentOnly")
df_adj[is.na(df_adj)] <- 0

# LocalArea
df_loc <- 
  file.path("results", "Navarro_DepletionApportionment_LocalArea_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  dplyr::select(WellNum, SegNum, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "depletion.prc.modflow")],
            by=c("SegNum", "WellNum")) %>% 
  transform(apportionment = "LocalArea")
df_loc[is.na(df_loc)] <- 0

# Dynamic - use 'Transient' glover model at end of simulation
df_dyn <- 
  file.path("results", "Depletion_Analytical_Transient_Dynamic_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(Time==max(Time) & analytical=="glover") %>% 
  dplyr::select(WellNum, SegNum, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "depletion.prc.modflow")],
            by=c("SegNum", "WellNum")) %>% 
  transform(apportionment = "Dynamic")
df_dyn[is.na(df_dyn)] <- 0

# WholeDomain
df_whl <- 
  file.path("results", "Navarro_DepletionApportionment_WholeDomain_AllMethods+Wells+Reaches.csv") %>% 
  read.csv(stringsAsFactors = F) %>% 
  dplyr::select(WellNum, SegNum, f.InvDistSq, f.WebSq, f.TPoly) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "depletion.prc.modflow")],
            by=c("SegNum", "WellNum")) %>% 
  transform(apportionment = "WholeDomain")
df_whl[is.na(df_whl)] <- 0

## combine into single melted data frame
df <- 
  rbind(df_adj, df_dyn, df_loc, df_whl) %>% 
  melt(id=c("WellNum", "SegNum", "depletion.prc.modflow", "apportionment"),
       value.name="depletion.prc", variable.name="method")

# check for NAs
sum(is.na(df))

# get rid of reaches with no depletion
df <- subset(df, depletion.prc.modflow > f.thres | depletion.prc > f.thres)

# check to make sure TPoly works: these should be the same
sum(df$method=="f.TPoly" & df$apportionment=="AdjacentOnly")
sum(df$method=="f.TPoly" & df$apportionment=="WholeDomain")

# calculate fit statistics
df.fit <- 
  df %>% 
  group_by(method, apportionment) %>% 
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

df.fit.big <- 
  df %>% 
  subset(depletion.prc > 0.05 | depletion.prc.modflow > 0.05) %>% 
  group_by(method, apportionment) %>% 
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

## make plots
# fit bar plot
df.fit %>% 
  ggplot(aes(x=method, y=KGE.overall)) +
  geom_bar(stat="identity") +
  facet_grid(~apportionment) +
  ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState_FitKGE_", stream_BC, "_", modflow_v, ".png")),
         width=8, height=8, units="in")

df.fit %>% 
  ggplot(aes(x=method, y=MSE.overall)) +
  geom_bar(stat="identity") +
  facet_grid(~apportionment) +
  ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState_FitMSE_", stream_BC, "_", modflow_v, ".png")),
         width=8, height=8, units="in")

df.fit.big %>% 
  ggplot(aes(x=method, y=KGE.overall)) +
  geom_bar(stat="identity") +
  facet_grid(~apportionment) +
  ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState_FitBigKGE_", stream_BC, "_", modflow_v, ".png")),
         width=8, height=8, units="in")

df.fit.big %>% 
  ggplot(aes(x=method, y=MSE.overall)) +
  geom_bar(stat="identity") +
  facet_grid(~apportionment) +
  ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState_FitBigMSE_", stream_BC, "_", modflow_v, ".png")),
         width=8, height=8, units="in")

# figure out limits
min.depletion <- min(c(min(df$depletion.prc.modflow), min(df$depletion.prc), 0))
max.depletion <- max(c(max(df$depletion.prc.modflow), max(df$depletion.prc), 1))

df %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point(shape=21) +
  stat_smooth(method="lm") +
  facet_grid(method~apportionment, labeller=as_labeller(c(labels.method, labels.apportionment))) +
  labs(title="Comparison of Depletion Apportionment Equations to Steady-State MODFLOW",
       subtitle=paste0("MODFLOW: ", modflow_v, "; Stream BC: ", stream_BC, "; 1 dot = 1 stream segment response to 1 pumping well")) +
  scale_x_continuous(name="Analytical Depletion Fraction", 
                     limits=c(0,1), breaks=seq(0,1,0.25), expand=c(0,0)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", 
                     limits=c(min.depletion,max.depletion), expand=c(0,0)) +
  ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState_FitScatters_", stream_BC, "_", modflow_v, ".png")),
         width=10, height=8, units="in")

df %>% 
  subset(depletion.prc > 0.05 | depletion.prc.modflow > 0.05) %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow)) +
  geom_abline(intercept=0, slope=1, color="gray65") +
  geom_point(shape=21) +
  stat_smooth(method="lm") +
  facet_grid(method~apportionment, labeller=as_labeller(c(labels.method, labels.apportionment))) +
  labs(title="Comparison of Depletion Apportionment Equations to Steady-State MODFLOW",
       subtitle=paste0("MODFLOW: ", modflow_v, "; Stream BC: ", stream_BC, "; 1 dot = 1 stream segment response to 1 pumping well")) +
  scale_x_continuous(name="Analytical Depletion Fraction", 
                     limits=c(0,1), breaks=seq(0,1,0.25), expand=c(0,0)) +
  scale_y_continuous(name="MODFLOW Depletion Fraction", 
                     limits=c(min.depletion,max.depletion), expand=c(0,0)) +
  ggsave(file.path("results", paste0("Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState_FitBigScatters_", stream_BC, "_", modflow_v, ".png")),
         width=10, height=8, units="in")

#### plots comparing among depletion apportionment (no MODFLOW)
df_TPoly <- 
  full_join(subset(df_adj[,c("WellNum", "SegNum", "f.TPoly")], f.TPoly != 0),
            subset(df_loc[,c("WellNum", "SegNum", "f.TPoly")], f.TPoly != 0),
            by=c("WellNum", "SegNum")) %>% 
  full_join(subset(df_dyn[,c("WellNum", "SegNum", "f.TPoly")], f.TPoly != 0),
            by=c("WellNum", "SegNum")) %>% 
  full_join(subset(df_whl[,c("WellNum", "SegNum", "f.TPoly")], f.TPoly != 0),
            by=c("WellNum", "SegNum")) %>% 
  set_colnames(c("WellNum", "SegNum", "AdjacentOnly", "LocalArea", "Dynamic", "WholeDomain")) %>% 
  subset(WellNum %in% df_dyn$WellNum)
df_TPoly[is.na(df_TPoly)] <- 0

ggplot(df_TPoly, aes(x=WholeDomain, y=LocalArea)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

ggplot(df_TPoly, aes(x=WholeDomain, y=AdjacentOnly)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

ggplot(df_TPoly, aes(x=WholeDomain, y=Dynamic)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

## compare web squared
df_WebSq <- 
  full_join(subset(df_adj[,c("WellNum", "SegNum", "f.WebSq")], f.WebSq != 0),
            subset(df_loc[,c("WellNum", "SegNum", "f.WebSq")], f.WebSq != 0),
            by=c("WellNum", "SegNum")) %>% 
  full_join(subset(df_dyn[,c("WellNum", "SegNum", "f.WebSq")], f.WebSq != 0),
            by=c("WellNum", "SegNum")) %>% 
  full_join(subset(df_whl[,c("WellNum", "SegNum", "f.WebSq")], f.WebSq != 0),
            by=c("WellNum", "SegNum")) %>% 
  set_colnames(c("WellNum", "SegNum", "AdjacentOnly", "LocalArea", "Dynamic", "WholeDomain")) %>% 
  subset(WellNum %in% df_dyn$WellNum)
df_WebSq[is.na(df_WebSq)] <- 0

ggplot(df_WebSq, aes(x=WholeDomain, y=LocalArea)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

ggplot(df_WebSq, aes(x=WholeDomain, y=AdjacentOnly)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")

ggplot(df_WebSq, aes(x=WholeDomain, y=Dynamic)) +
  geom_point() +
  geom_abline(intercept=0, slope=1, color="red")
