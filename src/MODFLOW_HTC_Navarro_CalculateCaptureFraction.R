## MODFLOW_HTC_Navarro_CalculateCaptureFraction.R
#' This script is intended to calculate capture fraction for MODFLOW output.
#' It needs to be run after the scripts MODFLOW_HTC_Navarro_***-SummarizeLeakage.py
#' 
#' Note that the output of MODFLOW_HTC_Navarro_***-SummarizeLeakage.py is too big
#' for GitHub so it is needs to be transferred from the ComputeCanada HTC manually.

source(file.path("src", "paths+packages.R"))

## choose stream boundary condition and modflow version
stream_BC <- "RIV"       # "RIV" or "SFR"
modflow_v <- "mfnwt"     # "mfnwt" or "mf2005"
timeType  <- "Transient" # "SteadyState" or "Transient"

## define which directory you are interested in
dir.runs <- file.path("modflow", "HTC", "Navarro", timeType, stream_BC, modflow_v)

## load MODFLOW output
# leakage
df.MODFLOW <- 
  file.path(dir.runs, paste0(stream_BC, "-SummarizeLeakage.csv")) %>% 
  read.csv(stringsAsFactors=F)

## if it's SFR: need to get SegNum, figure out Time
if (stream_BC=="SFR"){
  # get SegNum
  df.MODFLOW <- 
    read.table(file.path("modflow", "input", "isfr_SegmentData.txt"), 
               stringsAsFactors=F, header=T) %>% 
    left_join(df.MODFLOW, ., by="SFR_NSEG") %>% 
    dplyr::select(Qaquifer, WellNum, kstpkper, MNW_net, SegNum) %>% 
    set_colnames(c("leakage", "WellNum", "kstpkper", "MNW_net", "SegNum"))
  
  # number of days in each stress period
  days_in_sp <- 
    seq(1,12) %>% 
    days_in_month() %>% 
    rep(.,2) %>% 
    as.numeric()
  
  # extract timestep and stress period from kstpkper
  df.sp <- 
    df.MODFLOW$kstpkper %>%
    unique() %>%
    str_replace(., "\\(", "") %>% 
    str_replace(., "\\)", "") %>% 
    str_split(., ",") %>% 
    unlist(.) %>% 
    as.numeric() %>% 
    "+"(1) %>%          # python indexing is 0-based
    matrix(., nrow=length(unique(df.MODFLOW$kstpkper)), byrow=T) %>% 
    as.data.frame() %>% 
    cbind(unique(df.MODFLOW[,c("kstpkper")])) %>% 
    set_colnames(c("stp", "per", "kstpkper")) %>% 
    transform(Time = NaN)
  
  get_Time <- function(stp, per, days_in_sp){
    if (per==1){
      Time <- stp
    } else {
      Time <- sum(days_in_sp[1:(per-1)]) +stp
    }
  }
  
  for (i in 1:length(df.sp$Time)){
    df.sp$Time[i] <- get_Time(stp=df.sp$stp[i], per=df.sp$per[i], days_in_sp=days_in_sp)
  }
  
  df.MODFLOW <- 
    df.MODFLOW %>% 
    left_join(., df.sp[,c("kstpkper", "Time")], by="kstpkper") %>% 
    dplyr::select(Time, WellNum, SegNum, MNW_net, leakage)
  
}

# at this point: df.MODFLOW should have the columns:
#  Time, WellNum, SegNum, MNW_net, leakage

## calculate MODFLOW depletion
# make separate column for no pumping depletion
df.MODFLOW <- 
  df.MODFLOW %>% 
  subset(WellNum==0, select=c("SegNum", "leakage", "MNW_net", "Time")) %>% 
  set_colnames(c("SegNum", "leakage_NoPump", "MNW_NoPump", "Time")) %>% 
  left_join(subset(df.MODFLOW, WellNum != 0), ., by=c("SegNum", "Time")) %>% 
  arrange(Time, WellNum, SegNum)

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

## check if depletion.prc.modflow sums to 0
df.MODFLOW.sum <-
  df.MODFLOW %>% 
  subset((Time >= sum(days_in_month(seq(1,4))) + 1)) %>% 
  group_by(Time, WellNum) %>% 
  summarize(depletion.prc.total = sum(depletion.prc.modflow))

## save output
# only save time/well/segment combos where depletion occurs
f.thres <- 0.0001  # =0.01%
df.MODFLOW[,c("Qw_m3.d", "depletion_m3.d")] <- signif(df.MODFLOW[,c("Qw_m3.d", "depletion_m3.d")], 4)
df.MODFLOW %>% 
  subset(is.finite(depletion.prc.modflow) & abs(depletion.prc.modflow) > f.thres) %>% 
  dplyr::select(SegNum, WellNum, Time, Qw_m3.d, depletion_m3.d) %>% 
  write.csv(., file.path(dir.runs, "Depletion_MODFLOW.csv"), row.names=F)

p.depletion.sum <-
  ggplot(df.MODFLOW.sum, aes(x=Time, y=depletion.prc.total, group=WellNum)) +
  geom_line(alpha=0.25) +
  coord_cartesian(ylim=c(-1,1))
