## Navarro_CompareMODFLOWtoDepletionApportionment_Transient.R
#' This script is intended to compare estimated transient streamflow depletion
#' from MODFLOW (calculate from output of script MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
#' and geometric depletion apportionment methods (from script Navarro_Analytical_Transient.R).

source(file.path("src", "paths+packages.R"))

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## choose stream boundary condition and modflow version
stream_BC <- "RIV"    # "RIV" or "SFR"
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"
masked <- F           # use masked depletion apportionment? should only be T for SFR
if (stream_BC != "SFR" & masked){
  stop("Error: Using masked depletion apportionment results with RIV MODFLOW output")
}

## load MODFLOW results
df.MODFLOW <- 
  file.path("modflow","HTC", "Navarro", "Transient", stream_BC, modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d)

# if it's SFR: need to get SegNum
if (stream_BC=="SFR"){
  df.MODFLOW <- 
    read.table(file.path("modflow", "input", "isfr_SegmentData.txt"), 
               stringsAsFactors=F, header=T) %>% 
    left_join(df.MODFLOW, ., by="SFR_NSEG")
  
}

## load analytical results
df.analytical <- 
  file.path("modflow","HTC", "Navarro", "Transient", stream_BC, modflow_v, "Depletion_Analytical.csv") %>% 
  read.csv(stringsAsFactors=F)

## combine MODFLOW with analytical
df <- 
  full_join(df.analytical, df.MODFLOW[,c("SegNum", "WellNum", "Time", "depletion.prc.modflow")], by=c("SegNum", "WellNum", "Time")) %>% 
  melt(id=c("SegNum", "WellNum", "Time", "analytical", "depletion.prc.modflow"),
       value.name="depletion.prc", variable.name="method")

# missing values should be 0 (0s were filtered out in previous scripts)
df$depletion.prc.modflow[is.na(df$depletion.prc.modflow)] <- 0
df$analytical[is.na(df$depletion.prc)] <- "hunt"
tmp <- subset(df, is.na(depletion.prc))
tmp$analytical <- "glover"
df <- rbind(df, tmp)
df$depletion.prc[is.na(df$depletion.prc)] <- 0

# calculate sum of depletion.prc for each well/timestep/method
df.sum <-
  df %>% 
  group_by(Time, WellNum, analytical, method) %>% 
  summarize(Qf.MODFLOW = sum(depletion.prc.modflow),
            Qf.analytical = sum(depletion.prc))

## calculate fit statistics
df.fit <- 
  df %>% 
  group_by(analytical, method, Time) %>% 
  summarize(n.reach = sum(is.finite(depletion.prc)),
            #cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
            #bias = pbias(depletion.prc, depletion.prc.modflow),
            #R2 = R2(depletion.prc, depletion.prc.modflow),
            #MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
            #MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
            #MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
            #MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
            #MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
            #MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
            #MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
            KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>% 
  subset(n.reach > 3)

## make plots
# some lines/annotations
ts.pump.start <- sum(days_in_month(seq(1,4))) + 1

p.depletion.ts <-
  ggplot(df.fit, aes(x=Time, y=KGE.overall, color=method)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_vline(xintercept=ts.pump.start, linetype="dashed") +
  annotate("rect", 
           xmin=0, 
           xmax=sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  annotate("rect", 
           xmin=sum(days_in_month(seq(1,10)))+1, 
           xmax=365+sum(days_in_month(seq(1,4))), 
           ymin=-Inf, 
           ymax=Inf, 
           fill=col.gray, alpha=0.5) +
  geom_line() +
  facet_wrap(~analytical, ncol=1) +
  scale_x_continuous(name="Day of Simulation", expand=c(0,0)) +
  scale_y_continuous(name="Kling-Gupta Efficiency")

p.depletion.scatter <-
  subset(df, Time %in% seq(130,730,100)) %>% 
  ggplot(aes(x=depletion.prc, y=depletion.prc.modflow, color=method)) +
  geom_abline(intercept=0, slope=1, color=col.gray) +
  geom_point() +
  facet_grid(Time ~ analytical)
