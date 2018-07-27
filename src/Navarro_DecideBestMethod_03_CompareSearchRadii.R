## Navarro_DecideBestMethod_03_CompareSearchRadii.R
#' This script is intended to decide the best search radii based on one criterion:
#'   (3) Most accurate prediction of depletion in all stream reaches --> lowest MSE, highest KGE

source(file.path("src", "paths+packages.R"))

#### (0) Prep various inputs, parameters, etc
## which depletion apportionment equation(s) to compare? winners from Navarro_DecideBestMethod_01_CompareDepletionApportionment.R
methods.plot <- c("Qf.Web", "Qf.WebSq")  # options: c("Qf.InvDist", "Qf.InvDistSq", "Qf.Web", "Qf.WebSq", "Qf.TPoly")

## which analytical model to use? winner from Navarro_DecideBestMethod_02_CompareAnalytical.R
analytical.plot <- "hunt"  # options: "hunt" or "glover"

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

# which MODFLOW version? 
modflow_v <- "mfnwt"  # (no other option)

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

# (3) Most accurate prediction of depletion in all reaches -----------------------------

### only have to run this once:
# start.flag <- T
# for (timeType in c("Transient", "Intermittent")) {
# 
#   ## open MODFLOW results
#   df.MODFLOW.RIV <-
#     file.path("modflow","HTC", "Navarro", timeType, "RIV", modflow_v, "Depletion_MODFLOW.csv") %>%
#     read.csv(stringsAsFactors=F) %>%
#     #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
#     transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
#               stream_BC = "RIV",
#               stringsAsFactors=F) 
# 
#   df.MODFLOW.SFR <-
#     file.path("modflow","HTC", "Navarro", timeType, "SFR", modflow_v, "Depletion_MODFLOW.csv") %>%
#     read.csv(stringsAsFactors=F) %>%
#     #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
#     transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
#               stream_BC = "SFR",
#               stringsAsFactors=F)
# 
#   df.MODFLOW <- 
#     rbind(df.MODFLOW.RIV, df.MODFLOW.SFR) %>% 
#     subset(depletion.prc.modflow > f.thres)
#   df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
# 
#   for (apportionment_name in c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic", "Adjacent+Dynamic")) {
#     ## load analytical output
#     df.analytical <-
#       paste0("Depletion_Analytical_", timeType, "_", apportionment_name, "_AllMethods+Wells+Reaches.csv") %>%
#       file.path("results", .) %>%
#       read.csv(stringsAsFactors=F) %>%
#       subset(analytical %in% analytical.plot) %>% 
#       dplyr::select(c("SegNum", "WellNum", "Time", "analytical", methods.plot)) %>%
#       melt(id=c("SegNum", "WellNum", "Time", "analytical"),
#            value.name="depletion.prc", variable.name="method") %>% 
#       subset(depletion.prc > f.thres)
# 
#     ## Time has long decimals; round before merging to ensure time match
#     df.analytical$Time <- round(df.analytical$Time, 1)
# 
#     for (BC in c("SFR", "RIV")) {
#       for (m in unique(df.analytical$method)) {
#         for (a in unique(df.analytical$analytical)) {
#           # combine and calculate fit
#           df <-
#             df.MODFLOW %>%
#             subset(stream_BC == BC) %>%
#             left_join(subset(df.analytical, method==m & analytical==a),
#                       by=c("WellNum", "Time")) %>%
#             replace_na(list("analytical"=a, "method"=m, "stream_BC"=BC, "depletion.prc"=0, "depletion.prc.modflow" = 0)) %>%
#             transform(apportionment = apportionment_name,
#                       pump = timeType,
#                       stringsAsFactors=F) %>% 
#             group_by(stream_BC, pump, analytical, apportionment, method, Time) %>% 
#             summarize(MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
#                       KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))
#             
# 
#           if (start.flag) {
#             df.fit.all <- df
#             start.flag <- F
#           } else {
#             df.fit.all <- rbind(df.fit.all, df)
#           }
# 
#           # status update
#           print(paste(timeType, apportionment_name, BC, m, a, "complete"))
# 
#         }  # end of a loop
#       }  # end of m loop
#     }  # end of BC loop
#   }  # end of apportionment_name loop
# }  # end of timeType loop
# 
# ## save fit statistics!
# write.csv(df.fit.all, file.path("results", "Navarro_DecideBestMethod_03_CompareSearchRadii_fit-All.csv"),
#           row.names=F, quote=F)

## read in fit statistics
df.fit.all <- read.csv(file.path("results", "Navarro_DecideBestMethod_03_CompareSearchRadii_fit-All.csv"),
                       stringsAsFactors=F)
df.fit.all$apportionment <- factor(df.fit.all$apportionment, levels=c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic", "Adjacent+Dynamic"))

## plots
stream_BC_plot <- c("RIV")

p.MSE <- 
  df.fit.all %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=MSE.overall, color=apportionment)) +
  geom_line() +
  facet_wrap(pump+stream_BC+method ~ ., ncol=4) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="MSE, all wells and reaches") +
  theme(legend.position="bottom")

p.KGE <- 
  df.fit.all %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=apportionment)) +
  geom_line() +
  facet_wrap(pump+stream_BC+method ~ ., ncol=4) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE, all wells and reaches") +
  theme(legend.position="bottom")

ggsave(file.path("results", "Navarro_DecideBestMethod_03_CompareSearchRadii_fit-All.png"),
       grid.arrange(p.KGE, p.MSE, ncol=1),
       width=190, height=200, units="mm")
