## Navarro_DecideBestMethod_02_CompareAnalytical.R
#' This script is intended to decide the best analytical model based on two criteria:
#'   (2a) Most accurate prediction of depletion in most-affected reach --> lowest MSE, highest KGE
#'   (2b) Most accurate prediction of cumulative domain-wide depletion --> lowest MSE, highest KGE

source(file.path("src", "paths+packages.R"))

#### (0) Prep various inputs, parameters, etc
## which depletion apportionment equation(s) to compare? winners from Navarro_DecideBestMethod_01_CompareDepletionApportionment.R
methods.plot <- c("Qf.Web", "Qf.WebSq")

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

# which MODFLOW version? 
modflow_v <- "mfnwt"  # (no other option)

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

# (2a) Most accurate prediction of depletion in most-affected reach -----------------------------

## fit statistics are already calcualted from Navarro_DecideBestMethod_01_CompareDepletionApportionment.R
df.fit.match <- 
  read.csv(file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-Match.csv"),
           stringsAsFactors=F) %>% 
  subset(method %in% methods.plot)
df.fit.match$apportionment <- factor(df.fit.match$apportionment, levels=c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic"))

## plots
stream_BC_plot <- c("RIV")

df.fit.match %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=MSE.match, color=method, linetype=analytical)) +
  geom_line() +
  facet_grid(pump+stream_BC+method ~ apportionment) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="MSE of Most-Affected Reach") +
  scale_color_manual(values=pal.method.Qf, labels=labels.method.Qf) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_02_CompareAnalytical_fit-Match_MSE.png"),
         width=190, height=200, units="mm")

df.fit.match %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=KGE.match, color=method, linetype=analytical)) +
  geom_line() +
  facet_grid(pump+stream_BC+method ~ apportionment) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE of Most-Affected Reach") +
  scale_color_manual(values=pal.method.Qf, labels=labels.method.Qf) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_02_CompareAnalytical_fit-Match_KGE.png"),
         width=190, height=200, units="mm")

# (2b) Most accurate prediction of cumulative domain-wide depletion -----------------------------

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
#               stringsAsFactors=F) %>% 
#     group_by(WellNum, Time, stream_BC) %>% 
#     summarize(Qf.total.modflow = sum(depletion.prc.modflow))
#   
#   df.MODFLOW.SFR <- 
#     file.path("modflow","HTC", "Navarro", timeType, "SFR", modflow_v, "Depletion_MODFLOW.csv") %>% 
#     read.csv(stringsAsFactors=F) %>% 
#     #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
#     transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
#               stream_BC = "SFR",
#               stringsAsFactors=F) %>% 
#     group_by(WellNum, Time, stream_BC) %>% 
#     summarize(Qf.total.modflow = sum(depletion.prc.modflow))
#   
#   df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)
#   df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
#   
#   for (apportionment_name in c("LocalArea", "AdjacentOnly", "WholeDomain", "Dynamic")) {
#     ## load analytical output
#     df.analytical <- 
#       paste0("Depletion_Analytical_", timeType, "_", apportionment_name, "_AllMethods+Wells+Reaches.csv") %>% 
#       file.path("results", .) %>% 
#       read.csv(stringsAsFactors=F) %>% 
#       dplyr::select(c("SegNum", "WellNum", "Time", "analytical", methods.plot)) %>% 
#       melt(id=c("SegNum", "WellNum", "Time", "analytical"),
#            value.name="depletion.prc", variable.name="method") %>% 
#       group_by(WellNum, Time, analytical, method) %>% 
#       summarize(Qf.total.analytical = sum(depletion.prc))
#     
#     ## Time has long decimals; round before merging to ensure time match
#     df.analytical$Time <- round(df.analytical$Time, 1)
#     
#     for (BC in c("SFR", "RIV")) {
#       for (m in unique(df.analytical$method)) {
#         for (a in unique(df.analytical$analytical)) {
#           # combine
#           df.sum <- 
#             df.MODFLOW %>% 
#             subset(stream_BC == BC) %>% 
#             left_join(subset(df.analytical, method==m & analytical==a), 
#                       by=c("WellNum", "Time")) %>% 
#             replace_na(list("analytical"=a, "method"=m, "stream_BC"=BC, "Qf.total.analytical"=0, "Qf.total.modflow" = 0)) %>% 
#             transform(apportionment = apportionment_name,
#                       pump = timeType,
#                       stringsAsFactors=F)
#           
#           if (start.flag) {
#             df.sum.all <- df.sum
#             start.flag <- F
#           } else {
#             df.sum.all <- rbind(df.sum.all, df.sum)
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
# ## calculate fit statistics
# df.fit.sum <-
#   df.sum.all %>% 
#   group_by(stream_BC, pump, analytical, apportionment, method, Time) %>% 
#   summarize(MSE.sum = MSE(Qf.total.analytical, Qf.total.modflow),
#             KGE.sum = KGE(Qf.total.analytical, Qf.total.modflow, method="2012"))
# 
# ## save fit statistics!
# write.csv(df.fit.sum, file.path("results", "Navarro_DecideBestMethod_02_CompareAnalytical_fit-Sum.csv"),
#           row.names=F, quote=F)

## read in fit statistics
df.fit.sum <- read.csv(file.path("results", "Navarro_DecideBestMethod_02_CompareAnalytical_fit-Sum.csv"),
                       stringsAsFactors=F)
df.fit.sum$apportionment <- factor(df.fit.sum$apportionment, levels=c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic"))

## plots
stream_BC_plot <- c("RIV")

df.fit.sum %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=MSE.sum, color=method, linetype=analytical)) +
  geom_line() +
  facet_grid(pump+stream_BC+method ~ apportionment) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="MSE of Total Capture Fraction") +
  scale_color_manual(values=pal.method.Qf, labels=labels.method.Qf) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_02_CompareAnalytical_fit-Sum_MSE.png"),
         width=190, height=200, units="mm")

df.fit.sum %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=KGE.sum, color=method, linetype=analytical)) +
  geom_line() +
  facet_grid(pump+stream_BC+method ~ apportionment) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE of Total Capture Fraction") +
  scale_color_manual(values=pal.method.Qf, labels=labels.method.Qf) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_02_CompareAnalytical_fit-Sum_KGE.png"),
         width=190, height=200, units="mm")
