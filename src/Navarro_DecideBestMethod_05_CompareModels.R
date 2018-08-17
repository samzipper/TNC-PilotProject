## Navarro_DecideBestMethod_05_CompareModels.R
#' This script is intended to compare the output fromNavarro_DecideBestMethod_04_TweakParameters.R
#' to MODFLOW results.

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

#### (0) Prep various inputs, parameters, etc
## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

# which MODFLOW version? 
modflow_v <- "mfnwt"  # (no other option)

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

# Calculate transient fit statistics -----------------------------

# ## only have to run this once
# start.flag <- T
# for (timeType in c("Intermittent", "Transient")) {
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
#   #df.MODFLOW.SFR <-
#   #  file.path("modflow","HTC", "Navarro", timeType, "SFR", modflow_v, "Depletion_MODFLOW.csv") %>%
#   #  read.csv(stringsAsFactors=F) %>%
#   #  #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
#   #  transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
#   #            stream_BC = "SFR",
#   #            stringsAsFactors=F)
#   
#   #df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)
#   df.MODFLOW <- df.MODFLOW.RIV
#   df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
#   
#   ## for each stream_BC, WellNum, and Time find the most affected segment
#   df.MODFLOW.most <-
#     df.MODFLOW %>%
#     subset(depletion.prc.modflow > f.thres) %>%
#     #group_by(stream_BC, WellNum, Time) %>%
#     group_by(WellNum, Time) %>%
#     filter(depletion.prc.modflow==max(depletion.prc.modflow))
#   
#   ## Load output from Navarro_DecideBestMethod_04_TweakParameters.R
#   df.analytical <- read.csv(file.path("results", paste0("Navarro_DecideBestMethod_04_TweakParameters-", timeType, ".csv")),
#                             stringsAsFactors=F)
#   df.analytical$Time <- round(df.analytical$Time, 1)
#   colnames(df.analytical)[colnames(df.analytical)=="Qf"] <- "depletion.prc"
#   
#   ## for each stream_BC, WellNum, and Time find the most affected segment
#   df.analytical.max <-
#     df.analytical %>%
#     subset(depletion.prc > f.thres) %>%
#     group_by(Qf.thres, web.exp, WellNum, Time) %>%
#     filter(depletion.prc == max(depletion.prc)) %>%
#     dplyr::select(Time, WellNum, SegNum, Qf.thres, web.exp)
#   
#   for (qt in unique(df.analytical$Qf.thres)) {
#     for (w in unique(df.analytical$web.exp)) {
#       # combine max depletion segment info
#       df.max <-
#         df.MODFLOW.most %>%
#         dplyr::select(SegNum, WellNum, Time, depletion.prc.modflow) %>%
#         left_join(subset(df.analytical.max, web.exp==w & Qf.thres==qt),
#                   by=c("WellNum", "Time"), suffix=c(".modflow", ".analytical")) %>%
#         # add depletion in the most affected MODFLOW segment
#         left_join(subset(df.analytical,
#                          web.exp==w & Qf.thres==qt)[,c("WellNum", "Time", "web.exp", "Qf.thres", "SegNum", "depletion.prc")],
#                   by=c("WellNum", "Time", "web.exp", "Qf.thres", "SegNum.modflow"="SegNum")) %>%
#         replace_na(list("SegNum.analytical"=9999, "Qf.thres"=qt, "web.exp"=w, "depletion.prc" = 0)) %>%
#         transform(pump = timeType,
#                   stringsAsFactors=F)
#       
#       # combine and calculate fit
#       df <-
#         df.MODFLOW %>%
#         dplyr::select(WellNum, SegNum, Time, depletion.prc.modflow) %>%
#         left_join(subset(df.analytical,
#                          web.exp==w & Qf.thres==qt)[,c("WellNum", "Time", "SegNum", "Qf.thres", "web.exp", "depletion.prc")],
#                   by=c("WellNum", "Time", "SegNum")) %>%
#         replace_na(list("web.exp"=w, "Qf.thres"=qt,  "depletion.prc"=0, "depletion.prc.modflow" = 0)) %>%
#         group_by(Time, Qf.thres, web.exp) %>%
#         summarize(bias.overall = pbias(depletion.prc, depletion.prc.modflow),
#                   MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
#                   KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
#         transform(pump = timeType,
#                   stringsAsFactors=F)
#       
#       
#       if (start.flag) {
#         df.max.all <- df.max
#         df.fit.all <- df
#         start.flag <- F
#       } else {
#         df.max.all <- rbind(df.max.all, df.max)
#         df.fit.all <- rbind(df.fit.all, df)
#       }
#       
#       # status update
#       print(paste(timeType, w, qt, "complete"))
#       
#     }  # end of w loop
#   }  # end of qt loop
# }  # end of timeType loop
# 
# ## save fit statistics for all reach-well combos!
# write.csv(df.fit.all, file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-All.csv"),
#           row.names=F, quote=F)
# 
# ## calculate fit statistics
# df.fit.match <-
#   df.max.all %>%
#   group_by(Time, pump, Qf.thres, web.exp) %>%
#   summarize(n.reach = sum(is.finite(SegNum.modflow)),
#             n.match = sum(SegNum.modflow==SegNum.analytical),
#             bias.match = pbias(depletion.prc, depletion.prc.modflow),
#             MSE.match = MSE(depletion.prc, depletion.prc.modflow),
#             KGE.match = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
#   transform(prc.match = n.match/n.reach)
# 
# ## save fit statistics!
# write.csv(df.fit.match, file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match.csv"),
#           row.names=F, quote=F)

## load data and plot
df.fit.match <- read.csv(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match.csv"), 
                         stringsAsFactors=F)

df.fit.all <- read.csv(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-All.csv"), 
                       stringsAsFactors=F)

## plots - compare web.exp for a single Qf.thres
# best web.exp based on KGE at most affected reach
p.match.KGE <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot(aes(x=Time, y=KGE.match, color=web.exp, group=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE at Most Affected Reach") +
  scale_color_continuous(name="Web Exponent") +
  theme(legend.position="bottom") +
  NULL

p.match.KGE.best <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  group_by(Time, pump) %>% 
  filter(KGE.match == max(KGE.match)) %>% 
  ggplot(aes(x=Time, y=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Best Web Exponent based on\nKGE of Most Affected Reach")

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match_KGE.png"),
       grid.arrange(p.match.KGE, p.match.KGE.best, ncol=1),
       width=6, height=8, units="in")

# best web.exp based on bias at most affected reach
p.match.bias <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot(aes(x=Time, y=bias.match, color=web.exp, group=web.exp)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Bias at Most Affected Reach") +
  scale_color_continuous(name="Web Exponent") +
  theme(legend.position="bottom") +
  NULL

p.match.bias.best <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  group_by(Time, pump) %>% 
  filter(bias.match == min(abs(bias.match))) %>% 
  ggplot(aes(x=Time, y=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Best Web Exponent based on\nBias of Most Affected Reach")

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match_Bias.png"),
       grid.arrange(p.match.bias, p.match.bias.best, ncol=1),
       width=6, height=8, units="in")

# best web.exp based on MSE at most affected reach
p.match.MSE <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot(aes(x=Time, y=MSE.match, color=web.exp, group=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="MSE at Most Affected Reach") +
  scale_color_continuous(name="Web Exponent") +
  theme(legend.position="bottom") +
  NULL

p.match.MSE.best <-
  df.fit.match %>% 
  subset(Qf.thres == 0.01) %>% 
  group_by(Time, pump) %>% 
  filter(MSE.match == min(MSE.match)) %>% 
  ggplot(aes(x=Time, y=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Best Web Exponent based on\nMSE of Most Affected Reach")

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match_MSE.png"),
       grid.arrange(p.match.MSE, p.match.MSE.best, ncol=1),
       width=6, height=8, units="in")

# best web.exp based on KGE at all reaches
p.all.KGE <-
  df.fit.all %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=web.exp, group=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE, all Reaches") +
  scale_color_continuous(name="Web Exponent") +
  theme(legend.position="bottom") +
  NULL

p.all.KGE.best <-
  df.fit.all %>% 
  subset(Qf.thres == 0.01) %>% 
  group_by(Time, pump) %>% 
  filter(KGE.overall == max(KGE.overall)) %>% 
  ggplot(aes(x=Time, y=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Best Web Exponent based on\nKGE of All Reaches")

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-All_KGE.png"),
       grid.arrange(p.all.KGE, p.all.KGE.best, ncol=1),
       width=6, height=8, units="in")

# best web.exp based on bias at all reaches
p.all.bias <-
  df.fit.all %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot(aes(x=Time, y=bias.overall, color=web.exp, group=web.exp)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Bias at Most Affected Reach") +
  scale_color_continuous(name="Web Exponent") +
  theme(legend.position="bottom") +
  NULL

p.all.bias.best <-
  df.fit.all %>% 
  subset(Qf.thres == 0.01) %>% 
  group_by(Time, pump) %>% 
  filter(bias.overall == min(abs(bias.overall))) %>% 
  ggplot(aes(x=Time, y=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Best Web Exponent based on\nBias of Most Affected Reach")

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-All_Bias.png"),
       grid.arrange(p.all.bias, p.all.bias.best, ncol=1),
       width=6, height=8, units="in")

# best web.exp based on MSE at all reaches
p.all.MSE <-
  df.fit.all %>% 
  subset(Qf.thres == 0.01) %>% 
  ggplot(aes(x=Time, y=MSE.overall, color=web.exp, group=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="MSE, all Reaches") +
  scale_color_continuous(name="Web Exponent") +
  theme(legend.position="bottom") +
  NULL

p.all.MSE.best <-
  df.fit.all %>% 
  subset(Qf.thres == 0.01) %>% 
  group_by(Time, pump) %>% 
  filter(MSE.overall == min(MSE.overall)) %>% 
  ggplot(aes(x=Time, y=web.exp)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Best Web Exponent based on\nMSE of All Reaches")

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-All_MSE.png"),
       grid.arrange(p.all.MSE, p.all.MSE.best, ncol=1),
       width=6, height=8, units="in")

## now: for a given web.exp, choose best Qf.thres
# best web.exp based on KGE at all reaches
p.all.KGE.Qf <-
  df.fit.all %>% 
  subset(web.exp == 1.75) %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=factor(Qf.thres*100))) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE, all Reaches") +
  scale_color_discrete(name="Dynamic Depletion Threshold [%]") +
  theme(legend.position="bottom") +
  NULL

p.all.KGE.best.Qf <-
  df.fit.all %>% 
  subset(web.exp == 1.75) %>% 
  group_by(Time, pump) %>% 
  filter(KGE.overall == max(KGE.overall)) %>% 
  ggplot(aes(x=Time, y=Qf.thres*100)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_log10(name="Best Dynamic Depletion Threshold [%]\nbased on KGE of All Reaches") +
  NULL

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-All_Qf-KGE.png"),
       grid.arrange(p.all.KGE.Qf, p.all.KGE.best.Qf, ncol=1),
       width=6, height=8, units="in")

p.match.KGE.Qf <-
  df.fit.match %>% 
  subset(web.exp == 1.75) %>% 
  ggplot(aes(x=Time, y=KGE.match, color=factor(Qf.thres*100))) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE, Most Affected Reach") +
  scale_color_discrete(name="Dynamic Depletion Threshold [%]") +
  theme(legend.position="bottom") +
  NULL

p.match.KGE.best.Qf <-
  df.fit.match %>% 
  subset(web.exp == 1.75) %>% 
  group_by(Time, pump) %>% 
  filter(KGE.match == max(KGE.match)) %>% 
  ggplot(aes(x=Time, y=Qf.thres*100)) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_log10(name="Best Dynamic Depletion Threshold [%]\nbased on KGE of Most Affected Reach") +
  NULL

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match_Qf-KGE.png"),
       grid.arrange(p.match.KGE.Qf, p.match.KGE.best.Qf, ncol=1),
       width=6, height=8, units="in")

p.match.bias.Qf <-
  df.fit.match %>% 
  subset(web.exp == 1.75) %>% 
  ggplot(aes(x=Time, y=bias.match, color=factor(Qf.thres*100))) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  facet_grid(. ~ pump, scales="free_y") +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="Bias, Most Affected Reach") +
  scale_color_discrete(name="Dynamic Depletion Threshold [%]") +
  theme(legend.position="bottom") +
  NULL

p.match.bias.best.Qf <-
  df.fit.match %>% 
  subset(web.exp == 1.75) %>% 
  group_by(Time, pump) %>% 
  filter(bias.match == min(abs(bias.match))) %>% 
  ggplot(aes(x=Time, y=Qf.thres*100)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  facet_grid(. ~ pump) +
  scale_x_continuous(name="Time [days]") +
  scale_y_log10(name="Best Dynamic Depletion Threshold [%]\nbased on Bias of Most Affected Reach") +
  NULL

ggsave(file.path("results", "Navarro_DecideBestMethod_05_CompareModels_fit-Match_Qf-Bias.png"),
       grid.arrange(p.match.bias.Qf, p.match.bias.best.Qf, ncol=1),
       width=6, height=8, units="in")
