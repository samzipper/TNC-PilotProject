## Navarro_DecideBestMethod_01_CompareDepletionApportionment.R
#' This script is intended to decide the best depletion apportionment equation
#' under all analytical models and search radii based on two criteria:
#'   (1a) Best fit under steady-state conditions -->  lowest MSE, highest KGE
#'   (1b) Best at identifying the most-affected reach --> match percentage

source(file.path("src", "paths+packages.R"))

#### (0) Prep various inputs, parameters, etc
## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

# which MODFLOW version? 
modflow_v <- "mfnwt"  # (no other option)

## threshold for inclusion in fit statistics, etc
f.thres <- 0.001  # 0.1%

# (1a) Best fit under steady-state conditions -----------------------------

### only have to run this once
# ## First: process steady-state MODFLOW data
# # RIV
# df_MODFLOW_SS_RIV <-
#   file.path("modflow","HTC", "Navarro", "SteadyState", "RIV", modflow_v, "RIV-SummarizeLeakage.csv") %>%
#   read.csv(stringsAsFactors=F) %>%
#   transform(stream_BC = "RIV")
# 
# # SFR
# df_MODFLOW_SS_SFR <-
#   file.path("modflow","HTC", "Navarro", "SteadyState", "SFR", modflow_v, "SFR-SummarizeLeakage.csv") %>%
#   read.csv(stringsAsFactors=F) %>%
#   transform(stream_BC = "SFR") %>%
#   left_join(read.table(file.path("modflow", "input", "isfr_SegmentData.txt"),
#                        stringsAsFactors=F, header=T),
#             by="SFR_NSEG") %>%
#   dplyr::select(colnames(df_MODFLOW_SS_RIV))
# 
# # calculate MODFLOW depletion
# df_MODFLOW_SS_RIV <-
#   df_MODFLOW_SS_RIV %>%
#   subset(WellNum==0, select=c("SegNum", "leakage", "MNW_net")) %>%
#   set_colnames(c("SegNum", "leakage_NoPump", "MNW_NoPump")) %>%
#   left_join(subset(df_MODFLOW_SS_RIV, WellNum != 0), ., by=c("SegNum")) %>%
#   arrange(WellNum, SegNum) %>%
#   transform(depletion.prc.modflow = (leakage_NoPump - leakage)/Qw) %>%
#   subset(depletion.prc.modflow > f.thres)
# 
# df_MODFLOW_SS_SFR <-
#   df_MODFLOW_SS_SFR %>%
#   subset(WellNum==0, select=c("SegNum", "leakage", "MNW_net")) %>%
#   set_colnames(c("SegNum", "leakage_NoPump", "MNW_NoPump")) %>%
#   left_join(subset(df_MODFLOW_SS_SFR, WellNum != 0), ., by=c("SegNum")) %>%
#   arrange(WellNum, SegNum) %>%
#   transform(depletion.prc.modflow = (leakage_NoPump - leakage)/Qw) %>%
#   subset(depletion.prc.modflow > f.thres)
# 
# # combine MODFLOW
# df_MODFLOW_SS <-
#   rbind(df_MODFLOW_SS_RIV, df_MODFLOW_SS_SFR) %>%
#   dplyr::select(stream_BC, WellNum, SegNum, depletion.prc.modflow)
# 
# ## load depletion apportionment output
# # AdjacentOnly
# df_adj <-
#   file.path("results", "Navarro_DepletionApportionment_AdjacentOnly_AllMethods+Wells+Reaches.csv") %>%
#   read.csv(stringsAsFactors = F) %>%
#   dplyr::select(WellNum, SegNum, f.InvDist, f.InvDistSq, f.Web, f.WebSq, f.TPoly) %>%
#   transform(apportionment = "AdjacentOnly",
#             stream_BC = "RIV",
#             stringsAsFactors=F) %>%
#   rbind(.,.)
# df_adj$stream_BC[duplicated(df_adj)] <- "SFR"
# df_adj <-
#   full_join(df_adj, df_MODFLOW_SS, by=c("WellNum", "SegNum", "stream_BC")) %>%
#   replace_na(list("apportionment"="AdjacentOnly",
#                   "depletion.prc.modflow"=0,
#                   "f.InvDist"=0,
#                   "f.InvDistSq"=0,
#                   "f.Web"=0,
#                   "f.WebSq"=0,
#                   "f.TPoly"=0))
# 
# # LocalArea
# df_loc <-
#   file.path("results", "Navarro_DepletionApportionment_LocalArea_AllMethods+Wells+Reaches.csv") %>%
#   read.csv(stringsAsFactors = F) %>%
#   dplyr::select(WellNum, SegNum, f.InvDist, f.InvDistSq, f.Web, f.WebSq, f.TPoly) %>%
#   transform(apportionment = "LocalArea",
#             stream_BC = "RIV",
#             stringsAsFactors=F) %>%
#   rbind(.,.)
# df_loc$stream_BC[duplicated(df_loc)] <- "SFR"
# df_loc <-
#   full_join(df_loc, df_MODFLOW_SS, by=c("WellNum", "SegNum", "stream_BC")) %>%
#   replace_na(list("apportionment"="LocalArea",
#                   "depletion.prc.modflow"=0,
#                   "f.InvDist"=0,
#                   "f.InvDistSq"=0,
#                   "f.Web"=0,
#                   "f.WebSq"=0,
#                   "f.TPoly"=0))
# 
# # WholeDomain
# df_whl <-
#   file.path("results", "Navarro_DepletionApportionment_WholeDomain_AllMethods+Wells+Reaches.csv") %>%
#   read.csv(stringsAsFactors = F) %>%
#   dplyr::select(WellNum, SegNum, f.InvDist, f.InvDistSq, f.Web, f.WebSq, f.TPoly) %>%
#   transform(apportionment = "WholeDomain",
#             stream_BC = "RIV",
#             stringsAsFactors=F) %>%
#   rbind(.,.)
# df_whl$stream_BC[duplicated(df_whl)] <- "SFR"
# df_whl <-
#   full_join(df_whl, df_MODFLOW_SS, by=c("WellNum", "SegNum", "stream_BC")) %>%
#   replace_na(list("apportionment"="WholeDomain",
#                   "depletion.prc.modflow"=0,
#                   "f.InvDist"=0,
#                   "f.InvDistSq"=0,
#                   "f.Web"=0,
#                   "f.WebSq"=0,
#                   "f.TPoly"=0))
# 
# ## combine into single melted data frame
# df <-
#   rbind(df_adj, df_loc, df_whl) %>%
#   melt(id=c("WellNum", "SegNum", "depletion.prc.modflow", "apportionment", "stream_BC"),
#        value.name="depletion.prc", variable.name="method")
# 
# # check for NAs
# sum(is.na(df))
# 
# # get rid of reaches with no depletion
# df <- subset(df, (depletion.prc.modflow > f.thres) | (depletion.prc > f.thres))
# df <- df[((df$depletion.prc.modflow > f.thres) | (df$depletion.prc > f.thres)), ]
# 
# # check to make sure TPoly works: these should be the same
# sum(df$method=="f.TPoly" & df$apportionment=="AdjacentOnly" & df$stream_BC=="RIV")
# sum(df$method=="f.TPoly" & df$apportionment=="WholeDomain" & df$stream_BC=="RIV")
# 
# ## calculate fit statistics
# df.fit.SS <-
#   df %>%
#   group_by(stream_BC, method, apportionment) %>%
#   summarize(n.reach = sum(is.finite(depletion.prc)),
#             cor = cor(depletion.prc, depletion.prc.modflow, method="pearson"),
#             bias = pbias(depletion.prc, depletion.prc.modflow),
#             R2 = R2(depletion.prc, depletion.prc.modflow),
#             MSE.bias = MSE.bias(depletion.prc, depletion.prc.modflow),
#             MSE.var = MSE.var(depletion.prc, depletion.prc.modflow),
#             MSE.cor = MSE.cor(depletion.prc, depletion.prc.modflow),
#             MSE.bias.norm = MSE.bias.norm(depletion.prc, depletion.prc.modflow),
#             MSE.var.norm = MSE.var.norm(depletion.prc, depletion.prc.modflow),
#             MSE.cor.norm = MSE.cor.norm(depletion.prc, depletion.prc.modflow),
#             MSE.overall = MSE(depletion.prc, depletion.prc.modflow),
#             KGE.overall = KGE(depletion.prc, depletion.prc.modflow, method="2012"))
# 
# ## save fit statistics!
# write.csv(df.fit.SS, file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-SS.csv"),
#           row.names=F, quote=F)

## read in fit statistics
df.fit.SS <- read.csv(file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-SS.csv"),
                      stringsAsFactors=F)
df.fit.SS$apportionment <- factor(df.fit.SS$apportionment, levels=c("AdjacentOnly", "LocalArea", "WholeDomain"))

# get best method for each stream_BC and apportionment
df.best.MSE.SS <-
  df.fit.SS %>% 
  group_by(stream_BC, apportionment) %>% 
  filter(MSE.overall==min(MSE.overall)) %>% 
  transform(metric="MSE")

df.best.KGE.SS <-
  df.fit.SS %>% 
  group_by(stream_BC, apportionment) %>% 
  filter(KGE.overall==max(KGE.overall)) %>% 
  transform(metric="KGE")

df.best.overall.MSE.SS <-
  df.fit.SS %>% 
  group_by(stream_BC) %>% 
  filter(MSE.overall==min(MSE.overall)) %>% 
  transform(metric="MSE")

df.best.overall.KGE.SS <-
  df.fit.SS %>% 
  group_by(stream_BC) %>% 
  filter(KGE.overall==max(KGE.overall)) %>% 
  transform(metric="KGE")

## make plots
stream_BC_plot <- c("RIV")

# fit bar plot
p.fit.KGE.SS <- 
  df.fit.SS %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=method, y=KGE.overall)) +
  geom_bar(stat="identity") +
  geom_bar(data=subset(df.best.KGE.SS, stream_BC %in% stream_BC_plot), 
           fill="red", stat="identity") +
  geom_bar(data=subset(df.best.overall.KGE.SS,stream_BC %in% stream_BC_plot),
           fill="blue", stat="identity") +
  geom_hline(yintercept=0, color=col.gray) +
  scale_x_discrete(name="Depletion Apportionment Equation") +
  scale_y_continuous(name="KGE") +
  facet_wrap(~apportionment) 

p.fit.MSE.SS <- 
  df.fit.SS %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=method, y=MSE.overall)) +
  geom_bar(stat="identity") +
  geom_bar(data=subset(df.best.MSE.SS, stream_BC %in% stream_BC_plot),
           fill="red", stat="identity") +
  geom_bar(data=subset(df.best.overall.MSE.SS, stream_BC %in% stream_BC_plot),
           fill="blue", stat="identity") +
  scale_x_discrete(name="Depletion Apportionment Equation") +
  scale_y_continuous(name="MSE") +
  facet_wrap(~apportionment) 

ggsave(file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-SS-facets.png"),
       grid.arrange(p.fit.KGE.SS, p.fit.MSE.SS, ncol=1),
       width=190, height=200, units="mm")

df.fit.SS %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=method, y=MSE.overall, fill=apportionment)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="Depletion Apportionment Equation") +
  scale_y_continuous(name="MSE")

df.fit.SS %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=method, y=KGE.overall, fill=apportionment)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="Depletion Apportionment Equation") +
  scale_y_continuous(name="KGE")

df.fit.SS %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  dplyr::select(method, apportionment, KGE.overall, MSE.overall) %>% 
  melt(id=c("method", "apportionment"),
       variable.name="FitStatistic") %>%
  ggplot(aes(x=method, y=value, fill=apportionment)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_discrete(name="Depletion Apportionment Equation") +
  facet_wrap(~FitStatistic, scales="free_y", ncol=1) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-SS.png"),
         width=190, height=120, units="mm")

# (1b) Best at identifying the most-affected reach (transient) -----------------------------

### only have to run this once
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
#   df.MODFLOW <- rbind(df.MODFLOW.RIV, df.MODFLOW.SFR)
#   df.MODFLOW$Time <- round(df.MODFLOW$Time, 1)
# 
#   ## for each stream_BC, WellNum, and Time find the most affected segment
#   df.MODFLOW.most <-
#     df.MODFLOW %>%
#     subset(depletion.prc.modflow > f.thres) %>%
#     group_by(stream_BC, WellNum, Time) %>%
#     filter(depletion.prc.modflow==max(depletion.prc.modflow))
# 
#   for (apportionment_name in c("LocalArea", "AdjacentOnly", "WholeDomain", "Dynamic")) {
#     ## load analytical output
#     df.analytical <-
#       paste0("Depletion_Analytical_", timeType, "_", apportionment_name, "_AllMethods+Wells+Reaches.csv") %>%
#       file.path("results", .) %>%
#       read.csv(stringsAsFactors=F) %>%
#       dplyr::select(SegNum, WellNum, Time, analytical, Qf.InvDist, Qf.InvDistSq, Qf.Web, Qf.WebSq, Qf.TPoly) %>%
#       melt(id=c("SegNum", "WellNum", "Time", "analytical"),
#            value.name="depletion.prc", variable.name="method") %>%
#       subset(depletion.prc > f.thres)
# 
#     ## Time has long decimals; round before merging to ensure time match
#     df.analytical$Time <- round(df.analytical$Time, 1)
# 
#     df.analytical.max <-
#       df.analytical %>%
#       group_by(analytical, WellNum, Time, method) %>%
#       filter(depletion.prc==max(depletion.prc)) %>%
#       dplyr::select(analytical, WellNum, Time, method, SegNum)
# 
#     for (BC in c("SFR", "RIV")) {
#       for (m in unique(df.analytical$method)) {
#         for (a in unique(df.analytical.max$analytical)) {
#           # combine
#           df.max <-
#             df.MODFLOW.most %>%
#             subset(stream_BC == BC) %>%
#             dplyr::select(SegNum, WellNum, Time, stream_BC, depletion.prc.modflow) %>%
#             left_join(subset(df.analytical.max, method==m & analytical==a),
#                       by=c("WellNum", "Time"), suffix=c(".modflow", ".analytical")) %>%
#             # add depletion in the most affected MODFLOW segment
#             left_join(subset(df.analytical, method==m & analytical==a),
#                       by=c("WellNum", "Time", "SegNum.modflow"="SegNum", "analytical", "method")) %>%
#             replace_na(list("analytical"=a, "method"=m, "SegNum.analytical"=9999, "depletion.prc" = 0)) %>%
#             transform(apportionment = apportionment_name,
#                       pump = timeType)
# 
#           if (start.flag) {
#             df.max.all <- df.max
#             start.flag <- F
#           } else {
#             df.max.all <- rbind(df.max.all, df.max)
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
# df.fit.match <-
#   df.max.all %>%
#   group_by(stream_BC, pump, analytical, apportionment, method, Time) %>%
#   summarize(n.reach = sum(is.finite(SegNum.modflow)),
#             n.match = sum(SegNum.modflow==SegNum.analytical),
#             MSE.match = MSE(depletion.prc, depletion.prc.modflow),
#             KGE.match = KGE(depletion.prc, depletion.prc.modflow, method="2012")) %>%
#   transform(prc.match = n.match/n.reach)
# 
# ## save fit statistics!
# write.csv(df.fit.match, file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-Match.csv"),
#           row.names=F, quote=F)

## read in fit statistics
df.fit.match <- 
  read.csv(file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-Match.csv"),
           stringsAsFactors=F)
df.fit.match$apportionment <- factor(df.fit.match$apportionment, levels=c("AdjacentOnly", "LocalArea", "WholeDomain", "Dynamic"))

## plots
stream_BC_plot <- c("RIV")

df.fit.match %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=prc.match, color=method)) +
  geom_line() +
  facet_grid(pump+stream_BC+analytical ~ apportionment) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="% of Wells where Most-Affected Reach is Correctly Identified", 
                     limits=c(0,1), expand=c(0,0), breaks=seq(0,1,0.25)) +
  scale_color_manual(values=pal.method.Qf, labels=labels.method.Qf) +
  ggsave(file.path("results", "Navarro_DecideBestMethod_01_CompareDepletionApportionment_fit-Match_prc.png"),
         width=190, height=200, units="mm")

df.fit.match %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=MSE.match, color=method)) +
  geom_line() +
  facet_grid(pump+stream_BC+analytical ~ apportionment) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="MSE of Most-Affected Reach") +
  scale_color_manual(values=pal.method.Qf, labels=labels.method.Qf)

df.fit.match %>% 
  subset(stream_BC %in% stream_BC_plot) %>% 
  ggplot(aes(x=Time, y=KGE.match, color=method)) +
  geom_line() +
  facet_grid(pump+stream_BC+analytical ~ apportionment) +
  scale_x_continuous(name="Time [days]") +
  scale_y_continuous(name="KGE of Most-Affected Reach") +
  scale_color_manual(values=pal.method.Qf, labels=labels.method.Qf)
