## Figure_MODFLOWSpinUp.R
#' This figure plots net river leakage for the entire domain during the spin-up period.

source(file.path("src", "paths+packages.R"))

## load RIV and SFR output
df.RIV <- 
  file.path("modflow", "Navarro-Transient-SpinUp", "RIV", "mfnwt", "Navarro-Transient-SpinUp_SummarizeBudget.csv") %>% 
  read.csv() %>% 
  transform(stream_BC = "RIV")

df.SFR <- 
  file.path("modflow", "Navarro-Transient-SpinUp", "SFR", "mfnwt", "Navarro-Transient-SpinUp_SummarizeBudget.csv") %>% 
  read.csv() %>% 
  transform(stream_BC = "SFR")

df <- rbind(df.RIV[,c("datetime", "leakage", "stream_BC")], 
            df.SFR[,c("datetime", "leakage", "stream_BC")])
df$datetime <- ymd_hms(df$datetime)

# calculate days since start of simulation
df$days_since_start <- (df$datetime - min(df$datetime))/86400
df$year <- year(df$datetime)

# summarize to annual
df.ann <- 
  df %>% 
  group_by(stream_BC, year) %>% 
  summarize(leakage_max = max(leakage),
            leakage_min = min(leakage),
            leakage_range = leakage_max - leakage_min) %>% 
  transform(change_range = c(NaN, diff(leakage_range))) %>% 
  transform(change_range_prc = change_range/leakage_range)
df.ann[min(which(df.ann$stream_BC=="SFR")), c("change_range", "change_range_prc")] <- NaN

# monthly plot
p.mo.riv <- 
  ggplot(subset(df, stream_BC=="RIV"), aes(x=days_since_start/365, y=leakage)) +
  geom_line(color=col.cat.blu) +
  scale_x_continuous(name=NULL, 
                     limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0)) +
  scale_y_continuous(name="Net Stream Leakage [m3/d]")

p.mo.sfr <- 
  ggplot(subset(df, stream_BC=="SFR"), aes(x=days_since_start/365, y=leakage)) +
  geom_line(color=col.cat.blu) +
  scale_x_continuous(name=NULL, 
                     limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0)) +
  scale_y_continuous(name=NULL)

# annual plot
p.ann.change.riv <-
  ggplot(subset(df.ann, stream_BC=="RIV"), aes(x=(year-1971), y=change_range)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  geom_point() +
  scale_x_continuous(name=NULL, 
                     limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0)) +
  scale_y_continuous(name="Change in Stream Leakage Range\nfrom Previous Year [m3/d]")

p.ann.change.sfr <-
  ggplot(subset(df.ann, stream_BC=="SFR"), aes(x=(year-1971), y=change_range)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  geom_point() +
  scale_x_continuous(name=NULL, 
                     limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0)) +
  scale_y_continuous(name=NULL)

p.ann.change_prc.riv <-
  ggplot(subset(df.ann, stream_BC=="RIV"), aes(x=(year-1971), y=change_range_prc)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  geom_point() +
  scale_x_continuous(name="Year of Transient Spin-Up", 
                     limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0)) +
  scale_y_continuous(name="Change in Stream Leakage Range\nfrom Previous Year [% of range]", labels=scales::percent)

p.ann.change_prc.sfr <-
  ggplot(subset(df.ann, stream_BC=="SFR"), aes(x=(year-1971), y=change_range_prc)) +
  geom_hline(yintercept=0, color=col.gray) +
  geom_line() +
  geom_point() +
  scale_x_continuous(name="Year of Transient Spin-Up", 
                     limits=c(0,20), breaks=seq(0,20,5), expand=c(0,0)) +
  scale_y_continuous(name=NULL, labels=scales::percent) +
  coord_cartesian(ylim=c(-0.01, max(subset(df.ann, stream_BC=="SFR")$change_range_prc, na.rm=T)))

# align plots
p1.riv <- ggplotGrob(p.mo.riv + theme(plot.margin=margin(1,3,1,1, "mm")) + labs(title="RIV"))
p2.riv <- ggplotGrob(p.ann.change.riv + theme(plot.margin=margin(1,3,1,1, "mm")))
p3.riv <- ggplotGrob(p.ann.change_prc.riv + theme(plot.margin=margin(1,3,1,1, "mm")))
p.riv <- rbind(p1.riv, p2.riv, p3.riv, size="first")
p.riv$widths <- unit.pmax(p1.riv$widths, p2.riv$widths, p3.riv$widths)

p1.sfr <- ggplotGrob(p.mo.sfr + theme(plot.margin=margin(1,2,1,2, "mm")) + labs(title="SFR"))
p2.sfr <- ggplotGrob(p.ann.change.sfr + theme(plot.margin=margin(1,2,1,2, "mm")))
p3.sfr <- ggplotGrob(p.ann.change_prc.sfr + theme(plot.margin=margin(1,2,1,2, "mm")))
p.sfr <- rbind(p1.sfr, p2.sfr, p3.sfr, size="first")
p.sfr$widths <- unit.pmax(p1.sfr$widths, p2.sfr$widths, p3.sfr$widths)

ggsave(file.path("figures+tables", "Figure_MODFLOWSpinUp.png"),
       grid.arrange(p.riv, p.sfr, ncol=2, widths=c(1,0.91)),
       width=190, height=210, units="mm")
