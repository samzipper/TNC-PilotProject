## Navarro_Analytical_Transient_CompareMethods.R

source(file.path("src", "paths+packages.R"))
require(streamDepletr)

# what simulation?
timeType <- "Transient"  # "Transient" or "Intermittent"

## what is the pumping rate?
Qw <- -6*100*0.00378541  # [m3/d]

## threshold for analysis
f.thres <- 0.001

## choose modflow version
modflow_v <- "mfnwt"  # "mfnwt" or "mf2005"
masked <- F           # use masked depletion apportionment? should only be T for SFR
stream_BC <- "RIV"    # stream boundary condition to use for setting steady-state head (screen interval)
timeType  <- "Transient" # "Transient" or "Intermittent"

# which methods to analyze?
methods.plot <- c("Qf.InvDistSq", "Qf.WebSq", "Qf.TPoly")

#### (0) Prep spatial data

## load well locations
df.wel <- read.table(file.path("modflow", "input", "iwel.txt"), sep=" ", header=T)

## load output from steady-state, no pumping scenario - this is used to define the screen interval
## so that screen interval is consistent between MODFLOW and analytical
m.wte <- as.matrix(read.csv(file.path("modflow", "Navarro-SteadyState", stream_BC, modflow_v, "wte.csv")), header=F)

# grab steady-state head based on row/col (need to add 1 because python is 0-based indexing)
df.wel$wte_m <- m.wte[as.matrix(df.wel[,c("row", "col")])+1]

# make a spatial points data frame
xy <- df.wel[,c("lon", "lat")]
spdf.wel <- SpatialPointsDataFrame(coords = xy, data = df.wel,
                                   proj4string = CRS(crs.MODFLOW))

## load stream data - created in MODFLOW_Navarro_InputPrepData.R
shp.streams <- readOGR(dsn=file.path("modflow", "input"), layer="iriv")

# shapefile shortens names; rename them
names(shp.streams) <- c("OBJECTID", "REACHCODE", "TerminalPa", "lineLength_m", "TotDASqKM", "StreamOrde", 
                        "TerminalFl", "SLOPE", "FromNode", "ToNode", "SegNum")

# figure out which segments are part of Navarro
segs.navarro <- shp.streams@data$SegNum[shp.streams@data$TerminalPa==outlet.TerminalPa]

# domain boundary shapefile
shp <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU10_Navarro")
shp.UTM <- spTransform(shp, crs.MODFLOW)

shp.adj <- readOGR(dsn=file.path("data", "NHD", "WBD"), layer="WBDHU12_Navarro+Adjacent")
shp.adj.UTM <- spTransform(shp.adj, crs.MODFLOW)

## prep polygon boundaries for plots
df.basin <- tidy(shp.UTM)
df.basin.adj <- tidy(shp.adj.UTM)
df.riv <- tidy(shp.streams)

#### (1) Load MODFLOW results and figure out timesteps and wells for comparison 
####     (MODFLOW_HTC_Navarro_CalculateCaptureFraction.R)
df_MODFLOW_RIV <- 
  file.path("modflow","HTC", "Navarro", timeType, "RIV", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
  transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
            Time = round(Time, 1)) %>% 
  subset(depletion.prc.modflow > f.thres)

df_MODFLOW_SFR <- 
  file.path("modflow","HTC", "Navarro", timeType, "SFR", modflow_v, "Depletion_MODFLOW.csv") %>% 
  read.csv(stringsAsFactors=F) %>% 
  #  transform(depletion.prc.modflow = depletion_m3.d/Qw_m3.d,  # use net pumping rate from MODFLOW
  transform(depletion.prc.modflow = depletion_m3.d/Qw,          # use prescribed pumping rate
            Time = round(Time, 1)) %>% 
  subset(depletion.prc.modflow > f.thres)

# combine into one
df_MODFLOW <- 
  full_join(df_MODFLOW_RIV[,c("SegNum", "WellNum", "Time", "depletion.prc.modflow")], 
            df_MODFLOW_SFR[,c("SegNum", "WellNum", "Time", "depletion.prc.modflow")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  set_colnames(c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")) %>% 
  subset(SegNum %in% segs.navarro)
df_MODFLOW[is.na(df_MODFLOW)] <- 0

#### (2) Load analytical results and combine with MODFLOW
####     (Navarro_Analytical_[timeType].R)
# NoApportionment
df_no <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_NoApportionment_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  subset(Qf.NoApport > f.thres) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", "analytical", Qf.NoApport)) %>% 
  transform(Time = round(Time, 1))

# AdjacentOnly
df_adj_g <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_AdjacentOnly_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="glover"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="glover" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "glover", 
            apportionment = "AdjacentOnly") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

df_adj_h <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_AdjacentOnly_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="hunt"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="hunt" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "hunt", 
            apportionment = "AdjacentOnly") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

# LocalArea
df_loc_g <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_LocalArea_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="glover"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="glover" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "glover", 
            apportionment = "LocalArea") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

df_loc_h <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_LocalArea_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="hunt"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="hunt" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "hunt", 
            apportionment = "LocalArea") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

# Dynamic
df_dyn_g <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_Dynamic_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="glover"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="glover" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "glover", 
            apportionment = "Dynamic") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

df_dyn_h <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_Dynamic_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="hunt"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="hunt" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "hunt", 
            apportionment = "Dynamic") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

# WholeDomain
df_whl_g <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_WholeDomain_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="glover"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="glover" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "glover", 
            apportionment = "WholeDomain") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

df_whl_h <- 
  file.path("results", paste0("Depletion_Analytical_", timeType, "_WholeDomain_AllMethods+Wells+Reaches.csv")) %>% 
  read.csv(stringsAsFactors = F) %>% 
  subset(SegNum %in% segs.navarro) %>% 
  transform(Time = round(Time, 1)) %>% 
  full_join(subset(df_no, analytical=="hunt"), by=c("WellNum", "SegNum", "Time", "analytical")) %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0)) %>% 
  subset(analytical=="hunt" & (Qf.InvDistSq > f.thres | Qf.WebSq > f.thres | Qf.TPoly > f.thres | Qf.NoApport > f.thres)) %>% 
  dplyr::select(c("WellNum", "SegNum", "Time", methods.plot, "Qf.NoApport")) %>% 
  full_join(df_MODFLOW[,c("SegNum", "WellNum", "Time", "Qf.RIV", "Qf.SFR")],
            by=c("SegNum", "WellNum", "Time")) %>% 
  transform(analytical = "hunt", 
            apportionment = "WholeDomain") %>% 
  replace_na(list(Qf.NoApport = 0, Qf.TPoly = 0, Qf.InvDistSq = 0, Qf.WebSq = 0, Qf.RIV = 0, Qf.SFR = 0))

## combine into one data frame and melt into analytical and modflow columns
df <- 
  rbind(df_adj_g, df_adj_h, 
        df_dyn_g, df_dyn_h, 
        df_loc_g, df_loc_h, 
        df_whl_g, df_whl_h) %>% 
  melt(id=c("analytical", "WellNum", "SegNum", "Time", "Qf.RIV", "Qf.SFR", "apportionment"),
       value.name="depletion.prc", variable.name="method") %>% 
  melt(id=c("analytical", "WellNum", "SegNum", "Time", "apportionment", "method", "depletion.prc"),
       value.name="depletion.prc.modflow", variable.name="stream_BC") %>% 
  subset((depletion.prc > f.thres | depletion.prc.modflow > f.thres))

# check to make sure TPoly works: these should be the same
sum(df$method=="Qf.TPoly" & df$apportionment=="AdjacentOnly" & df$stream_BC=="Qf.RIV" & df$Time==max(df$Time) & df$analytical=="glover")
sum(df$method=="Qf.TPoly" & df$apportionment=="WholeDomain" & df$stream_BC=="Qf.RIV" & df$Time==max(df$Time) & df$analytical=="glover")

# calculate fit statistics
df.fit <- 
  df %>% 
  transform(stream_BC = gsub("Qf.", "", stream_BC)) %>% 
  group_by(method, apportionment, analytical, Time, stream_BC) %>% 
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
  transform(stream_BC = gsub("Qf.", "", stream_BC)) %>% 
  group_by(method, apportionment, analytical, Time, stream_BC) %>% 
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

df.fit %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=method)) +
  geom_line() +
  facet_grid(stream_BC+analytical~apportionment, 
             labeller=as_labeller(c(labels.apportionment, 
                                    c("RIV"="RIV", "SFR"="SFR"), 
                                    c("glover"="glover", "hunt"="hunt")))) +
  scale_color_manual(values=pal.method.Qf)

df.fit.big %>% 
  ggplot(aes(x=Time, y=KGE.overall, color=method)) +
  geom_line() +
  facet_grid(stream_BC+analytical~apportionment, 
             labeller=as_labeller(c(labels.apportionment, 
                                    c("RIV"="RIV", "SFR"="SFR"), 
                                    c("glover"="glover", "hunt"="hunt")))) +
  scale_color_manual(values=pal.method.Qf)
