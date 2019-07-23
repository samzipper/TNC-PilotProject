## Figure_StreamWidth.R

source("src/paths+packages.R")

## get data from USGS - output from script Navarro_StreamflowData.R
df <- read.csv(file.path("data", "Navarro_StreamWidthEstimates.csv"), skip=1)

## best fit function - from Excel
line_graph <- function(x) 9.7133*exp(0.0023*x)

## plots
ggplot(data.frame(x = c(0,max(df$Drainage.Area..km2.)))) +
  geom_point(data=df, aes(y=Width..m., x=Drainage.Area..km2.)) +
  scale_x_continuous(name = "Drainage Area [km^2]") +
  scale_y_continuous(name = "Stream Segment Width [m]") +
  stat_function(fun = line_graph, geom = "path", color = col.cat.blu) +
  ggsave(file.path("figures+tables", "ZipperEtAl_NavarroAnalyticalDepletionFunctions", "Figure_StreamWidth.png"),
         width=95, height=65, units="mm")

# add reeves et al estimates
df$width_m_CA <- line_graph(df$Drainage.Area..km2.)


reeves_width <-
  function(da_km2) {
    da_mi2 <- da_km2/(1.609344*1.609344)
    (3.28*(10^(0.522358*log10(da_mi2*1.6093*1.6093)) - 0.18786))*0.3048
  }


df$width_m_MI <- reeves_width(df$Drainage.Area..km2.)

ggplot(data.frame(x = c(0,max(df$Drainage.Area..km2.)))) +
  geom_point(data=df, aes(y=Width..m., x=Drainage.Area..km2.)) +
  scale_x_continuous(name = "Drainage Area [km^2]") +
  scale_y_continuous(name = "Stream Segment Width [m]") +
  stat_function(fun = line_graph, geom = "path", color = col.cat.blu) +
  stat_function(fun = reeves_width, geom = "path", color = col.cat.red)

MSE(df$width_m_CA, df$Width..m.)
MSE(df$width_m_MI, df$Width..m.)

R2(df$width_m_CA, df$Width..m.)
R2(df$width_m_MI, df$Width..m.)
