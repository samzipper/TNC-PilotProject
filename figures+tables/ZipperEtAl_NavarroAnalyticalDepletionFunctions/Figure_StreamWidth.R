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
