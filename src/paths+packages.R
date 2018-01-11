## paths+packages.R
# This script holds packages, paths, color schemes, etc which are shared among scripts.

## packages
require(ggplot2)
require(dplyr)
require(reshape2)
require(lubridate)
require(grid)
require(gridExtra)
require(hydrostats)
require(waterData)

## labels
labs.mo <- c("1"="Jan", "2"="Feb", "3"="Mar", "4"="Apr", "5"="May", "6"="Jun",
             "7"="Jul", "8"="Aug", "9"="Sep", "10"="Oct", "11"="Nov", "12"="Dec")

## ggplot theme
windowsFonts(Arial=windowsFont("TT Arial"))
theme_scz <- function(...){
  theme_bw(base_size=10, base_family="Arial") + 
    theme(
      text=element_text(color="black"),
      axis.title=element_text(face="bold", size=rel(1)),
      axis.text=element_text(size=rel(1)),
      strip.text=element_text(size=rel(1)),
      legend.title=element_text(face="bold", size=rel(1)),
      legend.text=element_text(size=rel(1)),
      panel.grid=element_blank())
}

## functions
# extract p-value from linear model fits
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
