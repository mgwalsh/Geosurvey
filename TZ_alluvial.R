# Alluvial diagram of TZ GeoSurvey data
# M. Walsh, September 2017

# install alluvial package
# require(devtools)
# install_github("mbojan/alluvial")

require(downloader)
require(alluvial)

# Data setup --------------------------------------------------------------
# Create a data folder in  your current working directory
dir.create("TZ_data", showWarnings=F)
setwd("./TZ_data")

# download GeoSurvey data
download("https://www.dropbox.com/s/5k4awl2se6s982y/TZ_geos_082317.csv.zip?raw=1", "TZ_geos_082317.csv.zip", mode="wb")
unzip("TZ_geos_082317.csv.zip", overwrite=T)
geos <- read.table("TZ_geos_082317.csv", header=T, sep=",")

geof <- as.data.frame(table(geos$BP, geos$CP, geos$WP))
colnames(geof) <- c("BP","CP","WP","Freq")
geof

# <alluvial> diagram ------------------------------------------------------
alluvial(geof[,1:3], freq=geof$Freq, border=NA, col=ifelse(geof$BP == "Y", "red", "gray"))
