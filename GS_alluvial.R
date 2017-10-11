# Alluvial diagrams of TZ, ET, NG & GH GeoSurvey data
# M. Walsh, September 2017

# install alluvial package
# require(devtools)
# install_github("mbojan/alluvial")

require(downloader)
require(alluvial)

# Data setup --------------------------------------------------------------
# Create a data folder in  your current working directory
dir.create("GS_data", showWarnings=F)
setwd("./GS_data")

# download GeoSurvey data
# Tanzania
download("https://www.dropbox.com/s/57kuxbkm5sv092a/TZ_geos_2017.csv.zip?raw=1", "TZ_geos_2017.csv.zip", mode="wb")
unzip("TZ_geos_2017.csv.zip", overwrite=T)
TZ_geos <- read.table("TZ_geos_2017.csv", header=T, sep=",")
TZ_geos <- TZ_geos[complete.cases(TZ_geos), ]

TZ_geof <- as.data.frame(table(TZ_geos$BP, TZ_geos$CP, TZ_geos$WP))
colnames(TZ_geof) <- c("BP","CP","WP","Freq")
TZ_geof

# Ethiopia
download("https://www.dropbox.com/s/mfw02vsrrit674z/ET_geos_2017.csv.zip?raw=1", "ET_geos_2017.csv.zip", mode="wb")
unzip("ET_geos_2017.csv.zip", overwrite=T)
ET_geos <- read.table("ET_geos_2017.csv", header=T, sep=",")
ET_geos <- ET_geos[complete.cases(ET_geos), ]

ET_geof <- as.data.frame(table(ET_geos$BP, ET_geos$CP, ET_geos$WP))
colnames(ET_geof) <- c("BP","CP","WP","Freq")
ET_geof

# Nigeria
download("https://www.dropbox.com/s/iytk6ljh16wwpds/NG_geos_2017.csv.zip?raw=1", "NG_geos_2017.csv.zip", mode="wb")
unzip("NG_geos_2017.csv.zip", overwrite=T)
NG_geos <- read.table("NG_geos_2017.csv", header=T, sep=",")
NG_geos <- NG_geos[complete.cases(NG_geos), ]

NG_geof <- as.data.frame(table(NG_geos$BP, NG_geos$CP, NG_geos$WP))
colnames(NG_geof) <- c("BP","CP","WP","Freq")
NG_geof

# Ghana
download("https://www.dropbox.com/s/s5b1katorq0s1mq/GH_geos_2017.csv.zip?raw=1", "GH_geos_2017.csv.zip", mode="wb")
unzip("GH_geos_2017.csv.zip", overwrite=T)
GH_geos <- read.table("GH_geos_2017.csv", header=T, sep=",")
GH_geos <- GH_geos[complete.cases(GH_geos), ]

GH_geof <- as.data.frame(table(GH_geos$BP, GH_geos$CP, GH_geos$WP))
colnames(GH_geof) <- c("BP","CP","WP","Freq")
GH_geof

# <alluvial> diagrams -----------------------------------------------------
par(mfrow=c(2,2), mar=c(5,5,1,1))
alluvial(TZ_geof[,1:3], freq=TZ_geof$Freq, border=NA, col=ifelse(TZ_geof$BP == "Y", "red", "gray"))
alluvial(ET_geof[,1:3], freq=ET_geof$Freq, border=NA, col=ifelse(ET_geof$BP == "Y", "red", "gray"))
alluvial(NG_geof[,1:3], freq=NG_geof$Freq, border=NA, col=ifelse(NG_geof$BP == "Y", "red", "gray"))
alluvial(GH_geof[,1:3], freq=GH_geof$Freq, border=NA, col=ifelse(GH_geof$BP == "Y", "red", "gray"))
dev.off()
