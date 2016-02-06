# Country-level GeoSurvey validation results
# M. Walsh, February 2016

# install alluvial package
# require(devtools)
# install_github("mbojan/alluvial")

require(downloader)
require(alluvial)

# Data setup --------------------------------------------------------------
# Create a data folder in  your current working directory
dir.create("GS_val", showWarnings=F)
setwd("./GS_val")

# Download
download("https://www.dropbox.com/s/vfvqdt78bmvckic/CID.csv?dl=0", "CID.csv", mode="wb")
gscid <- read.table("CID.csv", header=T, sep=",")

gscrp <- as.data.frame(table(gscid$CID, gscid$CP, gscid$BP, gscid$WP))
colnames(gscrp) <- c("Country","CP","BP","WP","Freq")

# <alluvial> diagram ------------------------------------------------------
alluvial(gscrp[,c(2:4,1)], freq=gscrp$Freq, border=NA,
         hide = gscrp$Freq < quantile(gscrp$Freq, 0.25),
         col=ifelse(gscrp$CP == "Yes", "red", "gray"))

