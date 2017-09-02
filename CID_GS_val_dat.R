# Country-level GeoSurvey validation results
# M. Walsh, February 2016

# install alluvial package
# require(devtools)
# install_github("mbojan/alluvial")

require(downloader)
require(alluvial)
require(circlize)

# Data setup --------------------------------------------------------------
# Create a data folder in  your current working directory
dir.create("GS_val", showWarnings=F)
setwd("./GS_val")

# Download
download("https://www.dropbox.com/s/vfvqdt78bmvckic/CID.csv?dl=0", "CID.csv", mode="wb")
gscid <- read.table("CID.csv", header=T, sep=",")

gscrp <- as.data.frame(table(gscid$CID, gscid$CP, gscid$BP, gscid$WP))
colnames(gscrp) <- c("Country","CP","BP","WP","Freq")

# Alluvial diagram --------------------------------------------------------
alluvial(gscrp[,c(1,2:4)], freq=gscrp$Freq, border=NA,
         hide = gscrp$Freq < quantile(gscrp$Freq, 0.25),
         col=ifelse(gscrp$CP == "Yes", "red", "gray"))

# Chord diagrams ----------------------------------------------------------
grid.col <- c(Yes="red", No="gray", NN="gray", NY="gray",YN="gray", YY="gray")

et <- gscid[which(gscid$CID=='ET'), ]
etcrp <- as.data.frame(table(et$CP, et$System))
chordDiagram(etcrp, grid.col = grid.col,
             annotationTrack = "grid", annotationTrackHeight = 0.15)

tz <- gscid[which(gscid$CID=='TZ'), ]
tzcrp <- as.data.frame(table(tz$CP, tz$System))
chordDiagram(tzcrp, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

gh <- gscid[which(gscid$CID=='GH'), ]
ghcrp <- as.data.frame(table(gh$CP, gh$System))
chordDiagram(ghcrp, grid.col = grid.col, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))


