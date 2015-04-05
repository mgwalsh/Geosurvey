#' Ensemble validation of Tanzania 1M GeoSurvey cropland and human settlement observations.
#' M.Walsh & J.Chen April 2015

#+ Required packages
# install.packages(c("downloader","raster","rgdal","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(dismo)

#+ Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("TZ_1MGS_data", showWarnings=F)
dat_dir <- "./TZ_1MGS_data"

# download Tanzania test data
download("https://www.dropbox.com/s/gfgjnrgllwqt79d/TZ_geos_123114.csv?dl=0", "./TZ_1MGS_data/TZ_geos_123114.csv", mode="wb")
geosv <- read.table(paste(dat_dir, "/TZ_geos_123114.csv", sep=""), header=T, sep=",")

# download Tanzania prediction grids (~21.1 Mb) and stack in raster
download("https://www.dropbox.com/s/w8l41t5muc1rr4j/TZ_1MQ_preds.zip?dl=0", "./TZ_1MGS_data/TZ_1MQ_preds.zip", mode="wb")
unzip("./TZ_1MGS_data/TZ_1MQ_preds.zip", exdir="./TZ_1MGS_data", overwrite=T)
glist <- list.files(path="./TZ_1MGS_data", pattern="tif", full.names=T)
grid <- stack(glist)

#+ Data setup --------------------------------------------------------------
# Project test data to grid CRS
geosv.proj <- as.data.frame(project(cbind(geosv$Lon, geosv$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geosv.proj) <- c("x","y")
geosv <- cbind(geosv, geosv.proj)
coordinates(geosv) <- ~x+y
projection(geosv) <- projection(grid)

# Extract gridded variables to test data observations
gsexv <- data.frame(coordinates(geosv), geosv$CRP, geosv$HSP, extract(grid, geosv))
gsexv <- na.omit(gsexv)
colnames(gsexv)[3:4] <- c("CRP", "HSP")

# 1 MGS classifier performance evaluation ---------------------------------
# Cropland boosting
gbmcrp <- subset(gsexv, CRP=="Y", select=c(CRP_gbm))
gbmcra <- subset(gsexv, CRP=="N", select=c(CRP_gbm))
gbmcrp.eval <- evaluate(p=gbmcrp[,1], a=gbmcra[,1]) ## calculate ROC's on test set <dismo>
gbmcrp.eval
plot(gbmcrp.eval, "ROC")

# Cropland neural network
nncrp <- subset(gsexv, CRP=="Y", select=c(CRP_nn))
nncra <- subset(gsexv, CRP=="N", select=c(CRP_nn))
nncrp.eval <- evaluate(p=nncrp[,1], a=nncra[,1]) ## calculate ROC's on test set <dismo>
nncrp.eval
plot(nncrp.eval, "ROC")

# Cropland random forest
rfcrp <- subset(gsexv, CRP=="Y", select=c(CRP_rf))
rfcra <- subset(gsexv, CRP=="N", select=c(CRP_rf))
rfcrp.eval <- evaluate(p=rfcrp[,1], a=rfcra[,1]) ## calculate ROC's on test set <dismo>
rfcrp.eval
plot(rfcrp.eval, "ROC")

# Cropland ensemble
enscrp <- subset(gsexv, CRP=="Y", select=c(CRP_ens))
enscra <- subset(gsexv, CRP=="N", select=c(CRP_ens))
enscrp.eval <- evaluate(p=enscrp[,1], a=enscra[,1]) ## calculate ROC's on test set <dismo>
enscrp.eval
plot(enscrp.eval, "ROC")

# Building/rural settlement boosting classifier
gbmhsp <- subset(gsexv, HSP=="Y", select=c(RSP_gbm))
gbmhsa <- subset(gsexv, HSP=="N", select=c(RSP_gbm))
gbmhsp.eval <- evaluate(p=gbmhsp[,1], a=gbmhsa[,1]) ## calculate ROC's on test set <dismo>
gbmhsp.eval
plot(gbmhsp.eval, "ROC")

# Building/rural settlement neural network classifier
nnhsp <- subset(gsexv, HSP=="Y", select=c(RSP_nn))
nnhsa <- subset(gsexv, HSP=="N", select=c(RSP_nn))
nnhsp.eval <- evaluate(p=nnhsp[,1], a=nnhsa[,1]) ## calculate ROC's on test set <dismo>
nnhsp.eval
plot(nnhsp.eval, "ROC")

# Building/rural settlement random forest classifier
rfhsp <- subset(gsexv, HSP=="Y", select=c(RSP_rf))
rfhsa <- subset(gsexv, HSP=="N", select=c(RSP_rf))
rfhsp.eval <- evaluate(p=rfhsp[,1], a=rfhsa[,1]) ## calculate ROC's on test set <dismo>
rfhsp.eval
plot(rfhsp.eval, "ROC")

# Building/rural settlement ensemble classifier
enshsp <- subset(gsexv, HSP=="Y", select=c(RSP_ens))
enshsa <- subset(gsexv, HSP=="N", select=c(RSP_ens))
enshsp.eval <- evaluate(p=enshsp[,1], a=enshsa[,1]) ## calculate ROC's on test set <dismo>
enshsp.eval
plot(enshsp.eval, "ROC")

#+ Local classifier stacking ----------------------------------------------





