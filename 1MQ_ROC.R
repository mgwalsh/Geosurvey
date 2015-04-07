#' Receiver-Operator Characteristics (ROC) evaluation of Africa-wide 1MQ GeoSurvey
#' cropland and human settlement predictions with additional GeoSurvey expert test data.
#' M.Walsh, J.Chen and A.Verlinden, April 2015

#+ Required packages
# install.packages(c("downloader","raster","rgdal","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(dismo)

#+ Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("1MQ_data", showWarnings=F)
dat_dir <- "./1MQ_data"

# download 1MQ test data
download("https://www.dropbox.com/s/pt86fr3ko379f8h/1MQ_validation_data.csv?dl=0", "./1MQ_data/1MQ_validation_data.csv", mode="wb")
geosv <- read.table(paste(dat_dir, "/1MQ_validation_data.csv", sep=""), header=T, sep=",")

# BIG download of prediction grids (~447 Mb) and stack in raster
download("https://www.dropbox.com/s/nulz4r3395mh7t5/1MQ_pred_grids.zip?dl=0", "./1MQ_data/1MQ_pred_grids.zip", mode="wb")
unzip("./1MQ_data/1MQ_pred_grids.zip", exdir="./1MQ_data", overwrite=T)
glist <- list.files(path="./1MQ_data", pattern="tif", full.names=T)
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

#+ 1MQ classifier performance evaluation ----------------------------------
# Cropland boosting classifier
gbmcrp <- subset(gsexv, CRPc=="Y", select=c(CRP_gbm))
gbmcra <- subset(gsexv, CRPc=="N", select=c(CRP_gbm))
gbmcrp.eval <- evaluate(p=gbmcrp[,1], a=gbmcra[,1]) ## calculate ROC's on test set <dismo>
gbmcrp.eval
plot(gbmcrp.eval, "ROC")

# Cropland neural network classifier
nncrp <- subset(gsexv, CRPc=="Y", select=c(CRP_nn))
nncra <- subset(gsexv, CRPc=="N", select=c(CRP_nn))
nncrp.eval <- evaluate(p=nncrp[,1], a=nncra[,1]) ## calculate ROC's on test set <dismo>
nncrp.eval
plot(nncrp.eval, "ROC")

# Cropland random forest classifier
rfcrp <- subset(gsexv, CRPc=="Y", select=c(CRP_rf))
rfcra <- subset(gsexv, CRPc=="N", select=c(CRP_rf))
rfcrp.eval <- evaluate(p=rfcrp[,1], a=rfcra[,1]) ## calculate ROC's on test set <dismo>
rfcrp.eval
plot(rfcrp.eval, "ROC")

# Cropland 1MQ ensemble classifier
enscrp <- subset(gsexv, CRPc=="Y", select=c(CRP_ens))
enscra <- subset(gsexv, CRPc=="N", select=c(CRP_ens))
enscrp.eval <- evaluate(p=enscrp[,1], a=enscra[,1]) ## calculate ROC's on test set <dismo>
enscrp.eval
plot(enscrp.eval, "ROC")
enscrp.thld <- threshold(enscrp.eval, "spec_sens") ## TPR+TNR threshold for classification
CRP_ens_mask <- grid$CRP_ens > enscrp.thld
plot(CRP_ens_mask, axes = F, legend = F)

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

# Building/rural settlement 1MQ ensemble classifier
enshsp <- subset(gsexv, HSP=="Y", select=c(RSP_ens))
enshsa <- subset(gsexv, HSP=="N", select=c(RSP_ens))
enshsp.eval <- evaluate(p=enshsp[,1], a=enshsa[,1]) ## calculate ROC's on test set <dismo>
enshsp.eval
plot(enshsp.eval, "ROC")
enshsp.thld <- threshold(enshsp.eval, "spec_sens") ## TPR+TNR threshold for classification
RSP_ens_mask <- grid$RSP_ens > enshsp.thld
plot(RSP_ens_mask, axes = F, legend = F)
