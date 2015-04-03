#' Ensemble validation of Tanzania 1M GeoSurvey cropland and human settlement observations.
#' M.Walsh & J.Chen April 2015

#+ Required packages
# install.packages(c("downloader","raster","rgdal","caret","pROC","glmnet","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(pROC)
require(glmnet)
require(dismo)

#+ Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("TZ_1MGS_data", showWarnings=F)
dat_dir <- "./TZ_1MGS_data"

# download Tanzania test data
download("https://www.dropbox.com/s/gfgjnrgllwqt79d/TZ_geos_123114.csv?dl=0", "./TZ_1MGS_data/TZ_geos_123114.csv", mode="wb")
geosv <- read.table(paste(dat_dir, "/TZ_geos_123114.csv", sep=""), header=T, sep=",")

# download Tanzania prediction grids (~18.6 Mb) and stack in raster
download("https://www.dropbox.com/s/semwhdlxu9a0863/TZ_preds.zip?dl=0", "./TZ_1MGS_data/TZ_preds.zip", mode="wb")
unzip("./TZ_1MGS_data/TZ_preds.zip", exdir="./TZ_1MGS_data", overwrite=T)
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

# Regularized ensemble weighting <glmnet> -------------------------------
# 10-fold CV
ens <- trainControl(method = "cv", number = 10, classProbs = T, summaryFunction = twoClassSummary)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.ens <- train(CRP ~ CRP_gbm + CRP_nn + CRP_rf, data = gsexv,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "ROC",
                 trControl = ens)
CRP.ens
crp.pred <- predict(CRP.ens, gsexv, type="prob")
crp.test <- cbind(gsexv, crp.pred)
crp <- subset(crp.test, CRP=="Y", select=c(Y))
cra <- subset(crp.test, CRP=="N", select=c(Y))
crp.eval <- evaluate(p=crp[,1], a=cra[,1]) ## calculate ROC's on test set <dismo>
crp.eval
plot(crp.eval, 'ROC') ## plot ROC curve
crp.thld <- threshold(crp.eval, 'spec_sens') ## TPR+TNR threshold for classification
crpens.pred <- predict(grid, CRP.ens, type="prob") ## spatial prediction
plot(1-crpens.pred, axes = F)
crpmask <- 1-crpens.pred > crp.thld
plot(crpmask, axes = F, legend = F)

# presence/absence of Buildings/Rural Settlements (HSP, present = Y, absent = N)
HSP.ens <- train(HSP ~ RSP_gbm + RSP_nn + RSP_rf, data = gsexv,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "ROC",
                 trControl = ens)
HSP.ens
hsp.pred <- predict(HSP.ens, gsexv, type="prob")
hsp.test <- cbind(gsexv, hsp.pred)
hsp <- subset(hsp.test, HSP=="Y", select=c(Y))
hsa <- subset(hsp.test, HSP=="N", select=c(Y))
hsp.eval <- evaluate(p=hsp[,1], a=hsa[,1]) ## calculate ROC's on test set <dismo>
hsp.eval
plot(hsp.eval, 'ROC') ## plot ROC curve
hsp.thld <- threshold(hsp.eval, 'spec_sens') ## TPR+TNR threshold for classification
hspens.pred <- predict(grid, HSP.ens, type="prob") ## spatial prediction
plot(1-hspens.pred, axes = F)
hspmask <- 1-hspens.pred > hsp.thld
plot(hspmask, axes = F, legend = F)

#+ Write spatial predictions -----------------------------------------------
# Create a "Results" folder in current working directory
dir.create("TZ_1MGS_results", showWarnings=F)

# Export Gtif's to "./TZ_results"
enspred <- stack(1-crpens.pred, crpmask, 1-hspens.pred, hspmask)
writeRaster(enspred, filename="./TZ_1MGS_results/TZ_1MGS_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
