#' Evaluation and local stacking of Tanzania 1M GeoSurvey cropland and human settlement
#' predictions with Tanzania GeoSurvey test data.
#' M.Walsh, J.Chen & A.Verlinden, April 2015

#+ Required packages
# install.packages(c("downloader","raster","rgdal","dismo","caret","glmnet")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(dismo)
require(caret)
require(glmnet)

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
# Cropland boosting classifier
gbmcrp <- subset(gsexv, CRP=="Y", select=c(CRP_gbm))
gbmcra <- subset(gsexv, CRP=="N", select=c(CRP_gbm))
gbmcrp.eval <- evaluate(p=gbmcrp[,1], a=gbmcra[,1]) ## calculate ROC's on test set <dismo>
gbmcrp.eval
plot(gbmcrp.eval, "ROC")

# Cropland neural network classifier
nncrp <- subset(gsexv, CRP=="Y", select=c(CRP_nn))
nncra <- subset(gsexv, CRP=="N", select=c(CRP_nn))
nncrp.eval <- evaluate(p=nncrp[,1], a=nncra[,1]) ## calculate ROC's on test set <dismo>
nncrp.eval
plot(nncrp.eval, "ROC")

# Cropland random forest classifier
rfcrp <- subset(gsexv, CRP=="Y", select=c(CRP_rf))
rfcra <- subset(gsexv, CRP=="N", select=c(CRP_rf))
rfcrp.eval <- evaluate(p=rfcrp[,1], a=rfcra[,1]) ## calculate ROC's on test set <dismo>
rfcrp.eval
plot(rfcrp.eval, "ROC")

# Cropland 1MGS ensemble classifier
enscrp <- subset(gsexv, CRP=="Y", select=c(CRP_ens))
enscra <- subset(gsexv, CRP=="N", select=c(CRP_ens))
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

# Building/rural settlement 1MGS ensemble classifier
enshsp <- subset(gsexv, HSP=="Y", select=c(RSP_ens))
enshsa <- subset(gsexv, HSP=="N", select=c(RSP_ens))
enshsp.eval <- evaluate(p=enshsp[,1], a=enshsa[,1]) ## calculate ROC's on test set <dismo>
enshsp.eval
plot(enshsp.eval, "ROC")
enshsp.thld <- threshold(enshsp.eval, "spec_sens") ## TPR+TNR threshold for classification
RSP_ens_mask <- grid$RSP_ens > enshsp.thld
plot(RSP_ens_mask, axes = F, legend = F)

#+ Local classifier (re)stacking ------------------------------------------
# 10-fold CV
lcs <- trainControl(method = "cv", number = 10, classProbs = T, summaryFunction = twoClassSummary)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.lcs <- train(CRP ~ CRP_gbm + CRP_nn + CRP_rf, data = gsexv,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "ROC",
                 trControl = lcs)
CRP.lcs
crp.pred <- predict(CRP.lcs, gsexv, type="prob")
crp.test <- cbind(gsexv, crp.pred)
lcscrp <- subset(crp.test, CRP=="Y", select=c(Y))
lcscra <- subset(crp.test, CRP=="N", select=c(Y))
lcscrp.eval <- evaluate(p=lcscrp[,1], a=lcscra[,1]) ## calculate ROC's on test set <dismo>
lcscrp.eval
plot(lcscrp.eval, "ROC") ## plot ROC curve
lcscrp.thld <- threshold(crp.eval, 'spec_sens') ## TPR+TNR threshold for classification
CRP_lcs <- predict(grid, CRP.lcs, type="prob") ## spatial prediction
plot(1-CRP_lcs, axes = F)
CRP_lcs_mask <- 1-CRP_lcs > lcscrp.thld
plot(CRP_lcs_mask, axes = F, legend = F)

# presence/absence of Buildings/rural settlements (HSP, present = Y, absent = N)
RSP.lcs <- train(HSP ~ RSP_gbm + RSP_nn + RSP_rf, data = gsexv,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "ROC",
                 trControl = lcs)
RSP.lcs
rsp.pred <- predict(RSP.lcs, gsexv, type="prob")
rsp.test <- cbind(gsexv, rsp.pred)
lcsrsp <- subset(rsp.test, HSP=="Y", select=c(Y))
lcsrsa <- subset(rsp.test, HSP=="N", select=c(Y))
lcsrsp.eval <- evaluate(p=lcsrsp[,1], a=lcsrsa[,1]) ## calculate ROC's on test set <dismo>
lcsrsp.eval
plot(lcsrsp.eval, "ROC") ## plot ROC curve
lcsrsp.thld <- threshold(lcsrsp.eval, 'spec_sens') ## TPR+TNR threshold for classification
RSP_lcs <- predict(grid, RSP.lcs, type="prob") ## spatial prediction
plot(1-RSP_lcs, axes = F)
RSP_lcs_mask <- 1-RSP_lcs > lcsrsp.thld
plot(RSP_lcs_mask, axes = F, legend = F)

#+ Write spatial predictions -----------------------------------------------
# Create a "Results" folder in current working directory
dir.create("TZ_1MGS_results", showWarnings=F)
LCS_pred <- stack(1-CRP_lcs, CRP_lcs_mask, 1-RSP_lcs, RSP_lcs_mask)
writeRaster(LCS_pred, filename="./TZ_1MGS_results/TZ_lcs_pred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)