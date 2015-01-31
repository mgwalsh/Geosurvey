#' Validation of ensemble predictions of Tanzania 1M GeoSurvey cropland and
#' human settlement observations.
#' M.Walsh, J.Chen & A.Verlinden, January 2015

#+ Required packages
# install.packages(c("downloader","raster","rgdal","caret","pROC","MASS","randomForest","gbm","nnet","glmnet","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(pROC)
require(MASS)
require(randomForest)
require(gbm)
require(nnet)
require(glmnet)
require(dismo)

#+ Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("TZ_1MGS_data", showWarnings=F)
dat_dir <- "./TZ_1MGS_data"

# download 1M GeoSurvey data
download("https://www.dropbox.com/s/jtwmf7ck9ebh7bm/1MGS_cleaned.csv.zip?dl=0", "./TZ_1MGS_data/1MGS_cleaned.csv.zip", mode="wb")
unzip("./TZ_1MGS_data/1MGS_cleaned.csv.zip", exdir="./TZ_1MGS_data", overwrite=T)
geos <- read.table(paste(dat_dir, "/1MGS_cleaned.csv", sep=""), header=T, sep=",")

# download Tanzania validation data
download("https://www.dropbox.com/s/gfgjnrgllwqt79d/TZ_geos_123114.csv?dl=0", "./TZ_1MGS_data/TZ_geos_123114.csv", mode="wb")
geosv <- read.table(paste(dat_dir, "/TZ_geos_123114.csv", sep=""), header=T, sep=",")

# download Tanzania Gtifs (~27.9 Mb) and stack in raster
download("https://www.dropbox.com/s/otiqe78s0kf1z1s/TZ_grids.zip?dl=0", "./TZ_1MGS_data/TZ_grids.zip", mode="wb")
unzip("./TZ_1MGS_data/TZ_grids.zip", exdir="./TZ_1MGS_data", overwrite=T)
glist <- list.files(path="./TZ_1MGS_data", pattern="tif", full.names=T)
grid <- stack(glist)

#+ Data setup --------------------------------------------------------------
# Project 1M GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grid)

# Project validation data to grid CRS
geosv.proj <- as.data.frame(project(cbind(geosv$Lon, geosv$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geosv.proj) <- c("x","y")
geosv <- cbind(geosv, geosv.proj)
coordinates(geosv) <- ~x+y
projection(geosv) <- projection(grid)

# Extract gridded variables to GeoSurvey locations
gsext <- data.frame(coordinates(geos), geos$CRP, geos$HSP, extract(grid, geos))
gsext <- na.omit(gsext)
colnames(gsext)[3:4] <- c("CRP", "HSP")
# write.csv(gsext, "TZ_1MGS.csv", row.names=F)

# Extract gridded variables to validation observations
gsexv <- data.frame(coordinates(geosv), geosv$CRP, geosv$HSP, extract(grid, geosv))
gsexv <- na.omit(gsext)
colnames(gsext)[3:4] <- c("CRP", "HSP")

# Assemble dataframes
# presence/absence of Cropland (CRP, present = Y, absent = N)
crpdat <- data.frame(gsext$CRP, gsext[,5:28])
colnames(crpdat)[1] <- "CRP"

# presence/absence of Buildings/Settlements (HSP, present = Y, absent = N)
hspdat <- data.frame(gsext$HSP, gsext[,5:28])
colnames(hspdat)[1] <- "HSP"

#+ Stepwise main effects GLM's <MASS> --------------------------------------
# 10-fold CV
step <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.glm <- train(CRP ~ ., data = crpdat,
                 family = binomial, 
                 method = "glmStepAIC",
                 metric = "Kappa",
                 trControl = step)
CRP.glm
crpglm.pred <- predict(grid, CRP.glm, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.glm <- train(HSP ~ ., data = hspdat,
                 family=binomial, 
                 method = "glmStepAIC",
                 metric = "Kappa",
                 trControl = step)
HSP.glm
hspglm.pred <- predict(grid, HSP.glm, type = "prob") ## spatial predictions

#+ Random forests <randomForest> -------------------------------------------
# out-of-bag predictions
oob <- trainControl(method = "oob")

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.rf <- train(CRP ~ ., data = crpdat,
                method = "rf",
                metric = "Kappa",
                trControl = oob)
CRP.rf
crprf.pred <- predict(grid, CRP.rf, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.rf <- train(HSP ~ ., data = hspdat,
                method = "rf",
                metric = "Kappa",
                trControl = oob)
HSP.rf
hsprf.pred <- predict(grid, HSP.rf, type = "prob") ## spatial predictions

#+ Gradient boosting <gbm> ------------------------------------------
# CV for training gbm's
gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.gbm <- train(CRP ~ ., data = crpdat,
                 method = "gbm",
                 metric = "Kappa",
                 trControl = gbm)
CRP.gbm
crpgbm.pred <- predict(grid, CRP.gbm, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.gbm <- train(HSP ~ ., data = hspdat,
                 method = "gbm",
                 metric = "Kappa",
                 trControl = gbm)
HSP.gbm
hspgbm.pred <- predict(grid, HSP.gbm, type = "prob") ## spatial predictions

#+ Neural nets <nnet> ------------------------------------------------------
# CV for training nnet's
nn <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.nn <- train(CRP ~ ., data = crpdat,
                method = "nnet",
                metric = "Kappa",
                trControl = nn)
CRP.nn
crpnn.pred <- predict(grid, CRP.nn, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.nn <- train(HSP ~ ., data = hspdat,
                method = "nnet",
                metric = "Kappa",
                trControl = nn)
HSP.nn
hspnn.pred <- predict(grid, HSP.nn, type = "prob") ## spatial predictions

#+ Plot predictions by GeoSurvey variables ---------------------------------
# Cropland prediction plots
crp.preds <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred)
names(crp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(crp.preds, axes = F)

# Settlement prediction plots
hsp.preds <- stack(1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(hsp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(hsp.preds, axes = F)

#+ Ensemble predictions <glm>, <rf>, <gbm>, <nnet> -------------------------
# Ensemble set-up
pred <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred,
              1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(pred) <- c("CRPglm","CRPrf","CRPgbm","CRPnn",
                 "HSPglm","HSPrf","HSPgbm","HSPnn")

# presence/absence of Cropland (CRP, present = Y, absent = N)
crpens <- data.frame(geosv$CRP, extract(pred, geosv))
crpens <- na.omit(crpens)
colnames(crpens)[1] <- "CRP"

# presence/absence of Buildings/Settlements (HSP, present = Y, absent = N)
hspens <- data.frame(geosv$HSP, extract(pred, geosv))
hspens <- na.omit(hspens)
colnames(hspens)[1] <- "HSP"

# Regularized ensemble weighting on the test set <glmnet>
# 10-fold CV
ens <- trainControl(method = "cv", number = 10, classProbs = T, summaryFunction = twoClassSummary)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.ens <- train(CRP ~ CRPglm + CRPrf + CRPgbm + CRPnn, data = crpens,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "ROC",
                 trControl = ens)
CRP.ens
crp.pred <- predict(CRP.ens, crpens, type="prob")
crp.test <- cbind(crpens, crp.pred)
crp <- subset(crp.test, CRP=="Y", select=c(Y))
cra <- subset(crp.test, CRP=="N", select=c(Y))
crp.eval <- evaluate(p=crp[,1], a=cra[,1]) ## calculate ROC's on test set <dismo>
crp.eval
plot(crp.eval, 'ROC') ## plot ROC curve
crp.thld <- threshold(crp.eval, 'spec_sens') ## TPR+TNR threshold for classification
crpens.pred <- predict(pred, CRP.ens, type="prob") ## spatial prediction
plot(1-crpens.pred, axes = F)
crpmask <- 1-crpens.pred > crp.thld
plot(crpmask, axes = F, legend = F)

# presence/absence of Buildings/Rural Settlements (HSP, present = Y, absent = N)
HSP.ens <- train(HSP ~ HSPglm + HSPrf + HSPgbm + HSPnn, data = hspens,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "ROC",
                 trControl = ens)
HSP.ens
hsp.pred <- predict(HSP.ens, hspens, type="prob")
hsp.test <- cbind(hspens, hsp.pred)
hsp <- subset(hsp.test, HSP=="Y", select=c(Y))
hsa <- subset(hsp.test, HSP=="N", select=c(Y))
hsp.eval <- evaluate(p=hsp[,1], a=hsa[,1]) ## calculate ROC's on test set <dismo>
hsp.eval
plot(hsp.eval, 'ROC') ## plot ROC curve
hsp.thld <- threshold(hsp.eval, 'spec_sens') ## TPR+TNR threshold for classification
hspens.pred <- predict(pred, HSP.ens, type="prob") ## spatial prediction
plot(1-hspens.pred, axes = F)
hspmask <- 1-hspens.pred > hsp.thld
plot(hspmask, axes = F, legend = F)

#+ Write spatial predictions -----------------------------------------------
# Create a "Results" folder in current working directory
dir.create("TZ_1MGS_results", showWarnings=F)

# Export Gtif's to "./TZ_results"
writeRaster(crp.preds, filename="./TZ_1MGS_results/TZ_1MGS_crpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(hsp.preds, filename="./TZ_1MGS_results/TZ_1MGS_hspreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
# Ensemble predictions
enspred <- stack(1-crpens.pred, crpmask, 1-hspens.pred, hspmask)
writeRaster(enspred, filename="./TZ_1MGS_results/TZ_1MGS_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
