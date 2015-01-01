# Validation of ensemble predictions of Tanzania 1M GeoSurvey cropland and
# human settlement observations.
# M.Walsh, J.Chen & A.Verlinden, January 2015

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","MASS","randomForest","gbm","nnet","glmnet","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(MASS)
require(randomForest)
require(gbm)
require(nnet)
require(glmnet)
require(dismo)

# Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("TZ_1MGS_data", showWarnings=F)
dat_dir <- "./TZ_1MGS_data"

# download 1M GeoSurvey data
download("https://www.dropbox.com/s/eq19mgnj86d4qo2/1MGS_123114.csv?dl=0", "./TZ_1MGS_data/1MGS_123114.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/1MGS_123114.csv", sep=""), header=T, sep=",")

# download Tanzania validation data
download("https://www.dropbox.com/s/gfgjnrgllwqt79d/TZ_geos_123114.csv?dl=0", "./TZ_1MGS_data/TZ_geos_123114.csv", mode="wb")
geosv <- read.table(paste(dat_dir, "/TZ_geos_123114.csv", sep=""), header=T, sep=",")

# download Tanzania Gtifs (~27.9 Mb) and stack in raster
download("https://www.dropbox.com/s/otiqe78s0kf1z1s/TZ_grids.zip?dl=0", "./TZ_1MGS_data/TZ_grids.zip", mode="wb")
unzip("./TZ_1MGS_data/TZ_grids.zip", exdir="./TZ_1MGS_data", overwrite=T)
glist <- list.files(path="./TZ_1MGS_data", pattern="tif", full.names=T)
grid <- stack(glist)

# Data setup --------------------------------------------------------------
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

# Extract gridded variables at GeoSurvey locations
geosgrid <- extract(grid, geos)
geosvgrid <- extract(grid, geosv)

# Assemble dataframes
# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP <- geos$CRP
crpdat <- cbind.data.frame(CRP, geosgrid)
crpdat <- na.omit(crpdat)

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
# note that this excludes large urban areas where MODIS fPAR = 0
HSP <- geos$HSP
hspdat <- cbind.data.frame(HSP, geosgrid)
hspdat <- na.omit(hspdat)

# Stepwise main effects GLM's <MASS> --------------------------------------
# 10-fold CV
step <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.glm <- train(CRP ~ ., data = crpdat,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = step)
crpglm.test <- predict(CRP.glm, crpdat) ## predict test-set
confusionMatrix(crpglm.test, crpdat$CRP, "Y") ## print validation summaries
crpglm.pred <- predict(grid, CRP.glm, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.glm <- train(HSP ~ ., data = hspdat,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
hspglm.test <- predict(HSP.glm, hspdat) ## predict test-set
confusionMatrix(hspglm.test, hspdat$HSP, "Y") ## print validation summaries
hspglm.pred <- predict(grid, HSP.glm, type = "prob") ## spatial predictions

# Random forests <randomForest> -------------------------------------------
# out-of-bag predictions
oob <- trainControl(method = "oob")

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.rf <- train(CRP ~ ., data = crpdat,
                method = "rf",
                trControl = oob)
crprf.test <- predict(CRP.rf, crpdat) ## predict test-set
confusionMatrix(crprf.test, crpdat$CRP, "Y") ## print validation summaries
crprf.pred <- predict(grid, CRP.rf, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.rf <- train(HSP ~ ., data = hspdat,
                method = "rf",
                trControl = oob)
hsprf.test <- predict(HSP.rf, hspdat) ## predict test-set
confusionMatrix(hsprf.test, hspdat$HSP, "Y") ## print validation summaries
hsprf.pred <- predict(grid, HSP.rf, type = "prob") ## spatial predictions

# Gradient boosting <gbm> ------------------------------------------
# CV for training gbm's
gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.gbm <- train(CRP ~ ., data = crpdat,
                 method = "gbm",
                 trControl = gbm)
crpgbm.test <- predict(CRP.gbm, crpdat) ## predict test-set
confusionMatrix(crpgbm.test, crpdat$CRP, "Y") ## print validation summaries
crpgbm.pred <- predict(grid, CRP.gbm, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.gbm <- train(HSP ~ ., data = hspdat,
                 method = "gbm",
                 trControl = gbm)
hspgbm.test <- predict(HSP.gbm, hspdat) ## predict test-set
confusionMatrix(hspgbm.test, hspdat$HSP, "Y") ## print validation summaries
hspgbm.pred <- predict(grid, HSP.gbm, type = "prob") ## spatial predictions

# Neural nets <nnet> ------------------------------------------------------
# CV for training nnet's
nn <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.nn <- train(CRP ~ ., data = crpdat,
                method = "nnet",
                trControl = nn)
crpnn.test <- predict(CRP.nn, crpdat) ## predict test-set
confusionMatrix(crpnn.test, crpdat$CRP, "Y") ## print validation summaries
crpnn.pred <- predict(grid, CRP.nn, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.nn <- train(HSP ~ ., data = hspdat,
                method = "nnet",
                trControl = nn)
hspnn.test <- predict(HSP.nn, hspdat) ## predict test-set
confusionMatrix(hspnn.test, hspdat$HSP, "Y") ## print validation summaries
hspnn.pred <- predict(grid, HSP.nn, type = "prob") ## spatial predictions

# Plot predictions by GeoSurvey variables ---------------------------------
# Cropland prediction plots
crp.preds <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred)
names(crp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(crp.preds, axes = F)

# Settlement prediction plots
hsp.preds <- stack(1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(hsp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(hsp.preds, axes = F)

# Ensemble predictions <glm>, <rf>, <gbm>, <nnet> --------------------------
# Ensemble set-up
pred <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred,
              1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(pred) <- c("CRPglm","CRPrf","CRPgbm","CRPnn",
                 "HSPglm","HSPrf","HSPgbm","HSPnn")
geosvpred <- extract(pred, geosv)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRPv <- geosv$CRP
crpens <- cbind.data.frame(CRPv, geosvpred)
crpens <- na.omit(crpens)

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSPv <- geosv$HSP
hspens <- cbind.data.frame(HSPv, geosvpred)
hspens <- na.omit(hspens)

# Regularized ensemble weighting on the test set <glmnet>
# 10-fold CV
ens <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.ens <- train(CRPv ~ CRPglm + CRPrf + CRPgbm + CRPnn, data = crpens,
                 family = "binomial", 
                 method = "glmnet",
                 trControl = ens)
crp.pred <- predict(CRP.ens, crpens, type="prob")
crp.test <- cbind(crpens, crp.pred)
crp <- subset(crp.test, CRPv=="Y", select=c(Y))
cra <- subset(crp.test, CRPv=="N", select=c(Y))
crp.eval <- evaluate(p=crp[,1], a=cra[,1]) ## calculate ROC's on test set <dismo>
crp.eval
plot(crp.eval, 'ROC') ## plot ROC curve
crp.thld <- threshold(crp.eval, 'spec_sens') ## TPR+TNR threshold for classification
crpens.pred <- predict(pred, CRP.ens, type="prob") ## spatial prediction
plot(1-crpens.pred, axes = F)
crpmask <- 1-crpens.pred > crp.thld
plot(crpmask, axes = F, legend = F)

# presence/absence of Buildings/Rural Settlements (HSP, present = Y, absent = N)
HSP.ens <- train(HSPv ~ HSPglm + HSPrf + HSPgbm + HSPnn, data = hspens,
                 family = "binomial", 
                 method = "glmnet",
                 trControl = ens)
hsp.pred <- predict(HSP.ens, hspens, type="prob")
hsp.test <- cbind(hspens, hsp.pred)
hsp <- subset(hsp.test, HSPv=="Y", select=c(Y))
hsa <- subset(hsp.test, HSPv=="N", select=c(Y))
hsp.eval <- evaluate(p=hsp[,1], a=hsa[,1]) ## calculate ROC's on test set <dismo>
hsp.eval
plot(hsp.eval, 'ROC') ## plot ROC curve
hsp.thld <- threshold(hsp.eval, 'spec_sens') ## TPR+TNR threshold for classification
hspens.pred <- predict(pred, HSP.ens, type="prob") ## spatial prediction
plot(1-hspens.pred, axes = F)
hspmask <- 1-hspens.pred > hsp.thld
plot(hspmask, axes = F, legend = F)

# Write spatial predictions -----------------------------------------------
# Create a "Results" folder in current working directory
dir.create("TZ_1MGS_results", showWarnings=F)

# Export Gtif's to "./TZ_results"
writeRaster(crp.preds, filename="./TZ_1MGS_results/TZ_1MGS_crpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(hsp.preds, filename="./TZ_1MGS_results/TZ_1MGS_hspreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
# Ensemble predictions
enspred <- stack(1-crpens.pred, crpmask, 1-hspens.pred, hspmask)
writeRaster(enspred, filename="./TZ_1MGS_results/TZ_1MGS_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)


