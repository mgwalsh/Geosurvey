# Ensemble machine learning predictions of Tanzania GeoSurvey cropland,
# woody vegetation cover and human settlement observations. 
# M. Walsh, November 2014

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","MASS",randomForest","gbm","nnet")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(MASS)
require(randomForest)
require(gbm)
require(nnet)

# Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("Data", showWarnings=F)
dat_dir <- "./Data"

# download GeoSurvey data
download("https://www.dropbox.com/s/03l4m4zjdi5mhyu/TZ_geos_1114.csv?dl=0", "./Data/TZ_geos_1114.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/TZ_geos_1114.csv", sep=""), header=T, sep=",")
geos <- na.omit(geos)

# download Tanzania Gtifs (~27 Mb) and stack in raster
download("https://www.dropbox.com/s/otiqe78s0kf1z1s/TZ_grids.zip?dl=0", "./Data/TZ_grids.zip", mode="wb")
unzip("./Data/TZ_grids.zip", exdir="./Data", overwrite=T)
glist <- list.files(path="./Data", pattern="tif", full.names=T)
grid <- stack(glist)

# Data setup --------------------------------------------------------------
# Project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grid)

# Extract gridded variables at GeoSurvey locations
geosgrid <- extract(grid, geos)

# Assemble dataframes
# presence/absence of Cropland (CRP, present = Yes, absent = No)
CRP <- geos$CRP
crpdat <- cbind.data.frame(CRP, geosgrid)
crpdat <- na.omit(crpdat)

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Yes, absent = No)
WCP <- geos$WCP
wcpdat <- cbind.data.frame(WCP, geosgrid)
wcpdat <- na.omit(wcpdat)

# presence/absence of Buildings/Human Settlements (HSP, present = Yes, absent = No)
# note that this excludes large urban areas where MODIS fPAR = 0
HSP <- geos$HSP
hspdat <- cbind.data.frame(HSP, geosgrid)
hspdat <- na.omit(hspdat)

# set train/test set randomization seed
set.seed(1385321)

# Split data into train and test sets ------------------------------------
# Cropland train/test split
crpIndex <- createDataPartition(crpdat$CRP, p = 0.75, list = FALSE, times = 1)
crpTrain <- crpdat[ crpIndex,]
crpTest  <- crpdat[-crpIndex,]

# Woody cover train/test split
wcpIndex <- createDataPartition(wcpdat$WCP, p = 0.75, list = FALSE, times = 1)
wcpTrain <- wcpdat[ wcpIndex,]
wcpTest  <- wcpdat[-wcpIndex,]

# Settlement train/test split
hspIndex <- createDataPartition(hspdat$HSP, p = 0.75, list = FALSE, times = 1)
hspTrain <- hspdat[ hspIndex,]
hspTest  <- hspdat[-hspIndex,]

# Stepwise main effects GLM's <MASS> --------------------------------------
# 10-fold CV
step <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Yes, absent = No)
CRP.glm <- train(CRP ~ ., data = crpTrain,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = step)
crpglm.test <- predict(CRP.glm, crpTest) ## predict test-set
confusionMatrix(crpglm.test, crpTest$CRP) ## print validation summaries
crpglm.pred <- predict(grid, CRP.glm, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Yes, absent = No)
WCP.glm <- train(WCP ~ ., data = wcpTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
wcpglm.test <- predict(WCP.glm, wcpTest) ## predict test-set
confusionMatrix(wcpglm.test, wcpTest$WCP) ## print validation summaries
wcpglm.pred <- predict(grid, WCP.glm, type = "prob") ## spatial predictions

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = Yes, absent = No)
HSP.glm <- train(HSP ~ ., data = hspTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
hspglm.test <- predict(HSP.glm, hspTest) ## predict test-set
confusionMatrix(hspglm.test, hspTest$HSP) ## print validation summaries
hspglm.pred <- predict(grid, HSP.glm, type = "prob") ## spatial predictions

# Plot <MASS> predictions
glmpreds <- stack(1-crpglm.pred, 1-wcpglm.pred, 1-hspglm.pred)
names(glmpreds) <- c("CRPglm", "WCPglm", "HSPglm")
plot(glmpreds, axes = F)

# Random forests <randomForest> -------------------------------------------
# out-of-bag predictions
oob <- trainControl(method = "oob")

# presence/absence of Cropland (CRP, present = Yes, absent = No)
CRP.rf <- train(CRP ~ ., data = crpTrain,
                method = "rf",
                trControl = oob)
crprf.test <- predict(CRP.rf, crpTest) ## predict test-set
confusionMatrix(crprf.test, crpTest$CRP) ## print validation summaries
crprf.pred <- predict(grid, CRP.rf, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Yes, absent = No)
WCP.rf <- train(WCP ~ ., data = wcpTrain,
                method = "rf",
                trControl = oob)
wcprf.test <- predict(WCP.rf, wcpTest) ## predict test-set
confusionMatrix(wcprf.test, wcpTest$WCP) ## print validation summaries
wcprf.pred <- predict(grid, WCP.rf, type = "prob") ## spatial predictions

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = Yes, absent = No)
HSP.rf <- train(HSP ~ ., data = hspTrain,
                method = "rf",
                trControl = oob)
hsprf.test <- predict(HSP.rf, hspTest) ## predict test-set
confusionMatrix(hsprf.test, hspTest$HSP) ## print validation summaries
hsprf.pred <- predict(grid, HSP.rf, type = "prob") ## spatial predictions

# Plot <randomForest> predictions
rfpreds <- stack(1-crprf.pred, 1-wcprf.pred, 1-hsprf.pred)
names(rfpreds) <- c("CRPrf", "WCPrf", "HSPrf")
plot(rfpreds, axes = F)

# Generalized boosting <gbm> ------------------------------------
# CV for training gbm's
gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# presence/absence of Cropland (CRP, present = Yes, absent = No)
CRP.gbm <- train(CRP ~ ., data = crpTrain,
                 method = "gbm",
                 trControl = gbm)
crpgbm.test <- predict(CRP.gbm, crpTest) ## predict test-set
confusionMatrix(crpgbm.test, crpTest$CRP) ## print validation summaries
crpgbm.pred <- predict(grid, CRP.gbm, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Yes, absent = No)
WCP.gbm <- train(WCP ~ ., data = wcpTrain,
                 method = "gbm",
                 trControl = gbm)
wcpgbm.test <- predict(WCP.gbm, wcpTest) ## predict test-set
confusionMatrix(wcpgbm.test, wcpTest$WCP) ## print validation summaries
wcpgbm.pred <- predict(grid, WCP.gbm, type = "prob") ## spatial predictions

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = Yes, absent = No)
HSP.gbm <- train(HSP ~ ., data = hspTrain,
                 method = "gbm",
                 trControl = gbm)
hspgbm.test <- predict(HSP.gbm, hspTest) ## predict test-set
confusionMatrix(hspgbm.test, hspTest$HSP) ## print validation summaries
hspgbm.pred <- predict(grid, HSP.gbm, type = "prob") ## spatial predictions

# Plot <gbm> predictions
gbmpreds <- stack(1-crpgbm.pred, 1-wcpgbm.pred, 1-hspgbm.pred)
names(gbmpreds) <- c("CRPgbm", "WCPgbm", "HSPgbm")
plot(gbmpreds, axes = F)

# Neural nets <nnet> ------------------------------------------------------
# CV for training nnet's
nn <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Yes, absent = No)
CRP.nn <- train(CRP ~ ., data = crpTrain,
                method = "nnet",
                trControl = nn)
crpnn.test <- predict(CRP.nn, crpTest) ## predict test-set
confusionMatrix(crpnn.test, crpTest$CRP) ## print validation summaries
crpnn.pred <- predict(grid, CRP.nn, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Yes, absent = No)
WCP.nn <- train(WCP ~ ., data = wcpTrain,
                method = "nnet",
                trControl = nn)
wcpnn.test <- predict(WCP.nn, wcpTest) ## predict test-set
confusionMatrix(wcpnn.test, wcpTest$WCP) ## print validation summaries
wcpnn.pred <- predict(grid, WCP.nn, type = "prob") ## spatial predictions

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = Yes, absent = No)
HSP.nn <- train(HSP ~ ., data = hspTrain,
                method = "nnet",
                trControl = nn)
hspnn.test <- predict(HSP.nn, hspTest) ## predict test-set
confusionMatrix(hspnn.test, hspTest$HSP) ## print validation summaries
hspnn.pred <- predict(grid, HSP.nn, type = "prob") ## spatial predictions

# Plot <nnet> predictions
nnpreds <- stack(1-crpnn.pred, 1-wcpnn.pred, 1-hspnn.pred)
names(nnpreds) <- c("CRPnn", "WCPnn", "HSPnn")
plot(nnpreds, axes = F)

# Plot predictions by GeoSurvey variables ---------------------------------
# Cropland prediction plots
crp.preds <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred)
names(crp.preds) <- c("glm","randomForest","gbm","nnet")
plot(crp.preds, axes = F)

# Woody vegetation cover >60% prediction plots
wcp.preds <- stack(1-wcpglm.pred, 1-wcprf.pred, 1-wcpgbm.pred, 1-wcpnn.pred)
names(wcp.preds) <- c("glm","randomForest","gbm","nnet")
plot(wcp.preds, axes = F)

# Human settlement prediction plots
hsp.preds <- stack(1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(hsp.preds) <- c("glm","randomForest","gbm","nnet")
plot(hsp.preds, axes = F)

# Ensemble predictions <glm>, <rf>, <gbm>, <nnet> --------------------------
# Ensemble set up
preds <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred,
               1-wcpglm.pred, 1-wcprf.pred, 1-wcpgbm.pred, 1-wcpnn.pred,
               1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(preds) <- c("CRPglm","CRPrf","CRPgbm","CRPnn",
                  "WCPglm","WCPrf","WCPgbm","WCPnn",
                  "HSPglm","HSPrf","HSPgbm","HSPnn")
expreds <- extract(preds, geos)
ensdat <- data.frame(geos[,4:6], expreds)
ensdat <- na.omit(ensdat)
ensIndex <- createDataPartition(ensdat$CRP, p = 0.75, list = FALSE, times = 1)
ensTrain <- ensdat[ensIndex,] ## replicate previous training set with prediction rasters
ensTest  <- ensdat[-ensIndex,] ## replicate previous test set with prediction rasters

# GLM weighted ensemble predictions for individual <glm>, <randomForest>, <gbm>
# & <nnet> model train sets
# presence/absence of Cropland (CRP, present = Yes, absent = N)
CRP.ens <- train(CRP ~ CRPglm + CRPrf + CRPgbm + CRPnn, data = ensTrain,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = step)
summary(CRP.ens) 
crpens.test <- predict(CRP.ens, ensTest) ## predict test-set
confusionMatrix(crpens.test, ensTest$CRP) ## print test set performance
crpens.pred <- predict(preds, CRP.ens, type="prob") ## predict grid
plot(1-crpens.pred, axes = F) ## plot gridded ensemble predictions

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Yes, absent = No)
WCP.ens <- train(WCP ~ WCPglm + WCPrf + WCPgbm + WCPnn, data = ensTrain,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = step)
summary(WCP.ens) 
wcpens.test <- predict(WCP.ens, ensTest) ## predict test-set
confusionMatrix(wcpens.test, ensTest$WCP) ## print test set performance
wcpens.pred <- predict(preds, WCP.ens, type="prob") ## predict grid
plot(1-wcpens.pred, axes = F) ## plot gridded ensemble predictions

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = Yes, absent = No)
HSP.ens <- train(HSP ~ HSPglm + HSPrf + HSPgbm + HSPnn, data = ensTrain,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = step)
summary(HSP.ens) 
hspens.test <- predict(HSP.ens, ensTest) ## predict test-set
confusionMatrix(hspens.test, ensTest$HSP) ## print test set performance
hspens.pred <- predict(preds, HSP.ens, type="prob") ## predict grid
plot(1-hspens.pred, axes = F) ## plot gridded ensemble predictions

# Plot ensemble predictions
enspred <- stack(1-crpens.pred, 1-wcpens.pred, 1-hspens.pred)
names(enspred) <- c("CRP", "WCP", "HSP")
plot(enspred, axes = F, main = "")

# Write spatial predictions -----------------------------------------------
# Create a "Results" folder in your current working directory
dir.create("Results", showWarnings=F)

# Export Gtif's to "./Results"
writeRaster(crp.preds, filename="./Results/TZ_crpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(wcp.preds, filename="./Results/TZ_wcpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(hsp.preds, filename="./Results/TZ_hspreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(enspred, filename="./Results/TZ_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
