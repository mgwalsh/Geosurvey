# Ensemble machine learning predictions of Tanzania GeoSurvey cropland,
# woody vegetation cover and human settlement observations. 
# M. Walsh, November 2014

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","MASS","randomForest","gbm","nnet","ROCR)), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(MASS)
require(randomForest)
require(gbm)
require(nnet)
require(ROCR)

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
# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP <- geos$CRP
crpdat <- cbind.data.frame(CRP, geosgrid)
crpdat <- na.omit(crpdat)

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
WCP <- geos$WCP
wcpdat <- cbind.data.frame(WCP, geosgrid)
wcpdat <- na.omit(wcpdat)

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
# note that this excludes large urban areas where MODIS fPAR = 0
HSP <- geos$HSP
hspdat <- cbind.data.frame(HSP, geosgrid)
hspdat <- na.omit(hspdat)

# set train/test set randomization seed
seed <- 1385321
set.seed(seed)

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

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.glm <- train(CRP ~ ., data = crpTrain,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = step)
crpglm.test <- predict(CRP.glm, crpTest) ## predict test-set
confusionMatrix(crpglm.test, crpTest$CRP, "Y") ## print validation summaries
crpglm.pred <- predict(grid, CRP.glm, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
WCP.glm <- train(WCP ~ ., data = wcpTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
wcpglm.test <- predict(WCP.glm, wcpTest) ## predict test-set
confusionMatrix(wcpglm.test, wcpTest$WCP, "Y") ## print validation summaries
wcpglm.pred <- predict(grid, WCP.glm, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.glm <- train(HSP ~ ., data = hspTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
hspglm.test <- predict(HSP.glm, hspTest) ## predict test-set
confusionMatrix(hspglm.test, hspTest$HSP, "Y") ## print validation summaries
hspglm.pred <- predict(grid, HSP.glm, type = "prob") ## spatial predictions

# Plot <MASS> predictions
glmpreds <- stack(1-crpglm.pred, 1-wcpglm.pred, 1-hspglm.pred)
names(glmpreds) <- c("CRPglm", "WCPglm", "HSPglm")
plot(glmpreds, axes = F)

# Random forests <randomForest> -------------------------------------------
# out-of-bag predictions
oob <- trainControl(method = "oob")

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.rf <- train(CRP ~ ., data = crpTrain,
                method = "rf",
                trControl = oob)
crprf.test <- predict(CRP.rf, crpTest) ## predict test-set
confusionMatrix(crprf.test, crpTest$CRP, "Y") ## print validation summaries
crprf.pred <- predict(grid, CRP.rf, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
WCP.rf <- train(WCP ~ ., data = wcpTrain,
                method = "rf",
                trControl = oob)
wcprf.test <- predict(WCP.rf, wcpTest) ## predict test-set
confusionMatrix(wcprf.test, wcpTest$WCP, "Y") ## print validation summaries
wcprf.pred <- predict(grid, WCP.rf, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.rf <- train(HSP ~ ., data = hspTrain,
                method = "rf",
                trControl = oob)
hsprf.test <- predict(HSP.rf, hspTest) ## predict test-set
confusionMatrix(hsprf.test, hspTest$HSP, "Y") ## print validation summaries
hsprf.pred <- predict(grid, HSP.rf, type = "prob") ## spatial predictions

# Plot <randomForest> predictions
rfpreds <- stack(1-crprf.pred, 1-wcprf.pred, 1-hsprf.pred)
names(rfpreds) <- c("CRPrf", "WCPrf", "HSPrf")
plot(rfpreds, axes = F)

# Gradient boosting <gbm> ------------------------------------------
# CV for training gbm's
gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.gbm <- train(CRP ~ ., data = crpTrain,
                 method = "gbm",
                 trControl = gbm)
crpgbm.test <- predict(CRP.gbm, crpTest) ## predict test-set
confusionMatrix(crpgbm.test, crpTest$CRP, "Y") ## print validation summaries
crpgbm.pred <- predict(grid, CRP.gbm, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
WCP.gbm <- train(WCP ~ ., data = wcpTrain,
                 method = "gbm",
                 trControl = gbm)
wcpgbm.test <- predict(WCP.gbm, wcpTest) ## predict test-set
confusionMatrix(wcpgbm.test, wcpTest$WCP, "Y") ## print validation summaries
wcpgbm.pred <- predict(grid, WCP.gbm, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.gbm <- train(HSP ~ ., data = hspTrain,
                 method = "gbm",
                 trControl = gbm)
hspgbm.test <- predict(HSP.gbm, hspTest) ## predict test-set
confusionMatrix(hspgbm.test, hspTest$HSP, "Y") ## print validation summaries
hspgbm.pred <- predict(grid, HSP.gbm, type = "prob") ## spatial predictions

# Plot <gbm> predictions
gbmpreds <- stack(1-crpgbm.pred, 1-wcpgbm.pred, 1-hspgbm.pred)
names(gbmpreds) <- c("CRPgbm", "WCPgbm", "HSPgbm")
plot(gbmpreds, axes = F)

# Neural nets <nnet> ------------------------------------------------------
# CV for training nnet's
nn <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.nn <- train(CRP ~ ., data = crpTrain,
                method = "nnet",
                trControl = nn)
crpnn.test <- predict(CRP.nn, crpTest) ## predict test-set
confusionMatrix(crpnn.test, crpTest$CRP, "Y") ## print validation summaries
crpnn.pred <- predict(grid, CRP.nn, type = "prob") ## spatial predictions

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
WCP.nn <- train(WCP ~ ., data = wcpTrain,
                method = "nnet",
                trControl = nn)
wcpnn.test <- predict(WCP.nn, wcpTest) ## predict test-set
confusionMatrix(wcpnn.test, wcpTest$WCP, "Y") ## print validation summaries
wcpnn.pred <- predict(grid, WCP.nn, type = "prob") ## spatial predictions

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.nn <- train(HSP ~ ., data = hspTrain,
                method = "nnet",
                trControl = nn)
hspnn.test <- predict(HSP.nn, hspTest) ## predict test-set
confusionMatrix(hspnn.test, hspTest$HSP, "Y") ## print validation summaries
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
# Ensemble set-up
pred <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred,
              1-wcpglm.pred, 1-wcprf.pred, 1-wcpgbm.pred, 1-wcpnn.pred,
              1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(pred) <- c("CRPglm","CRPrf","CRPgbm","CRPnn",
                 "WCPglm","WCPrf","WCPgbm","WCPnn",
                 "HSPglm","HSPrf","HSPgbm","HSPnn")
geospred <- extract(pred, geos)

# presence/absence of Cropland (CRP, present = Y, absent = N)
crpens <- cbind.data.frame(CRP, geospred)
crpens <- na.omit(crpens)
crpensTest <- crpens[-crpIndex,] ## replicate previous test set

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
wcpens <- cbind.data.frame(WCP, geospred)
wcpens <- na.omit(wcpens)
wcpensTest <- wcpens[-wcpIndex,] ## replicate previous test set

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
hspens <- cbind.data.frame(HSP, geospred)
hspens <- na.omit(hspens)
hspensTest <- hspens[-hspIndex,] ## replicate previous test set

# GLM-based ensemble weighting on test set
# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.ens <- glm(CRP ~ CRPglm + CRPrf + CRPgbm + CRPnn, data=crpensTest,
               family = binomial(link="logit"))
summary(CRP.ens)
crpens.pred <- predict(pred, CRP.ens, type="response")
plot(crpens.pred, axes = F)

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
WCP.ens <- glm(WCP ~ WCPglm + WCPrf + WCPgbm + WCPnn, data=wcpensTest,
               family = binomial(link="logit"))
summary(WCP.ens)
wcpens.pred <- predict(pred, WCP.ens, type="response")
plot(wcpens.pred, axes = F)

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.ens <- glm(HSP ~ HSPglm + HSPrf + HSPgbm + HSPnn, data=hspensTest,
               family = binomial(link="logit"))
summary(HSP.ens)
hspens.pred <- predict(pred, HSP.ens, type="response")
plot(hspens.pred, axes = F)

# Plot ensemble predictions
enspred <- stack(crpens.pred, wcpens.pred, hspens.pred)
names(enspred) <- c("CRP", "WCP", "HSP")
plot(enspred, axes = F, main = "")

# Receiver/Operator curves of ensemble predictions on test set -------------
# Cropland ensemble predictions
crpprob <- predict(CRP.ens, crpensTest, type="response")
crppred <- prediction(crpprob, crpensTest$CRP)
crproc <- performance(crppred, "tpr", "fpr")
plot(crproc)
crpsens <- performance(crppred, "sens")
crpspec <- performance(crppred, "spec")
plot(crpsens, xlab = "p(CRP = Y)", col="blue", ylab = "Sensitivity & Specificity")
plot(crpspec, col="red", add = T)

# Woody vegetation cover >60% ensemble predictions
wcpprob <- predict(WCP.ens, wcpensTest, type="response")
wcppred <- prediction(wcpprob, wcpensTest$WCP)
wcproc <- performance(wcppred, "tpr", "fpr")
plot(wcproc)
wcpsens <- performance(wcppred, "sens")
wcpspec <- performance(wcppred, "spec")
plot(wcpsens, xlab = "p(WCP = Y)", col="blue", ylab = "Sensitivity & Specificity")
plot(wcpspec, col="red", add = T)

# Building/Human settlement ensemble predictions
hspprob <- predict(HSP.ens, hspensTest, type="response")
hsppred <- prediction(hspprob, hspensTest$HSP)
hsproc <- performance(hsppred, "tpr", "fpr")
plot(hsproc)
hspsens <- performance(hsppred, "sens")
hspspec <- performance(hsppred, "spec")
plot(hspsens, xlab = "p(HSP = Y)", col="blue", ylab = "Sensitivity & Specificity")
plot(hspspec, col="red", add = T)

# Write spatial predictions -----------------------------------------------
# Create a "Results" folder in your current working directory
dir.create("Results", showWarnings=F)

# Export Gtif's to "./Results"
writeRaster(crp.preds, filename="./Results/TZ_crpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(wcp.preds, filename="./Results/TZ_wcpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(hsp.preds, filename="./Results/TZ_hspreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(enspred, filename="./Results/TZ_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
