#' Ensemble predictions of Tanzania GeoSurvey cropland,
#' woody vegetation cover and human settlement observations. 
#' M. Walsh, November 2014

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

#+ Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("TZ_data", showWarnings=F)
dat_dir <- "./TZ_data"

# download GeoSurvey data
download("https://www.dropbox.com/s/gfgjnrgllwqt79d/TZ_geos_123114.csv?dl=0", "./TZ_data/TZ_geos_123114.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/TZ_geos_123114.csv", sep=""), header=T, sep=",")
geos <- na.omit(geos)

# download Tanzania Gtifs (~27.9 Mb) and stack in raster
download("https://www.dropbox.com/s/otiqe78s0kf1z1s/TZ_grids.zip?dl=0", "./TZ_data/TZ_grids.zip", mode="wb")
unzip("./TZ_data/TZ_grids.zip", exdir="./TZ_data", overwrite=T)
glist <- list.files(path="./TZ_data", pattern="tif", full.names=T)
grid <- stack(glist)

#+ Data setup --------------------------------------------------------------
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

#+ Split data into train and test sets ------------------------------------
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

#+ Stepwise main effects GLM's <MASS> --------------------------------------
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

#+ Random forests <randomForest> -------------------------------------------
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

#+ Gradient boosting <gbm> ------------------------------------------
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

#+ Neural nets <nnet> ------------------------------------------------------
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

#+ Plot predictions by GeoSurvey variables ---------------------------------
# Cropland prediction plots
crp.preds <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred)
names(crp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(crp.preds, axes = F)

# Woody vegetation cover >60% prediction plots
wcp.preds <- stack(1-wcpglm.pred, 1-wcprf.pred, 1-wcpgbm.pred, 1-wcpnn.pred)
names(wcp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(wcp.preds, axes = F)

# Settlement prediction plots
hsp.preds <- stack(1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(hsp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(hsp.preds, axes = F)

#+ Ensemble predictions <glm>, <rf>, <gbm>, <nnet> --------------------------
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

# Regularized ensemble weighting on the test set <glmnet>
# 10-fold CV
ens <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.ens <- train(CRP ~ CRPglm + CRPrf + CRPgbm + CRPnn, data = crpensTest,
                 family = "binomial", 
                 method = "glmnet",
                 trControl = ens)
crp.pred <- predict(CRP.ens, crpensTest, type="prob")
crp.test <- cbind(crpensTest, crp.pred)
crp <- subset(crp.test, CRP=="Y", select=c(Y))
cra <- subset(crp.test, CRP=="N", select=c(Y))
crp.eval <- evaluate(p=crp[,1], a=cra[,1]) ## calculate ROC's on test set <dismo>
crp.eval
plot(crp.eval, 'ROC') ## plot ROC curve
crp.thld <- threshold(crp.eval, 'spec_sens') ## TPR+TNR threshold for classification
crpens.pred <- predict(pred, CRP.ens, type="prob") ## spatial prediction
crpmask <- 1-crpens.pred > crp.thld
plot(crpmask, axes = F, legend = F)

# presence/absence of Woody Vegetation Cover >60% (WCP, present = Y, absent = N)
WCP.ens <- train(WCP ~ WCPglm + WCPrf + WCPgbm + WCPnn, data = wcpensTest,
                 family = "binomial", 
                 method = "glmnet",
                 trControl = ens)
wcp.pred <- predict(WCP.ens, wcpensTest, type="prob")
wcp.test <- cbind(wcpensTest, wcp.pred)
wcp <- subset(wcp.test, WCP=="Y", select=c(Y))
wca <- subset(wcp.test, WCP=="N", select=c(Y))
wcp.eval <- evaluate(p=wcp[,1], a=wca[,1]) ## calculate ROC's on test set <dismo>
wcp.eval
plot(wcp.eval, 'ROC') ## plot ROC curve
wcp.thld <- threshold(wcp.eval, 'spec_sens') ## TPR+TNR threshold for classification
wcpens.pred <- predict(pred, WCP.ens, type="prob") ## spatial prediction
wcpmask <- 1-wcpens.pred > wcp.thld
plot(wcpmask, axes = F, legend = F)

# presence/absence of Buildings/Rural Settlements (HSP, present = Y, absent = N)
HSP.ens <- train(HSP ~ HSPglm + HSPrf + HSPgbm + HSPnn, data = hspensTest,
                 family = "binomial", 
                 method = "glmnet",
                 trControl = ens)
hsp.pred <- predict(HSP.ens, hspensTest, type="prob")
hsp.test <- cbind(hspensTest, hsp.pred)
hsp <- subset(hsp.test, HSP=="Y", select=c(Y))
hsa <- subset(hsp.test, HSP=="N", select=c(Y))
hsp.eval <- evaluate(p=hsp[,1], a=hsa[,1]) ## calculate ROC's on test set <dismo>
hsp.eval
plot(hsp.eval, 'ROC') ## plot ROC curve
hsp.thld <- threshold(hsp.eval, 'spec_sens') ## TPR+TNR threshold for classification
hspens.pred <- predict(pred, HSP.ens, type="prob") ## spatial prediction
hspmask <- 1-hspens.pred > hsp.thld
plot(hspmask, axes = F, legend = F)

#+ Write spatial predictions -----------------------------------------------
# Create a "Results" folder in current working directory
dir.create("TZ_results", showWarnings=F)

# Export Gtif's to "./TZ_results"
writeRaster(crp.preds, filename="./TZ_results/TZ_crpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(wcp.preds, filename="./TZ_results/TZ_wcpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(hsp.preds, filename="./TZ_results/TZ_hspreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
# Ensemble predictions
enspred <- stack(1-crpens.pred, crpmask, 1-wcpens.pred, wcpmask, 1-hspens.pred, hspmask)
writeRaster(enspred, filename="./TZ_results/TZ_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
