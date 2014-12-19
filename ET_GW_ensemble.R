# Ensemble predictions of Ethiopia Geo-Wiki cropland,
# and human settlement observations.
# Cropland & Human Settlement survey data courtesy http://www.geo-wiki.org/download-data
# M. Walsh, December 2014

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","MASS","randomForest","gbm","nnet","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(MASS)
require(randomForest)
require(gbm)
require(nnet)
require(dismo)

# Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("ET_data", showWarnings=F)
dat_dir <- "./ET_data"

# download Geo-Wiki data
download("https://www.dropbox.com/s/qkgluhy31bhhsl8/ET_geow_31214.csv?dl=0", "./ET_data/ET_geow_31214.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/ET_geow_31214.csv", sep=""), header=T, sep=",")
geos <- na.omit(geos)

# download Ethiopia Gtifs (~35.5 Mb) and stack in raster
download("https://www.dropbox.com/s/xgwxukuj2q9dgbf/ET_grids.zip?dl=0", "./ET_Data/ET_grids.zip", mode="wb")
unzip("./ET_data/ET_grids.zip", exdir="./ET_data", overwrite=T)
glist <- list.files(path="./ET_data", pattern="tif", full.names=T)
grid <- stack(glist)

# Data setup --------------------------------------------------------------
# Project Geo-Wiki coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grid)

# Extract gridded variables at Geo-Wiki locations
geosgrid <- extract(grid, geos)

# Assemble dataframes
# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP <- geos$CRP
crpdat <- cbind.data.frame(CRP, geosgrid)
crpdat <- na.omit(crpdat)

# presence/absence of Human Settlements (HSP, present = Y, absent = N)
# note that this excludes large urban areas where MODIS fPAR = 0
HSP <- geos$HSP
hspdat <- cbind.data.frame(HSP, geosgrid)
hspdat <- na.omit(hspdat)

# Split data into train and test sets ------------------------------------
# set train/test set randomization seed
seed <- 1385321
set.seed(seed)

# Cropland train/test split
crpIndex <- createDataPartition(crpdat$CRP, p = 0.75, list = FALSE, times = 1)
crpTrain <- crpdat[ crpIndex,]
crpTest  <- crpdat[-crpIndex,]

# Human Settlement train/test split
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

# presence/absence of Human Settlements (HSP, present = Y, absent = N)
HSP.glm <- train(HSP ~ ., data = hspTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
hspglm.test <- predict(HSP.glm, hspTest) ## predict test-set
confusionMatrix(hspglm.test, hspTest$HSP, "Y") ## print validation summaries
hspglm.pred <- predict(grid, HSP.glm, type = "prob") ## spatial predictions

# <MASS> predictions
glmpreds <- stack(1-crpglm.pred, 1-hspglm.pred)
names(glmpreds) <- c("CRPglm", "HSPglm")

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

# presence/absence of Human Settlements (HSP, present = Y, absent = N)
HSP.rf <- train(HSP ~ ., data = hspTrain,
                method = "rf",
                trControl = oob)
hsprf.test <- predict(HSP.rf, hspTest) ## predict test-set
confusionMatrix(hsprf.test, hspTest$HSP, "Y") ## print validation summaries
hsprf.pred <- predict(grid, HSP.rf, type = "prob") ## spatial predictions

# <randomForest> predictions
rfpreds <- stack(1-crprf.pred, 1-hsprf.pred)
names(rfpreds) <- c("CRPrf", "HSPrf")

# Gradient boosting <gbm> -----------------------------------------------
# CV for training gbm's
gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.gbm <- train(CRP ~ ., data = crpTrain,
                 method = "gbm",
                 trControl = gbm)
crpgbm.test <- predict(CRP.gbm, crpTest) ## predict test-set
confusionMatrix(crpgbm.test, crpTest$CRP, "Y") ## print validation summaries
crpgbm.pred <- predict(grid, CRP.gbm, type = "prob") ## spatial predictions

# presence/absence of Human Settlements (HSP, present = Y, absent = N)
HSP.gbm <- train(HSP ~ ., data = hspTrain,
                 method = "gbm",
                 trControl = gbm)
hspgbm.test <- predict(HSP.gbm, hspTest) ## predict test-set
confusionMatrix(hspgbm.test, hspTest$HSP, "Y") ## print validation summaries
hspgbm.pred <- predict(grid, HSP.gbm, type = "prob") ## spatial predictions

# <gbm> predictions
gbmpreds <- stack(1-crpgbm.pred, 1-hspgbm.pred)
names(gbmpreds) <- c("CRPgbm", "HSPgbm")

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

# presence/absence of Human Settlements (HSP, present = Y, absent = N)
HSP.nn <- train(HSP ~ ., data = hspTrain,
                method = "nnet",
                trControl = nn)
hspnn.test <- predict(HSP.nn, hspTest) ## predict test-set
confusionMatrix(hspnn.test, hspTest$HSP, "Y") ## print validation summaries
hspnn.pred <- predict(grid, HSP.nn, type = "prob") ## spatial predictions

# <nnet> predictions
nnpreds <- stack(1-crpnn.pred, 1-hspnn.pred)
names(nnpreds) <- c("CRPnn", "HSPnn")

# Plot predictions by Geo-Wiki variables -----------------------------------
# Cropland prediction plots
crp.preds <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred)
names(crp.preds) <- c("glm","randomForest","gbm","nnet")
plot(crp.preds, axes = F)

# Human settlement prediction plots
hsp.preds <- stack(1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(hsp.preds) <- c("glm","randomForest","gbm","nnet")
plot(hsp.preds, axes = F)

# Ensemble predictions <glm>, <rf>, <gbm>, <nnet> --------------------------
# Ensemble set-up
pred <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred, 1-crpnn.pred,
              1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(pred) <- c("CRPglm","CRPrf","CRPgbm","CRPnn",
                 "HSPglm","HSPrf","HSPgbm","HSPnn")
geospred <- extract(pred, geos)

# presence/absence of Cropland (CRP, present = Y, absent = N)
crpens <- cbind.data.frame(CRP, geospred)
crpens <- na.omit(crpens)
crpensTest <- crpens[-crpIndex,] ## replicate previous test set

# presence/absence of Human Settlements (HSP, present = Y, absent = N)
hspens <- cbind.data.frame(HSP, geospred)
hspens <- na.omit(hspens)
hspensTest <- hspens[-hspIndex,] ## replicate previous test set

# GLM-based ensemble weighting on the test set
# 10-fold CV
ens <- trainControl(method = "cv", number = 10)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.ens <- train(CRP ~ CRPglm + CRPrf + CRPgbm + CRPnn, data = crpensTest,
                 family = binomial, 
                 method = "glm",
                 trControl = ens)
summary(CRP.ens)
crp.pred <- predict(CRP.ens, crpensTest, type="prob")
crp.test <- cbind(crpensTest, crp.pred)
crp <- subset(crp.test, CRP=="Y", select=c(Y))
cra <- subset(crp.test, CRP=="N", select=c(Y))
crp.eval <- evaluate(p=crp[,1], a=cra[,1])
crp.eval
threshold(crp.eval)
plot(crp.eval, 'ROC')
crpens.pred <- predict(pred, CRP.ens, type="prob")
plot(1-crpens.pred, axes = F)

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
HSP.ens <- train(HSP ~ HSPglm + HSPrf + HSPgbm + HSPnn, data = hspensTest,
                 family = binomial, 
                 method = "glm",
                 trControl = ens)
summary(HSP.ens)
hsp.pred <- predict(HSP.ens, hspensTest, type="prob")
hsp.test <- cbind(hspensTest, hsp.pred)
hsp <- subset(hsp.test, HSP=="Y", select=c(Y))
hsa <- subset(hsp.test, HSP=="N", select=c(Y))
hsp.eval <- evaluate(p=hsp[,1], a=hsa[,1])
hsp.eval
threshold(hsp.eval)
plot(hsp.eval, 'ROC')
hspens.pred <- predict(pred, HSP.ens, type="prob")
plot(1-hspens.pred, axes = F)

# Write spatial predictions -----------------------------------------------
# Create a "Results" folder in current working directory
dir.create("ET_results", showWarnings=F)

# Export Gtif's to "./ET_results"
# Individual model predictions
writeRaster(crp.preds, filename="./ET_results/ET_crpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(hsp.preds, filename="./ET_results/ET_hspreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
# Ensemble predictions
enspred <- stack(crpens.pred, hspens.pred)
writeRaster(enspred, filename="./ET_results/ET_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)






