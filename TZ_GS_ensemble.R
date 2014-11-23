# Ensemble machine learning predictions of Tanzania GeoSurvey cropland,
# woody vegetation cover and human settlement observations. 
# M. Walsh, November 2014

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","MASS","rpart","randomForest")), dependencies=TRUE)
require(downloader)
require(raster)
require(caret)
require(rgdal)
require(MASS)
require(randomForest)
require(gbm)

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
# presence/absence of Cropland (CRP, present = P, absent = A)
CRP <- geos$CRP
crpdat <- cbind.data.frame(CRP, geosgrid)
crpdat <- na.omit(crpdat)

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = P, absent = A)
WCP <- geos$WCP
wcpdat <- cbind.data.frame(WCP, geosgrid)
wcpdat <- na.omit(wcpdat)

# presence/absence of Buildings/Human Settlements (HSP, present = P, absent = A)
# note that this excludes urban areas where MODIS fPAR = 0
HSP <- geos$HSP
hspdat <- cbind.data.frame(HSP, geosgrid)
hspdat <- na.omit(hspdat)

# Split data into train and test sets ------------------------------------
set.seed(1385321) ## very important to set this as it will affect all subsequent calc's !!!

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

# presence/absence of Cropland (CRP, present = P, absent = A)
CRP.glm <- train(CRP ~ ., data = crpTrain,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = step)
crpglm.test <- predict(CRP.glm, crpTest) ## predict test-set
confusionMatrix(crpglm.test, crpTest$CRP) ## print validation summaries
crpglm.pred <- predict(grid, CRP.glm, type = "prob") ## spatial predictions
plot(1-crpglm.pred) ## map presence

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = P, absent = A)
WCP.glm <- train(WCP ~ ., data = wcpTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
wcpglm.test <- predict(WCP.glm, wcpTest) ## predict test-set
confusionMatrix(wcpglm.test, wcpTest$WCP) ## print validation summaries
wcpglm.pred <- predict(grid, WCP.glm, type = "prob") ## spatial predictions
plot(1-wcpglm.pred) ## map presence

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = P, absent = A)
HSP.glm <- train(HSP ~ ., data = hspTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = step)
hspglm.test <- predict(HSP.glm, hspTest) ## predict test-set
confusionMatrix(hspglm.test, hspTest$HSP) ## print validation summaries
hspglm.pred <- predict(grid, HSP.glm, type = "prob") ## spatial predictions
plot(1-hspglm.pred) ## map presence

# Random forests <randomForest> -------------------------------------------
# out-of-bag CV
oob <- trainControl(method = "oob")

# presence/absence of Cropland (CRP, present = P, absent = A)
CRP.rf <- train(CRP ~ ., data = crpTrain,
                method = "rf",
                trControl = oob)
crprf.test <- predict(CRP.rf, crpTest) ## predict test-set
confusionMatrix(crprf.test, crpTest$CRP) ## print validation summaries
crprf.pred <- predict(grid, CRP.rf, type = "prob") ## spatial predictions
plot(1-crprf.pred) ## map presence

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = P, absent = A)
WCP.rf <- train(WCP ~ ., data = wcpTrain,
                method = "rf",
                trControl = oob)
wcprf.test <- predict(WCP.rf, wcpTest) ## predict test-set
confusionMatrix(wcprf.test, wcpTest$WCP) ## print validation summaries
wcprf.pred <- predict(grid, WCP.rf, type = "prob") ## spatial predictions
plot(1-wcprf.pred) ## map presence

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = P, absent = A)
HSP.rf <- train(HSP ~ ., data = hspTrain,
                method = "rf",
                trControl = oob)
hsprf.test <- predict(HSP.rf, hspTest) ## predict test-set
confusionMatrix(hsprf.test, hspTest$HSP) ## print validation summaries
hsprf.pred <- predict(grid, HSP.rf, type = "prob") ## spatial predictions
plot(1-hsprf.pred) ## map presence

# Generalized boosted regression <gbm> ------------------------------------
# CV for training gbm's
gbm <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

# presence/absence of Cropland (CRP, present = P, absent = A)
CRP.gbm <- train(CRP ~ ., data = crpTrain,
                 method = "gbm",
                 trControl = gbm)
crpgbm.test <- predict(CRP.gbm, crpTest) ## predict test-set
confusionMatrix(crpgbm.test, crpTest$CRP) ## print validation summaries
crpgbm.pred <- predict(grid, CRP.gbm, type = "prob") ## spatial predictions
plot(1-crpgbm.pred) ## map presence

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = P, absent = A)
WCP.gbm <- train(WCP ~ ., data = wcpTrain,
                 method = "gbm",
                 trControl = gbm)
wcpgbm.test <- predict(WCP.gbm, wcpTest) ## predict test-set
confusionMatrix(wcpgbm.test, wcpTest$WCP) ## print validation summaries
wcpgbm.pred <- predict(grid, WCP.gbm, type = "prob") ## spatial predictions
plot(1-wcpgbm.pred) ## map presence

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = P, absent = A)
HSP.gbm <- train(HSP ~ ., data = hspTrain,
                 method = "gbm",
                 trControl = gbm)
hspgbm.test <- predict(HSP.gbm, hspTest) ## predict test-set
confusionMatrix(hspgbm.test, hspTest$HSP) ## print validation summaries
hspgbm.pred <- predict(grid, HSP.gbm, type = "prob") ## spatial predictions
plot(1-hspgbm.pred) ## map presence

# Ensemble predictions <glm>, <rf>, <gbm> ---------------------------------
# Ensemble set up
preds <- stack(1-crpglm.pred, 1-crprf.pred, 1-crpgbm.pred,
               1-wcpglm.pred, 1-wcprf.pred, 1-wcpgbm.pred,
               1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred)
names(preds) <- c("CRPglm","CRPrf","CRPgbm","WCPglm","WCPrf","WCPgbm","HSPglm","HSPrf","HSPgbm")
plot(preds)
expreds <- extract(preds, geos)
ensdat <- data.frame(geos[,4:6], expreds)
ensIndex <- createDataPartition(ensdat$HSP, p = 0.75, list = FALSE, times = 1)
ensTest  <- ensdat[-ensIndex,] ## previous set.seed important here to replicate test set

# GLM based weighting on test set & spatial predictions
# presence/absence of Cropland (CRP, present = P, absent = A)
CRP.ens <- glm(CRP~CRPglm+CRPrf+CRPgbm, family=binomial(link="logit"), data=ensTest)
crpens.pred <- predict(preds, CRP.ens, type="response")
plot(crpens.pred)

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = P, absent = A)
WCP.ens <- glm(WCP~WCPglm+WCPrf+WCPgbm, family=binomial(link="logit"), data=ensTest)
wcpens.pred <- predict(preds, WCP.ens, type="response")
plot(wcpens.pred)

# presence/absence of (rural) Buildings/Human Settlements (HSP, present = P, absent = A)
HSP.ens <- glm(HSP~HSPglm+HSPrf+HSPgbm, family=binomial(link="logit"), data=ensTest)
hspens.pred <- predict(preds, HSP.ens, type="response")
plot(hspens.pred)

# Write spatial predictions -----------------------------------------------
# Create a "Results" folder in your current working directory
dir.create("Results", showWarnings=F)

# Export Gtif's to "./Results"
writeRaster(preds, filename="./Results/TZ_gspreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
enspred <- stack(crpens.pred, wcpens.pred, hspens.pred)
names(enspred) <- c("CRPens", "WCPens", "HSPens")
writeRaster(enspred, filename="./Results/TZ_enspred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
