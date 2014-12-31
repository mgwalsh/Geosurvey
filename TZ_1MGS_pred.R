# Ensemble predictions of Tanzania 1M GeoSurvey cropland and human settlement observations 
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
dir.create("TZ_data", showWarnings=F)
dat_dir <- "./TZ_data"

# download GeoSurvey data
download("https://www.dropbox.com/s/eq19mgnj86d4qo2/1MGS_123114.csv?dl=0", "./TZ_data/1MGS_123114.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/1MGS_123114.csv", sep=""), header=T, sep=",")
geos <- na.omit(geos)

# download Tanzania Gtifs (~27.9 Mb) and stack in raster
download("https://www.dropbox.com/s/otiqe78s0kf1z1s/TZ_grids.zip?dl=0", "./TZ_data/TZ_grids.zip", mode="wb")
unzip("./TZ_data/TZ_grids.zip", exdir="./TZ_data", overwrite=T)
glist <- list.files(path="./TZ_data", pattern="tif", full.names=T)
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

# presence/absence of Buildings/Human Settlements (HSP, present = Y, absent = N)
# note that this excludes large urban areas where MODIS fPAR = 0
HSP <- geos$HSP
hspdat <- cbind.data.frame(HSP, geosgrid)
hspdat <- na.omit(hspdat)

# set randomization seed
seed <- 1385321
set.seed(seed)

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

# Human settlement prediction plots
hsp.preds <- stack(1-hspglm.pred, 1-hsprf.pred, 1-hspgbm.pred, 1-hspnn.pred)
names(hsp.preds) <- c("glmStepAIC","randomForest","gbm","nnet")
plot(hsp.preds, axes = F)

