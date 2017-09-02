#' Ensemble predictions of Marsabit GeoSurvey settlement, bare area & woody vegetation cover observations. 
#' M. Walsh, August 2016

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","randomForest","gbm","deepnet","glmnet","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(doParallel)
require(randomForest)
require(gbm)
require(deepnet)
require(glmnet)
require(dismo)

# Data downloads -----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("MS_data", showWarnings=F)
setwd("./MS_data")

# download GeoSurvey data
download("https://www.dropbox.com/s/yupshubof01r8fy/Marsabit_GS.csv.zip?dl=0", "Marsabit_GS.csv.zip", mode="wb")
unzip("Marsabit_GS.csv.zip", overwrite=T)
geos <- read.table("Marsabit_GS.csv", header=T, sep=",")

# Download grids
download("https://www.dropbox.com/s/awtcarx0l89ft3y/Marsabit_grids.zip?dl=0", "Marsabit_grids.zip", mode="wb")
unzip("Marsabit_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Data setup ---------------------------------------------------------------
# Project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# Extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)

# Assemble dataframes
# presence/absence of Buildings/Settlements (RSP, present = Y, absent = N)
RSP <- geos$RSP
rspdat <- cbind.data.frame(RSP, geosgrid)
rspdat <- na.omit(rspdat)

# presence/absence of >60% bare (uvegetated) area (BAP, present = Y, absent = N)
BAP <- geos$BAP
bapdat <- cbind.data.frame(BAP, geosgrid)
bapdat <- na.omit(bapdat)

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Y, absent = N)
WCP <- geos$WCP
wcpdat <- cbind.data.frame(WCP, geosgrid)
wcpdat <- na.omit(wcpdat)

# Split data into train and test sets -------------------------------------
# set train/test set randomization seed
seed <- 1385321
set.seed(seed)

# Settlement train/test split
rspIndex <- createDataPartition(rspdat$RSP, p = 2/3, list = FALSE, times = 1)
rspTrain <- rspdat[ rspIndex,]
rspTest  <- rspdat[-rspIndex,]

# Bare area train/test split
bapIndex <- createDataPartition(bapdat$BAP, p = 2/3, list = FALSE, times = 1)
bapTrain <- bapdat[ bapIndex,]
bapTest  <- bapdat[-bapIndex,]

# Woody cover train/test split
wcpIndex <- createDataPartition(wcpdat$WCP, p = 2/3, list = FALSE, times = 1)
wcpTrain <- wcpdat[ wcpIndex,]
wcpTest  <- wcpdat[-wcpIndex,]

# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Random forests <randomForest> -------------------------------------------
# presence/absence of Buildings/Settlements (RSP, present = Y, absent = N)
tc <- trainControl(method = "oob", allowParallel = TRUE)
tg <- expand.grid(mtry=seq(20, 200, by=10))
RSP.rf <- train(RSP ~ ., data = rspTrain, 
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(RSP.rf)
RSP.imp <- varImp(RSP.rf, useModel = FALSE)
plot(RSP.imp, top=23)
RSP_rf <- predict(grids, RSP.rf, type="prob")

# presence/absence of >40% (vegetated) area (VCP, present = Y, absent = N)
tc <- trainControl(method = "oob", allowParallel = TRUE)
tg <- expand.grid(mtry=seq(20, 200, by=10))
VCP.rf <- train(BAP ~ ., data = bapTrain, 
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(VCP.rf)
VCP.imp <- varImp(VCP.rf, useModel = FALSE)
plot(VCP.imp, top=23)
VCP_rf <- predict(grids, VCP.rf, type="prob")

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = Y, absent = N)
tc <- trainControl(method = "oob", allowParallel = TRUE)
tg <- expand.grid(mtry=seq(20, 200, by=10))
WCP.rf <- train(WCP ~ ., data = wcpTrain, 
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc)
print(WCP.rf)
WCP.imp <- varImp(WCP.rf, useModel = FALSE)
plot(WCP.imp, top=23)
WCP_rf <- predict(grids, WCP.rf, type="prob")

# Gradient boosting <gbm> -------------------------------------------------
# presence/absence of Buildings/Settlements (RSP, present = Y, absent = N)
tc <- trainControl(method = "cv", number=50, allowParallel = TRUE)
RSP.gbm <- train(RSP ~ ., data = rspTrain, 
                 method = "gbm", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                        .interaction.depth = 3,
                                        .shrinkage = 0.1,
                                        .n.minobsinnode = 100))
print(RSP.gbm)
RSP.imp <- varImp(RSP.gbm)
plot(RSP.imp, top=23)

# presence/absence of >40% (vegetated) area (VCP, present = Y, absent = N)
VCP.gbm <- train(BAP ~ ., data = bapTrain, 
                 method = "gbm", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 tuneGrid = expand.grid(.n.trees=seq(50,500,by=50), 
                                        .interaction.depth = 3,
                                        .shrinkage = 0.1,
                                        .n.minobsinnode = 100))
print(VCP.gbm)
VCP.imp <- varImp(VCP.gbm)
plot(VCP.imp, top=23)


RF.preds <- stack(1-RSP_rf, BAP_rf, 1-WCP_rf)
names(RF.preds) <- c("RSP","BAP","WCP")
plot(RF.preds, axes = F)

# Write spatial predictions -----------------------------------------------
# Create a "Results" folder in current working directory
dir.create("MS_results", showWarnings=F)

# Export Gtif's to "./MS_results"
writeRaster(RF.preds, filename="./MS_results/MS_RFpreds.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

