#' Ensemble predictions of Marsabit GeoSurvey settlement observations. 
#' M. Walsh, August 2016

# Required packages
# install.packages(c("downloader","raster","rgdal","plyr","caret","randomForest","gbm","nnet","glmnet","ROSE","dismo")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(caret)
require(doParallel)
require(randomForest)
require(gbm)
require(nnet)
require(plyr)
require(ipred)
require(glmnet)
require(ROSE)
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

# Split data into train and test sets -------------------------------------
# set train/test set randomization seed
seed <- 1385321
set.seed(seed)

# Settlement train/test split
rspIndex <- createDataPartition(rspdat$RSP, p = 3/4, list = FALSE, times = 1)
rspTrain <- rspdat[ rspIndex,]
rspTest  <- rspdat[-rspIndex,]
prop.table(table(rspTrain$RSP))

# balance the data with SMOTE <ROSE> for model training
rspROSE <- ROSE(RSP ~ ., data = rspTrain, seed = seed)$data
table(rspROSE$RSP)

# Random forests <randomForest> -------------------------------------------
# Start foreach to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

tc <- trainControl(method = "oob", 
                   classProbs = TRUE, 
                   summaryFunction = twoClassSummary,
                   allowParallel = TRUE)
tg <- expand.grid(mtry=seq(2, 20, by=2))

# imbalanced training data
set.seed(seed)
RSP.rf <- train(RSP ~ ., data = rspTrain, 
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc,
                metric = "Kappa")
print(RSP.rf)
RSP.imp <- varImp(RSP.rf, useModel = FALSE)
plot(RSP.imp, top=23)
RSP_rf <- predict(grids, RSP.rf, type="prob")
plot(1-RSP_rf, axes = F)

# balanced data <ROSE>
set.seed(seed)
RSP.rfB <- train(RSP ~ ., data = rspROSE, 
                preProc = c("center", "scale"),
                method = "rf",
                ntree = 501,
                tuneGrid = tg,
                trControl = tc,
                metric = "Kappa")
print(RSP.rfB)
RSP.imp <- varImp(RSP.rfB, useModel = FALSE)
plot(RSP.imp, top=23)
RSP_rfB <- predict(grids, RSP.rfB, type="prob")
plot(1-RSP_rfB, axes = F)

# Gradient boosting <gbm> -------------------------------------------------
tc <- trainControl(method = "repeatedcv", repeats=10,
                   classProbs = TRUE, 
                   summaryFunction = twoClassSummary,
                   allowParallel = TRUE)

# imbalanced training data
set.seed(seed)
RSP.gbm <- train(RSP ~ ., data = rspTrain, 
                 method = "gbm", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 metric = "ROC",
                 tuneGrid = expand.grid(.n.trees = seq(50,500,by=50), 
                                        .interaction.depth = 3,
                                        .shrinkage = 0.1,
                                        .n.minobsinnode = 100))
print(RSP.gbm)
RSP.imp <- varImp(RSP.gbm)
plot(RSP.imp, top=23)
RSP_gbm <- predict(grids, RSP.gbm, type="prob")
plot(1-RSP_gbm, axes = F)

# balanced data <ROSE>
set.seed(seed)
RSP.gbB <- train(RSP ~ ., data = rspROSE, 
                 method = "gbm", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 metric = "ROC",
                 tuneGrid = expand.grid(.n.trees = seq(50,500,by=50), 
                                        .interaction.depth = 3,
                                        .shrinkage = 0.1,
                                        .n.minobsinnode = 100))
print(RSP.gbB)
RSP.imp <- varImp(RSP.gbB)
plot(RSP.imp, top=23)
RSP_gbB <- predict(grids, RSP.gbB, type="prob")
plot(1-RSP_gbB, axes = F)

# Neural network <nnet> ---------------------------------------------------
tc <- trainControl(method = "repeatedcv", repeats=10,
                   classProbs = TRUE, 
                   summaryFunction = twoClassSummary,
                   allowParallel = TRUE)

# imbalanced training data
set.seed(seed)
RSP.nn <- train(RSP ~ ., data = rspTrain, 
                method = "nnet", 
                preProc = c("center", "scale"), 
                trControl = tc,
                metric = "ROC")
print(RSP.nn)
RSP.imp <- varImp(RSP.nn, useModel = FALSE)
plot(RSP.imp, top=23)
RSP_nn <- predict(grids, RSP.nn, type="prob")
plot(1-RSP_nn, axes = F)

# balanced data <ROSE>
set.seed(seed)
RSP.nnB <- train(RSP ~ ., data = rspROSE, 
                method = "nnet", 
                preProc = c("center", "scale"), 
                trControl = tc,
                metric = "ROC")
print(RSP.nnB)
RSP.imp <- varImp(RSP.nnB, useModel = FALSE)
plot(RSP.imp, top=23)
RSP_nnB <- predict(grids, RSP.nnB, type="prob")
plot(1-RSP_nnB, axes = F)

# Bagged CART <ipred> -----------------------------------------------------
tc <- trainControl(method = "repeatedcv", repeats=10,
                   classProbs = TRUE, 
                   summaryFunction = twoClassSummary,
                   allowParallel = TRUE)

# imbalanced training data
set.seed(seed)
RSP.tb <- train(RSP ~ ., data = rspTrain, 
                method = "treebag", 
                preProc = c("center", "scale"),
                trControl = tc,
                nbagg = 50,
                metric = "ROC")
print(RSP.tb)
RSP.imp <- varImp(RSP.tb, useModel = FALSE)
plot(RSP.imp, top=23)
RSP_tb <- predict(grids, RSP.tb, type="prob")
plot(1-RSP_tb, axes = F)

# balanced data <ROSE>
set.seed(seed)
RSP.tbB <- train(RSP ~ ., data = rspROSE, 
                method = "treebag", 
                preProc = c("center", "scale"),
                trControl = tc,
                nbagg = 50,
                metric = "ROC")
print(RSP.tbB)
RSP.imp <- varImp(RSP.tbB, useModel = FALSE)
plot(RSP.imp, top=23)
RSP_tbB <- predict(grids, RSP.tbB, type="prob")
plot(1-RSP_tbB, axes = F)

stopCluster(mc)


