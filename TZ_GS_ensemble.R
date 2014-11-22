# Ensemble predictions of Tanzania GeoSurvey cropland, woody vegetation cover
# and human settlement observations. 
# M. Walsh, November 2014

# Required packages
# install.packages(c("downloader","raster","rgdal","caret","MASS","rpart","randomForest")), dependencies=TRUE)
require(downloader)
require(raster)
require(caret)
require(rgdal)
require(MASS)
require(rpart)
require(randomForest)

# Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("Data", showWarnings=F)
dat_dir <- "./Data"

# download GeoSurvey data
download("https://www.dropbox.com/s/03l4m4zjdi5mhyu/TZ_geos_1114.csv?dl=0", "./Data/TZ_geos_1114.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/TZ_geos_1114.csv", sep=""), header=T, sep=",")
geos <- na.omit(geos)

# download Tanzania Gtifs (~27.2 Mb) and stack in raster
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
HSP <- geos$HSP
hspdat <- cbind.data.frame(HSP, geosgrid)
hspdat <- na.omit(hspdat)

# Split data into train and test sets ------------------------------------
set.seed(1385321)

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

# Cross-validation setup --------------------------------------------------
# 10-fold
cv10 <- trainControl(method = "cv", number = 10)

# Classification models ---------------------------------------------------
# Stepwise main effects GLM's
CRP.glm <- train(CRP ~ ., data = crpTrain,
                 family = binomial, 
                 method = "glmStepAIC",
                 trControl = cv10)
crpglm.test <- predict(CRP.glm, crpTest, type = "prob")
crpglm.pred <- predict(grid, CRP.glm, type = "prob")
plot(1-crpglm.pred)

WCP.glm <- train(WCP ~ ., data = wcpTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = cv10)
wcpglm.test <- predict(WCP.glm, wcpTest, type = "prob")
wcpglm.pred <- predict(grid, WCP.glm, type = "prob")
plot(1-wcpglm.pred)

HSP.glm <- train(HSP ~ ., data = hspTrain,
                 family=binomial, 
                 method = "glmStepAIC",
                 trControl = cv10)
hspglm.test <- predict(HSP.glm, hspTest, type = "prob")
hspglm.pred <- predict(grid, HSP.glm, type = "prob")
plot(1-hspglm.pred)

