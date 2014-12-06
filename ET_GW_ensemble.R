# Ensemble machine learning predictions of Ethiopia Geo-Wiki cropland,
# and human settlement observations.
# Cropland & Human Settlement data courtesy http://www.geo-wiki.org/download-data
# M. Walsh, December 2014

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
dir.create("ET_data", showWarnings=F)
dat_dir <- "./ET_data"

# download Geo-Wiki data
download("https://www.dropbox.com/s/qkgluhy31bhhsl8/ET_geow_31214.csv?dl=0", "./ET_data/ET_geow_31214.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/ET_geow_31214.csv", sep=""), header=T, sep=",")
geos <- na.omit(geos)

# download Ethiopia Gtifs (~34.6 Mb) and stack in raster
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

