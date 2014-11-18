# Ensemble predictions of Tanzania GeoSurvey cropland, woody vegetation cover
# and human settlement data. 
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
download("https://www.dropbox.com/s/8dgzphvwn0ls4sq/TZ_geos_0914.csv?dl=0", "./Data/TZ_geos_0914.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/TZ_geos_0914.csv", sep=""), header=T, sep=",")

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

# Extract gridded variables at GeoSurvey locations with <raster>
geosgrid <- extract(grid, geos)

# Assemble dataframes
# presence/absence of Cropland (CRP, present = 1, absent = 0)
CRP <- geos$CRP
crpdat <- data.frame(cbind(CRP, geosgrid))
crpdat <- na.omit(crpdat)

# presence/absence of Woody Vegetation Cover of >60% (WCP, present = 1, absent = 0)
WCP <- geos$WCP
wcpdat <- data.frame(cbind(WCP, geosgrid))
wcpdat <- na.omit(wcpdat)

# presence/absence of buildings/Human Settlements (HSP, present = 1, absent = 0)
HSP <- geos$HSP
hspdat <- data.frame(cbind(HSP, geosgrid))
hspdat <- na.omit(hspdat)

# Split data into train and test sets (with <caret>) ----------------------

