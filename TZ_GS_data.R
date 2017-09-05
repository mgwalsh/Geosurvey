# Tanzania GeoSurvey 250m resolution data setup 
# M. Walsh, September 2017

# required packages
# install.packages(c("downloader","rgdal","raster","caret")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(caret)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("TZ_GS250", showWarnings=F)
setwd("./TZ_GS250")

# download GeoSurvey data
download("https://www.dropbox.com/s/5k4awl2se6s982y/TZ_geos_082317.csv.zip?raw=1", "TZ_geos_082317.csv.zip", mode="wb")
unzip("TZ_geos_082317.csv.zip", overwrite=T)
geos <- read.table("TZ_geos_082317.csv", header=T, sep=",")

# download Tanzania Gtifs and stack in raster (note this is a big 500+ Mb download)
download("https://www.dropbox.com/s/16ifi3ucqxcaf5j/TZ_GRIDS250m.zip?raw=1", "TZ_GRIDS250m.zip", mode="wb")
unzip("TZ_GRIDS250m.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Data setup ---------------------------------------------------------------
# project GeoSurvey coords to grid CRS
geos.proj <- as.data.frame(project(cbind(geos$lon, geos$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.proj) <- c("x","y")
geos <- cbind(geos, geos.proj)
coordinates(geos) <- ~x+y
projection(geos) <- projection(grids)

# extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)

# Extract gridded variables at GeoSurvey locations
geosgrid <- extract(grids, geos)

# Assemble dataframes
# presence/absence of Buildings (BP, present = Y, absent = N)
BP <- geos$BP
bpdat <- cbind.data.frame(BP, geosgrid)
bpdat <- na.omit(bpdat)

# presence/absence of Cropland (CP, present = Y, absent = N)
CP <- geos$CP
cpdat <- cbind.data.frame(CP, geosgrid)
cpdat <- na.omit(cpdat)

# presence/absence of Woody Vegetation Cover >60% (WP, present = Y, absent = N)
WP <- geos$WP
wpdat <- cbind.data.frame(WP, geosgrid)
wpdat <- na.omit(wpdat)

# Split data into train and test sets -------------------------------------
# set train/test set randomization seed
seed <- 1385321
set.seed(seed)

# Buildings train/test split
bpIndex <- createDataPartition(bpdat$BP, p = 2/3, list = FALSE, times = 1)
bpTrain <- bpdat[ bpIndex,]
bpTest  <- bpdat[-bpIndex,]

# Cropland train/test split
cpIndex <- createDataPartition(cpdat$CP, p = 2/3, list = FALSE, times = 1)
cpTrain <- cpdat[ cpIndex,]
cpTest  <- cpdat[-cpIndex,]

# Woody cover train/test split
wpIndex <- createDataPartition(wpdat$WP, p = 2/3, list = FALSE, times = 1)
wpTrain <- wpdat[ wpIndex,]
wpTest  <- wpdat[-wpIndex,]
