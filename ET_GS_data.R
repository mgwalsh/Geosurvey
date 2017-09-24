# Ethiopia GeoSurvey 250m resolution data setup 
# M. Walsh, September 2017

# Required packages
# install.packages(c("downloader","rgdal","raster","caret")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
  require(caret)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("ET_GS250", showWarnings=F)
setwd("./ET_GS250")

# download GeoSurvey data
download("https://www.dropbox.com/s/lcn8ntiaa1vdn1t/ET_geos_200917.csv.zip?raw=1", "ET_geos_200917.csv.zip", mode="wb")
unzip("ET_geos_200917.csv.zip", overwrite=T)
geos <- read.table("ET_geos_200917.csv", header=T, sep=",")

# download Ethiopia Gtifs and stack in raster (note this is a big 900+ Mb download)
download("https://www.dropbox.com/s/iqix6sn66w04jo0/ET_250m_2017.zip?raw=1", "ET_250m_2017.zip", mode="wb")
unzip("ET_250m_2017.zip", overwrite=T)
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

# Assemble dataframes
# presence/absence of Buildings (BP, present = Y, absent = N)
BP <- geos$BP
bpdat <- as.data.frame(cbind(BP, geosgrid))
bpdat <- na.omit(bpdat)
write.csv(bpdat, "BP.csv", row.names = FALSE)

# presence/absence of Cropland (CP, present = Y, absent = N)
CP <- geos$CP
cpdat <- as.data.frame(cbind(CP, geosgrid))
cpdat <- na.omit(cpdat)
write.csv(cpdat, "CP.csv", row.names = FALSE)

# presence/absence of Woody Vegetation Cover >60% (WP, present = Y, absent = N)
WP <- geos$WP
wpdat <- as.data.frame(cbind(WP, geosgrid))
wpdat <- na.omit(wpdat)
write.csv(wpdat, "WP.csv", row.names = FALSE)

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
