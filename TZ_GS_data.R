# Tanzania GeoSurvey 250m resolution data setup 
# M. Walsh, September 2017

# Required packages
# install.packages(c("downloader","rgdal","raster")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
})

# Data downloads -----------------------------------------------------------
# set working directory
dir.create("TZ_GS250", showWarnings=F)
setwd("./TZ_GS250")

# download GeoSurvey data
download("https://www.dropbox.com/s/57kuxbkm5sv092a/TZ_geos_2017.csv.zip?raw=1", "TZ_geos_2017.csv.zip", mode="wb")
unzip("TZ_geos_2017.csv.zip", overwrite=T)
geos <- read.table("TZ_geos_2017.csv", header=T, sep=",")

# download Tanzania Gtifs and stack in raster (note this is a big 550+ Mb download)
download("https://www.dropbox.com/s/pshrtvjf7navegu/TZ_250m_2017.zip?raw=1", "TZ_250m_2017.zip", mode="wb")
unzip("TZ_250m_2017.zip", overwrite=T)
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
gsdat <- as.data.frame(cbind(geos, geosgrid)) 
gsdat <- na.omit(gsdat)

# Write output file -------------------------------------------------------
write.csv(gsdat, "gsdat.csv", row.names = FALSE)
