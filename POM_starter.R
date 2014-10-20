# AfSIS sentinel site environmental envelopes / presence-only models (POM)
# M. Walsh, October 2014

# Set local working directory e.g.
dat_dir <- "/Users/markuswalsh/Documents/LDSF/Af_grids_1k"
setwd(dat_dir)

# Required packages
# install.packages(c("downloader","proj4","raster")), dependencies=TRUE)
require(downloader)
require(proj4)
require(raster)

# Data downloads ----------------------------------------------------------

# AfSIS sentinel site GPS locations
download("https://www.dropbox.com/s/ptw7sv2s7oxpm9x/AfSIS_GPS.csv?dl=0", "AfSIS_GPS.csv", mode="wb")
geos <- read.table("AfSIS_GPS.csv", header=T, sep=",")

# Africa PC grid download (~295 Mb)
download("https://www.dropbox.com/s/jhzqx0a4f90owfq/Af_PC_1k.zip?dl=0", "Af_PC_1k.zip", mode="wb")
unzip("Af_PC_1k.zip", overwrite=T)

grid.list <- c("PC1.tif","PC2.tif","PC3.tif","PC4.tif")
grids <- stack(grid.list)
plot(grids)
  
# Overlay points and grids ------------------------------------------------

# Project sentinel site locations to Africa LAEA from LonLat
geos.laea <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.laea) <- c("x","y")
geos <- cbind(geos, geos.laea)
coordinates(geos) = ~x+y
proj4string(geos) = CRS("+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs")

# Overlay
test <- extract(grids, geos)




