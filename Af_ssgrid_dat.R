# AfSIS sentinel site remote sensing data
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

# Africa PC grid download (~325 Mb)
download("https://www.dropbox.com/s/jhzqx0a4f90owfq/Af_PC_1k.zip?dl=0", "Af_PC_1k.zip", mode="wb")
unzip("Af_PC_1k.zip", overwrite=T)

# LAEA CRS & grid ID's ----------------------------------------------------

# Project to Africa LAEA from LonLat
geos.laea <- as.data.frame(project(cbind(geos$Lon, geos$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geos.laea) <- c("x","y")
geos <- cbind(geos, geos.laea)

# Generate AfSIS grid cell ID's (GID)
res.pixel <- 1000
xgid <- ceiling(abs(geos$x)/res.pixel)
ygid <- ceiling(abs(geos$y)/res.pixel)
gidx <- ifelse(geos$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(geos$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
geos.gid <- cbind(geos, GID)
coordinates(geos.gid) = ~x+y
proj4string(geos.gid) = CRS("+proj=laea +datum=WGS84 +ellps=WGS84 +lat_0=5 +lon_0=20 +no_defs")

# Grid overlay ------------------------------------------------------------

grid.list <- c("PC1.tif","PC2.tif","PC3.tif")

for (i in 1:length(grid.list)){
  print(paste("extracting", grid.list[i]))
  grid.cov <- raster(grid.list[i]) 
  geos.gid@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
    x = grid.cov, 
    y = geos.gid,
    method = "simple")
}
ssgrid <- as.data.frame(geos.gid)
