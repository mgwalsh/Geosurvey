# Tanzania geosurvey data setup
# M. Walsh, K. Peng, J. Mulvey, A. Verlinden (Sep. 2014)

# Set working directory e.g.
dat_dir <- "/Users/markuswalsh/Documents/Tanzania/TZ_Geosurvey"
setwd(dat_dir)

# Required packages
# install.packages(c("downloader","proj4","raster")), dependencies=TRUE)
require(downloader)
require(proj4)
require(raster)

# Data downloads ----------------------------------------------------------

# Geosurvey data
download("https://www.dropbox.com/s/8dgzphvwn0ls4sq/TZ_geos_0914.csv?dl=0", "TZ_geos_0914.csv", mode="wb")
geos <- read.table("TZ_geos_0914.csv", header=T, sep=",")

# Tanzania Gtifs (~35.5 Mb)
download("https://www.dropbox.com/s/0ucx0oqx8c4lpej/TZ_grids.zip?dl=0", "TZ_grids.zip", mode="wb")
unzip("TZ_grids.zip", overwrite=T)

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

grid.list <- c("BSANs.tif","BSASs.tif","BSAVs.tif","CTIs.tif","ELEVs.tif","EVIAs.tif","LSTDs.tif","LSTNs.tif","REF1s.tif","REF2s.tif","REF3s.tif","REF7s.tif","RELIs.tif","TMAPs.tif","TMFIs.tif")

for (i in 1:length(grid.list)){
  print(paste("extracting", grid.list[i]))
  grid.cov <- raster(grid.list[i]) 
  geos.gid@data[strsplit(grid.list[i], split=".tif")[[1]]] <- extract(
    x = grid.cov, 
    y = geos.gid,
    method = "simple")
}
TZ_GS_data <- as.data.frame(geos.gid)

# Write csv ---------------------------------------------------------------

write.csv(TZ_GS_data, "TZ_geosurvey.csv")






