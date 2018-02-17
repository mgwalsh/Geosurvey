# GeoSurvey building count data setup 
# M. Walsh & J. Chen, February 2018

# Required packages
# install.packages(c("downloader","jsonlite")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(downloader)
  require(jsonlite)
  require(leaflet)
})

# Data downloads ----------------------------------------------------------
# set working directory
dir.create("test", showWarnings=F)
setwd("./test")

# download GeoSurvey data
download("https://www.dropbox.com/s/qyz01a6dbtae06p/export-malawi-buildings.csv.zip?raw=1", "export-malawi-buildings.csv.zip", mode="wb")
unzip("export-malawi-buildings.csv.zip", overwrite=T)
geos <-  read.csv("export-malawi-buildings.csv", stringsAsFactors = F)
is.na(geos$building_loc) <- geos$building == "No"
bp <- geos[ which(geos$building == "Yes"), ]

# Count number of buildings per quadrat -----------------------------------
n <- rep(NA, nrow(bp))
for(i in 1:nrow(bp)) {
  t <- fromJSON(bp$building_loc[i])
  n[i] <- nrow(t$features)
}
n ## vector of number of buildings per quadrat with buildings

# Write file --------------------------------------------------------------
bp <- cbind(bp, n)
write.csv(bp, "building_count.csv", row.names = F)

# GeoSurvey map widget ----------------------------------------------------
# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addCircleMarkers(geos$Longitude, geos$Latitude, clusterOptions = markerClusterOptions())
w ## plot widget 

# save widget
saveWidget(w, 'MW_GS.html', selfcontained = T)
