# Incoming GeoSurvey data check's
# M. Walsh, W. Wu, J. Chen, Dec. 2014

# Required packages
# install.packages(c("downloader","caret")), dependencies = T)
require(downloader)
require(caret)

# Load current GeoSurvey data ---------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("1MGS_data", showWarnings=F)
dat_dir <- "./1MGS_data"

# download GeoSurvey data
download("https://www.dropbox.com/s/a9tvwodid2g3ijd/1MGS_18_22_Dec_14.csv?dl=0", "./1MGS_data/1MGS_18_22_Dec_14.csv", mode="wb")
geos <- read.table(paste(dat_dir, "/1MGS_18_22_Dec_14.csv", sep=""), header=T, sep=",")

# Generate check sample ---------------------------------------------------
# set check set randomization seed
seed <- 1385321
set.seed(seed)

# Check split
checkIndex <- createDataPartition(geos$CP, p = 0.99, list = FALSE, times = 1)
checkTest  <- geos[-checkIndex, ]

# Write checkTest csv -----------------------------------------------------
dir.create("1MGS_check", showWarnings=F)









