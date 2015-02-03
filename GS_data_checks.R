# Incoming GeoSurvey data check's
# M. Walsh, J. Chen, W. Wu, February 2015

# Required packages
# install.packages(c("downloader","caret")), dependencies = T)
require(downloader)
require(caret)

# Load current GeoSurvey data ---------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("1MGS_data", showWarnings=F)
dat_dir <- "./1MGS_data"

# download GeoSurvey data
download("https://www.dropbox.com/s/jtwmf7ck9ebh7bm/1MGS_cleaned.csv.zip?dl=0", "./1MGS_data/1MGS_cleaned.csv.zip", mode="wb")
unzip("./1MGS_data/1MGS_cleaned.csv.zip", exdir="./1MGS_data", overwrite=T)
geos <- read.table(paste(dat_dir, "/1MGS_cleaned.csv", sep=""), header=T, sep=",")

# Generate check sample ---------------------------------------------------
# set check sample randomization seed
seed <- 1385321
set.seed(seed)

# Check sample split
checkIndex <- createDataPartition(geos$User, list=F, p=0.999, times=1)
checkTest  <- geos[-checkIndex, ]
checkSample <- checkTest[ ,2:3]

# Write checkSample csv ---------------------------------------------------
dir.create("1MGS_check", showWarnings=F)
write.csv(checkTest, "./1MGS_check/Check_test.csv", row.names=F)
write.csv(checkSample, "./1MGS_check/Check_sample.csv", row.names=F)









