#' Geolocation duplicate checks on 1M GeoSurvey cropland and human settlement observations.
#' J.Chen & M.Walsh, January 2015

# Required packages
require(downloader)

#+ Data download ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("1MGS_check", showWarnings=F)
dat_dir <- "./1MGS_check"

download("https://www.dropbox.com/s/joxjr1rz32uk2rm/2015-01-14_18-00-01_784908.zip?dl=0", "./1MGS_check/2015-01-14_18-00-01_784908.zip", mode="wb")
unzip("./1MGS_check/2015-01-14_18-00-01_784908.zip", exdir="./1MGS_check", junkpaths=T, overwrite=T)
geosurvey_a1 <- read.csv(paste(dat_dir, "/6_1m-point-survey-a1.csv", sep=""), stringsAsFactors=FALSE)
geosurvey_a2 <- read.csv(paste(dat_dir, "/7_1m-point-survey-a2.csv", sep=""), stringsAsFactors=FALSE)
geosurvey_total <- rbind(geosurvey_a1[ ,1:6], geosurvey_a2[ ,1:6])

#+ Calculate number of submissions per day --------------------------------
day <- do.call("rbind", do.call("rbind", lapply(geosurvey_total$Time, strsplit, split=" ")))[,1]
submis <- aggregate(geosurvey_total$Time, by=list(day), length)
names(submis) <- c("Date", "Submissions")
submis$Date <- as.Date(submis$Date)
submis

# Plot total number of submissions per day
plot(Submissions/1000~Date, submis, type="b", ylab="Total submissions ('000's)")

#+ Calculate number of duplicates per day ---------------------------------
geosurvey_duplicates <- duplicated(geosurvey_total[ ,2:3])
geodups <- aggregate(geosurvey_duplicates, by=list(day), sum)
names(geodups) <- c("Date", "Duplicates")
geodups$Date <- as.Date(geodups$Date)
geodups

# Plot number of duplicate observations per day
plot(Duplicates~Date, geodups, type="b", ylab="Total Geo-duplicates")

#+ Calculate number of duplicates per user --------------------------------
user_dups <- aggregate(geosurvey_duplicates, by=list(geosurvey_total$User), sum)
names(user_dups) <- c("Email", "Duplicates")
user_dups

#+ Remove coordinate duplicates -------------------------------------------
dupindex <- which(duplicated(geosurvey_total[ ,2:3]))
geos_nodups <- geosurvey_total[-dupindex, ]
nrow(geos_nodups)

# Write cleaned file
write.csv(geos_nodups, paste(dat_dir, "/geos_cleaned.csv", sep=""), row.names=F)
