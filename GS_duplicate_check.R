# Duplicate checks on 1M GeoSurvey cropland and human settlement observations.
# J. Chen, January 2015

# Required packages
require(downloader)

# Data download -----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("1MGS_check", showWarnings=F)
dat_dir <- "./1MGS_check"

download("https://www.dropbox.com/s/trxtidmym7slpdk/2015-01-05_21-28-06_062585.zip?dl=0", "./1MGS_check/2015-01-05_21-28-06_062585.zip", mode="wb")
unzip("./1MGS_check/2015-01-05_21-28-06_062585.zip", exdir="./1MGS_check", junkpaths=T, overwrite=T)
geosurvey_a1 <- read.csv(paste(dat_dir, "/6_1m-point-survey-a1.csv", sep=""), stringsAsFactors=FALSE)
geosurvey_a2 <- read.csv(paste(dat_dir, "/7_1m-point-survey-a2.csv", sep=""), stringsAsFactors=FALSE)
geosurvey_total <- rbind(geosurvey_a1, geosurvey_a2)

# Calculate number of submissions per day ---------------------------------
day <- do.call("rbind", do.call("rbind", lapply(geosurvey_total$Time, strsplit, split=" ")))[,1]
time_counts <- aggregate(geosurvey_total$Time, by=list(day), length)
names(time_counts) <- c("day", "countsperuser")

# Plot number of submissions per day
plot(1:length(c(18:31, 1:5)), time_counts$countsperuser,  axes=FALSE, xlab="day", ylab="submissions per user", type="b")
axis(1, at=1:length(c(18:31, 1:5)), labels = FALSE)
text(1:length(c(18:31, 1:5)), par("usr")[3] - 0.2, labels = time_counts$day, srt = 90, pos = 1, xpd = TRUE, cex=0.7)
axis(2)

# Calculate number of duplicates per day ----------------------------------
geosurvey_duplicates <- duplicated(geosurvey_total)
day_duplicates <- aggregate(geosurvey_duplicates, by=list(day), sum)

# Plot number of duplicates per day
plot(1:length(c(18:31, 1:5)),day_duplicates[,2],  axes=FALSE, xlab="Day", ylab="Number of duplicates", type="b")
axis(1, at=1:length(c(18:31, 1:5)), labels = FALSE)
text(1:length(c(18:31, 1:5)), par("usr")[3] - 0.2, labels = time_counts$day, srt = 90, pos = 1.5, offset=1.5, xpd = TRUE, cex=0.7)
axis(2)

# Calculate number of duplicates per user ---------------------------------
user_duplicates <- aggregate(geosurvey_duplicates, by=list(geosurvey_total$User), sum)
users <- do.call("rbind", strsplit(user_duplicates[,1], split="@"))[, 1]

# Plot number of duplicates per user
plot(1:length(user_duplicates[,1]),user_duplicates[,2],  axes=FALSE, xlab="Worker", ylab="Number of duplicates", type="b")
axis(1, at=1:length(user_duplicates[,1]), labels = FALSE)
text(1:length(user_duplicates[,1]), par("usr")[3] - 0.2, labels = users, srt = 90, pos = 1.5, offset=1.5, xpd = TRUE, cex=0.7)
axis(2)
abline(h=0)