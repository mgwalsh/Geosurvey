#' Geosurvey duplicate and accuracy screening checks on cropland and human settlement observations.
#' J. Chen & M. Walsh, February 2015

# Required packages
require(downloader)
require(arm)

#+ Data download ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("1MGS_data", showWarnings=F)
dat_dir <- "./1MGS_data"

# 1MQ GeoSurvey
download("https://www.dropbox.com/s/tnc1wwmi8a4h6b1/2015-01-20_18-00-02_248261.zip?dl=0", "./1MGS_data/2015-01-20_18-00-02_248261.zip", mode="wb")
unzip("./1MGS_data/2015-01-20_18-00-02_248261.zip", exdir="./1MGS_data", junkpaths=T, overwrite=T)
geosurvey_a1 <- read.csv(paste(dat_dir, "/6_1m-point-survey-a1.csv", sep=""), stringsAsFactors=FALSE)
geosurvey_a2 <- read.csv(paste(dat_dir, "/7_1m-point-survey-a2.csv", sep=""), stringsAsFactors=FALSE)
geosurvey_total <- rbind(geosurvey_a1[ ,1:6], geosurvey_a2[ ,1:6])

# Validation dataset
download("https://www.dropbox.com/s/pt86fr3ko379f8h/1MQ_validation_data.csv?dl=0", "./1MGS_data/1MQ_validation_data.csv", mode ="wb")
geosv <- read.table(paste(dat_dir, "/1MQ_validation_data.csv", sep=""), header=T, sep=",")

#+ Remove coordinate duplicates -------------------------------------------
dupindex <- which(duplicated(geosurvey_total[ ,2:3]))
dupindex_rev <- which(duplicated(geosurvey_total[ ,2:3], fromLast=TRUE))
dupindex_total <- unique(c(dupindex, dupindex_rev))
geos_nodups <- geosurvey_total[-dupindex, ]
nrow(geos_nodups)

#+ Geosurveyor (GS) accuracy assessment -----------------------------------
# Cropland observations
geov_crp <- geosv[, c(1:4, 6, 8)]
# NA values in 1MQ are "Don't know", which should also be compared to expert answers to calculate accuracy rates
geov_crp[,5] <- ifelse(is.na(geov_crp[,5]), "NA", geov_crp[,5])
geov_crp[,6] <- ifelse(is.na(geov_crp[,6]), "NA", geov_crp[,6])

# Cropland random effects model, GS vs expert
CRPcn <- as.numeric(factor(geov_crp$CRPc))
CRPn <- as.numeric(factor(geov_crp$CRP))
CRP_diff <- ifelse(CRPcn - CRPn==0, 1, 0)
geosv_narm <- geov_crp[!is.na(CRP_diff), ]
CRP_diff <- na.omit(CRP_diff)
CRP.glmer <- glmer(CRP_diff~1+(1|geosv_narm$User), family=binomial)
display(CRP.glmer)

# If a GS's estimated agreement rate is <50%, the data for that GS are removed from further analyses
coef(CRP.glmer)
nosample_CRP <- rownames(coef(CRP.glmer)[[1]])[coef(CRP.glmer)[[1]][,1]<0]
crp_data <- geos_nodups[!geos_nodups$User%in%nosample_CRP, ]
nrow(crp_data)

# Building / rural settlement observations
geov_hs <- geosv[, c(1:4, 5, 7)]
# NA values in 1MQ are "Don't know", which should also be compared to expert answers to calculate accuracy rates
geov_hs[,5] <- ifelse(is.na(geov_hs[,5]), "NA", geov_hs[,5])
geov_hs[,6] <- ifelse(is.na(geov_hs[,6]), "NA", geov_hs[,6])

# Building / rural settlement random effects model, GS vs expert
HSPcn <- as.numeric(factor(geov_hs$HSPc))
HSPn <- as.numeric(factor(geov_hs$HSP))
HSP_diff <- ifelse(HSPcn - HSPn==0, 1, 0)
geosv_narm <- geosv[!is.na(HSP_diff), ]
HSP_diff <- na.omit(HSP_diff)
HSP.glmer <- glmer(HSP_diff~1+(1|geosv_narm$User), family=binomial)
display(HSP.glmer)

# If a GS's estimated agreement rate is <50%, the data for that GS are removed from further analyses
coef(HSP.glmer)
nosample_HSP <- rownames(coef(HSP.glmer)[[1]])[coef(HSP.glmer)[[1]][,1]<0]
hsp_data <- geos_nodups[!geos_nodups$User%in%nosample_HSP, ]
nrow(hsp_data)

#+ Write cleaned files -----------------------------------------------------
write.csv(crp_data[ ,1:4,6], paste(dat_dir, "/CRP_cleaned.csv", sep=""), row.names=F)
write.csv(hsp_data[ ,1:5], paste(dat_dir, "/HSP_cleaned.csv", sep=""), row.names=F)
