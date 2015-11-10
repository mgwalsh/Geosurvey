#' Evaluation and local elastic net stacking of Tanzania 1MQ GeoSurvey cropland and 
#' human settlement predictions with additional Tanzania GeoSurvey test data.
#' M.Walsh & J.Chen, April 2015

#+ Required packages
# install.packages(c("downloader","raster","rgdal","dismo","caret","glmnet")), dependencies=TRUE)
require(downloader)
require(raster)
require(rgdal)
require(dismo)
require(caret)
require(glmnet)

#+ Data downloads ----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("TZ_1MQ_data", showWarnings=F)
dat_dir <- "./TZ_1MQ_data"

# download Tanzania test data
download("https://www.dropbox.com/s/gfgjnrgllwqt79d/TZ_geos_123114.csv?dl=0", "./TZ_1MQ_data/TZ_geos_123114.csv", mode="wb")
geosv <- read.table(paste(dat_dir, "/TZ_geos_123114.csv", sep=""), header=T, sep=",")

# download Tanzania 1MQ prediction grids (~21.1 Mb) and stack in raster
download("https://www.dropbox.com/s/w8l41t5muc1rr4j/TZ_1MQ_preds.zip?dl=0", "./TZ_1MQ_data/TZ_1MQ_preds.zip", mode="wb")
unzip("./TZ_1MQ_data/TZ_1MQ_preds.zip", exdir="./TZ_1MQ_data", overwrite=T)
glist <- list.files(path="./TZ_1MQ_data", pattern="tif", full.names=T)
grid <- stack(glist)

#+ Data setup --------------------------------------------------------------
# Project test data to grid CRS
geosv.proj <- as.data.frame(project(cbind(geosv$Lon, geosv$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geosv.proj) <- c("x","y")
geosv <- cbind(geosv, geosv.proj)
coordinates(geosv) <- ~x+y
projection(geosv) <- projection(grid)

# Extract gridded variables to test data observations
gsexv <- data.frame(coordinates(geosv), geosv$CRP, geosv$HSP, extract(grid, geosv))
gsexv <- na.omit(gsexv)
colnames(gsexv)[3:4] <- c("CRP", "HSP")

#+ 1MQ classifier performance evaluation ----------------------------------
# Cropland boosting classifier
gbmcrp <- subset(gsexv, CRP=="Y", select=c(CRP_gbm))
gbmcra <- subset(gsexv, CRP=="N", select=c(CRP_gbm))
gbmcrp.eval <- evaluate(p=gbmcrp[,1], a=gbmcra[,1]) ## calculate ROC's on test set <dismo>
gbmcrp.eval
plot(gbmcrp.eval, "ROC")

# Cropland neural network classifier
nncrp <- subset(gsexv, CRP=="Y", select=c(CRP_nn))
nncra <- subset(gsexv, CRP=="N", select=c(CRP_nn))
nncrp.eval <- evaluate(p=nncrp[,1], a=nncra[,1]) ## calculate ROC's on test set <dismo>
nncrp.eval
plot(nncrp.eval, "ROC")

# Cropland random forest classifier
rfcrp <- subset(gsexv, CRP=="Y", select=c(CRP_rf))
rfcra <- subset(gsexv, CRP=="N", select=c(CRP_rf))
rfcrp.eval <- evaluate(p=rfcrp[,1], a=rfcra[,1]) ## calculate ROC's on test set <dismo>
rfcrp.eval
plot(rfcrp.eval, "ROC")

# Cropland 1MQ ensemble classifier
enscrp <- subset(gsexv, CRP=="Y", select=c(CRP_ens))
enscra <- subset(gsexv, CRP=="N", select=c(CRP_ens))
enscrp.eval <- evaluate(p=enscrp[,1], a=enscra[,1]) ## calculate ROC's on test set <dismo>
enscrp.eval
plot(enscrp.eval, "ROC")
enscrp.thld <- threshold(enscrp.eval, "spec_sens") ## TPR+TNR threshold for classification
CRP_ens_mask <- grid$CRP_ens > enscrp.thld
plot(CRP_ens_mask, axes = F, legend = F)

# Building/rural settlement boosting classifier
gbmhsp <- subset(gsexv, HSP=="Y", select=c(RSP_gbm))
gbmhsa <- subset(gsexv, HSP=="N", select=c(RSP_gbm))
gbmhsp.eval <- evaluate(p=gbmhsp[,1], a=gbmhsa[,1]) ## calculate ROC's on test set <dismo>
gbmhsp.eval
plot(gbmhsp.eval, "ROC")

# Building/rural settlement neural network classifier
nnhsp <- subset(gsexv, HSP=="Y", select=c(RSP_nn))
nnhsa <- subset(gsexv, HSP=="N", select=c(RSP_nn))
nnhsp.eval <- evaluate(p=nnhsp[,1], a=nnhsa[,1]) ## calculate ROC's on test set <dismo>
nnhsp.eval
plot(nnhsp.eval, "ROC")

# Building/rural settlement random forest classifier
rfhsp <- subset(gsexv, HSP=="Y", select=c(RSP_rf))
rfhsa <- subset(gsexv, HSP=="N", select=c(RSP_rf))
rfhsp.eval <- evaluate(p=rfhsp[,1], a=rfhsa[,1]) ## calculate ROC's on test set <dismo>
rfhsp.eval
plot(rfhsp.eval, "ROC")

# Building/rural settlement 1MQ ensemble classifier
enshsp <- subset(gsexv, HSP=="Y", select=c(RSP_ens))
enshsa <- subset(gsexv, HSP=="N", select=c(RSP_ens))
enshsp.eval <- evaluate(p=enshsp[,1], a=enshsa[,1]) ## calculate ROC's on test set <dismo>
enshsp.eval
plot(enshsp.eval, "ROC")
enshsp.thld <- threshold(enshsp.eval, "spec_sens") ## TPR+TNR threshold for classification
RSP_ens_mask <- grid$RSP_ens > enshsp.thld
plot(RSP_ens_mask, axes = F, legend = F)

#+ Local classifier (re)stacking ------------------------------------------
# 10-fold CV
lcs <- trainControl(method = "cv", number = 10, classProbs = T)

# presence/absence of Cropland (CRP, present = Y, absent = N)
CRP.lcs <- train(CRP ~ CRP_gbm + CRP_nn + CRP_rf, data = gsexv,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "Accuracy",
                 trControl = lcs)
CRP.lcs
crp.pred <- predict(CRP.lcs, gsexv, type="prob")
crp.test <- cbind(gsexv, crp.pred)
lcscrp <- subset(crp.test, CRP=="Y", select=c(Y))
lcscra <- subset(crp.test, CRP=="N", select=c(Y))
lcscrp.eval <- evaluate(p=lcscrp[,1], a=lcscra[,1]) ## calculate ROC's on test set <dismo>
lcscrp.eval
plot(lcscrp.eval, "ROC") ## plot ROC curve
lcscrp.thld <- threshold(lcscrp.eval, 'spec_sens') ## max(sens+spec) threshold for classification
CRP_lcs <- predict(grid, CRP.lcs, type="prob") ## spatial prediction
plot(1-CRP_lcs, axes = F)
CRP_lcs_mask <- 1-CRP_lcs > lcscrp.thld
plot(CRP_lcs_mask, axes = F, legend = F)

# Cropland mask likelihood ratios
crp.cm <- table(crp.test$Y > lcscrp.thld, crp.test$CRP)
interval <- 0.95
a <- crp.cm[2, 2] ## tp
b <- crp.cm[2, 1] ## fp
c <- crp.cm[1, 2] ## fn
d <- crp.cm[1, 1] ## tn
crp.prev <- (a+c)/(a+b+c+d) ## prevalence = pretest probability
crp.spec <- d/(b+d) ## classifier specificity
crp.sens <- a/(a+c) ## classifier sensitivity 
lr.pos <- crp.sens/(1 - crp.spec) ## positive likelihood ratio (LR+)
sigma2.pos <- (1/a) - (1/(a+c)) + (1/b) - (1/(b+d))
lower.pos <- lr.pos * exp(-qnorm(1-((1-interval)/2))*sqrt(sigma2.pos)) ## lower LR+ 95% PI
upper.pos <- lr.pos * exp(qnorm(1-((1-interval)/2))*sqrt(sigma2.pos)) ## upper LR+ 95% PI
list(lower.pos=lower.pos, lr.pos=lr.pos, upper.pos=upper.pos)
lr.neg <- (1 - sens)/spec ## negative likelihood ratio (LR-)
sigma2.neg <- (1/c) - (1/(a+c)) + (1/d) - (1/(b+d))
lower.neg <- lr.neg * exp(-qnorm(1-((1-interval)/2))*sqrt(sigma2.neg)) ## lower LR- 95% PI
upper.neg <- lr.neg * exp(qnorm(1-((1-interval)/2))*sqrt(sigma2.neg)) ## upper LR- 95% PI
list(lower.neg=lower.neg, lr.neg=lr.neg, upper.neg=upper.neg)

# presence/absence of Buildings/rural settlements (HSP, present = Y, absent = N)
RSP.lcs <- train(HSP ~ RSP_gbm + RSP_nn + RSP_rf, data = gsexv,
                 family = "binomial", 
                 method = "glmnet",
                 metric = "Accuracy",
                 trControl = lcs)
RSP.lcs
rsp.pred <- predict(RSP.lcs, gsexv, type="prob")
rsp.test <- cbind(gsexv, rsp.pred)
lcsrsp <- subset(rsp.test, HSP=="Y", select=c(Y))
lcsrsa <- subset(rsp.test, HSP=="N", select=c(Y))
lcsrsp.eval <- evaluate(p=lcsrsp[,1], a=lcsrsa[,1]) ## calculate ROC's on test set <dismo>
lcsrsp.eval
plot(lcsrsp.eval, "ROC") ## plot ROC curve
lcsrsp.thld <- threshold(lcsrsp.eval, 'spec_sens') ## max(sens+spec) threshold for classification
RSP_lcs <- predict(grid, RSP.lcs, type="prob") ## spatial prediction
plot(1-RSP_lcs, axes = F)
RSP_lcs_mask <- 1-RSP_lcs > lcsrsp.thld
plot(RSP_lcs_mask, axes = F, legend = F)

#+ Write spatial predictions -----------------------------------------------
dir.create("TZ_1MQ_results", showWarnings=F)
CRP_lcs_pred <- stack(CRP_ens_mask, 1-CRP_lcs, CRP_lcs_mask)
RSP_lcs_pred <- stack(RSP_ens_mask, 1-RSP_lcs, RSP_lcs_mask)
writeRaster(CRP_lcs_pred, filename="./TZ_1MQ_results/TZ_crp_pred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
writeRaster(RSP_lcs_pred, filename="./TZ_1MQ_results/TZ_rsp_pred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
