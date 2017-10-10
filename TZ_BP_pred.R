# Stacked predictions of Tanzania presence/absence of buildings observations
# M. Walsh, October 2017

# Required packages
# install.packages(c("devtools","caret","plyr","doParallel")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(plyr)
  require(doParallel)
  require(dismo)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Geosurvey/blob/master/TZ_GS_data.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/blob/master/TZ_GS_data.R"
# source_url(SourceURL)
rm(list=setdiff(ls(), c("gsdat","grids","glist"))) ## scrub extraneous objects in memory

# set calibration/validation set randomization seed
seed <- 1385321
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$BP, p = 4/5, list = FALSE, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
cp_cal <- gs_cal$BP ## Buildings present? (Y/N)

# Raster calibration features
gf_cal <- gs_cal[,7:37] ## grid features

# Random forest <randomForest> --------------------------------------------
require(randomForest)

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
tg <- expand.grid(mtry=seq(1, 10, by=1))
BP.rf <- train(gf_cal, cp_cal,
               preProc = c("center","scale"),
               method = "rf",
               ntree = 501,
               metric = "ROC",
               tuneGrid = tg,
               trControl = tc)

# model outputs & predictions
print(BP.rf) ## ROC's accross tuning parameters
plot(varImp(BP.rf)) ## relative variable importance
confusionMatrix(BP.rf) ## cross-validation performance
bprf.pred <- predict(grids, BP.rf, type = "prob") ## spatial predictions

stopCluster(mc)

# Generalized boosting <gbm> ----------------------------------------------
require(gbm)

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE, summaryFunction = twoClassSummary,
                   allowParallel = T)

# model training
BP.gb <- train(gf_cal, cp_cal, 
               method = "gbm", 
               preProc = c("center", "scale"),
               trControl = tc,
               metric = "ROC")

# model outputs & predictions
print(BP.gb) ## ROC's accross tuning parameters
plot(varImp(BP.gb)) ## relative variable importance
confusionMatrix(BP.gb) ## cross-validation performance
bpgb.pred <- predict(grids, BP.gb, type = "prob") ## spatial predictions

stopCluster(mc)

# Neural network <nnet> ---------------------------------------------------
require(nnet)

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
BP.nn <- train(gf_cal, cp_cal, 
               method = "nnet",
               preProc = c("center","scale"), 
               trControl = tc,
               metric ="ROC")

# model outputs & predictions
print(BP.nn) ## ROC's accross tuning parameters
plot(varImp(BP.nn)) ## relative variable importance
confusionMatrix(BP.nn) ## cross-validation performance
bpnn.pred <- predict(grids, BP.nn, type = "prob") ## spatial predictions

stopCluster(mc)

# Regularized regression <glmnet> -----------------------------------------
require(glmnet)

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats=5, classProbs = TRUE,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
BP.rr <- train(gf_cal, cp_cal, 
               method = "glmnet",
               family = "binomial",
               preProc = c("center","scale"), 
               trControl = tc,
               metric ="ROC")

# model outputs & predictions
print(BP.rr) ## ROC's accross tuning parameters
plot(varImp(BP.rr)) ## relative variable importance
confusionMatrix(BP.rr) ## cross-validation performance
bprr.pred <- predict(grids, BP.rr, type = "prob") ## spatial predictions

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
preds <- stack(1-bprf.pred, 1-bpgb.pred, 1-bpnn.pred, 1-bprr.pred)
names(preds) <- c("rf","gb", "nn","rr")
plot(preds, axes=F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
cp_val <- gspred$BP ## subset validation labels
gf_val <- gspred[,38:41] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats = 5, classProbs = TRUE, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
BP.st <- train(gf_val, cp_val,
               method = "glmnet",
               family = "binomial",
               metric = "ROC",
               trControl = tc)

# model outputs & predictions
print(BP.st)
confusionMatrix(BP.st)
plot(varImp(BP.st))
bpst.pred <- predict(preds, BP.st, type = "prob") ## spatial predictions
plot(1-bpst.pred, axes=F)

stopCluster(mc)

# Receiver-operator characteristics ---------------------------------------
cp_pre <- predict(BP.st, gf_val, type="prob")
cp_val <- cbind(cp_val, cp_pre)
cpp <- subset(cp_val, cp_val=="Y", select=c(Y))
cpa <- subset(cp_val, cp_val=="N", select=c(Y))
cp_eval <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC's on test set
plot(cp_eval, 'ROC') ## plot ROC curve

# Generate building mask --------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on ROC
r <- matrix(c(0, t[,2], 0, t[,2], 1, 1), ncol=3, byrow=TRUE) ## set threshold value <spec_sens>
mask <- reclassify(1-bpst.pred, r) ## reclassify stacked predictions
plot(mask, axes=F)

# Write prediction files --------------------------------------------------
bppreds <- stack(preds, 1-bpst.pred, mask)
names(bppreds) <- c("bprf","bpgb","bpnn","bprr","bpst","bpmk")
writeRaster(bppreds, filename="./Results/TZ_bppreds_2017.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)
