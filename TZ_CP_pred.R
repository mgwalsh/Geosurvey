# Stacked predictions of Tanzania cropland observations
# M. Walsh, September 2017

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
gsIndex <- createDataPartition(gsdat$CP, p = 4/5, list = FALSE, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
cp_cal <- gs_cal$CP ## Croplands present? (Y/N)

# Raster calibration features
gf_cal <- gs_cal[,7:44] ## grid covariates

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
CP.rf <- train(gf_cal, cp_cal,
               preProc = c("center","scale"),
               method = "rf",
               ntree = 501,
               metric = "ROC",
               tuneGrid = tg,
               trControl = tc)

# model outputs & predictions
print(CP.rf) ## ROC's accross tuning parameters
plot(varImp(CP.rf)) ## relative variable importance
confusionMatrix(CP.rf) ## cross-validation performance
cprf.pred <- predict(grids, CP.rf, type = "prob") ## spatial predictions

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
CP.gb <- train(gf_cal, cp_cal, 
               method = "gbm", 
               preProc = c("center", "scale"),
               trControl = tc,
               metric = "ROC")

# model outputs & predictions
print(CP.gb) ## ROC's accross tuning parameters
plot(varImp(CP.gb)) ## relative variable importance
confusionMatrix(CP.gb) ## cross-validation performance
cpgb.pred <- predict(grids, CP.gb, type = "prob") ## spatial predictions

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
CP.nn <- train(gf_cal, cp_cal, 
               method = "nnet",
               preProc = c("center","scale"), 
               trControl = tc,
               metric ="ROC")

# model outputs & predictions
print(CP.nn) ## ROC's accross tuning parameters
plot(varImp(CP.nn)) ## relative variable importance
confusionMatrix(CP.nn) ## cross-validation performance
cpnn.pred <- predict(grids, CP.nn, type = "prob") ## spatial predictions

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
CP.rr <- train(gf_cal, cp_cal, 
               method = "glmnet",
               family = "binomial",
               preProc = c("center","scale"), 
               trControl = tc,
               metric ="ROC")

# model outputs & predictions
print(CP.rr) ## ROC's accross tuning parameters
plot(varImp(CP.rr)) ## relative variable importance
confusionMatrix(CP.rr) ## cross-validation performance
cprr.pred <- predict(grids, CP.rr, type = "prob") ## spatial predictions

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
preds <- stack(1-cprf.pred, 1-cpgb.pred, 1-cpnn.pred, 1-cprr.pred)
names(preds) <- c("rf","gb", "nn","rr")
plot(preds, axes=F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
cp_val <- gspred$CP ## subset validation labels
gf_val <- gspred[,46:49] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats = 5, classProbs = TRUE, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
CP.st <- train(gf_val, cp_val,
               method = "glmnet",
               family = "binomial",
               metric = "ROC",
               trControl = tc)

# model outputs & predictions
print(CP.st)
confusionMatrix(CP.st)
plot(varImp(CP.st))
cpst.pred <- predict(preds, CP.st, type = "prob") ## spatial predictions
plot(1-cpst.pred, axes=F)

stopCluster(mc)

# Receiver-operator characteristics ---------------------------------------
cp_pre <- predict(CP.st, gf_val, type="prob")
cp_val <- cbind(cp_val, cp_pre)
cpp <- subset(cp_val, cp_val=="Y", select=c(Y))
cpa <- subset(cp_val, cp_val=="N", select=c(Y))
cp_eval <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC on test set
plot(cp_eval, 'ROC') ## plot ROC curve

# Generate cropland mask --------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on ROC
r <- matrix(c(0, t[,2], 0, t[,2], 1, 1), ncol=3, byrow=TRUE) ## set threshold value <spec_sens>
mask <- reclassify(1-cpst.pred, r) ## reclassify stacked predictions
plot(mask, axes=F)

# Write prediction files --------------------------------------------------
cppreds <- stack(preds, 1-cpst.pred, mask)
names(cppreds) <- c("cprf","cpgb","cpnn","cprr","cpst","cpmk")
writeRaster(cppreds, filename="./Results/TZ_cppreds_2017.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Prediction map widget ---------------------------------------------------
require(leaflet)
require(htmlwidgets)

# ensemble prediction map 
pred <- 1-cpst.pred ## GeoSurvey ensemble probability

# set color pallet
pal <- colorBin("Greens", domain = 0:1) 

# render map
w <- leaflet() %>% 
      addTiles() %>% # default basemap: OSM
      addRasterImage(pred, colors = pal, opacity = 0.5) %>%
      addLegend(pal = pal, values = values(pred), title = "Cropland prob")
w ## plot widget 

# save widget
saveWidget(w, 'TZ_CP_prob.html', selfcontained = T)
