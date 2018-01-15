# Stacked predictions of Tanzania cropland observations
# M. Walsh, September 2017

# Required packages
# install.packages(c("devtools","caret","MASS","randomForest","gbm","nnet","glmnet","plyr","doParallel","dismo")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(MASS)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(glmnet)
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
seed <- 123581
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$CP, p = 4/5, list = FALSE, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
cp_cal <- gs_cal$CP ## Croplands present? (Y/N)

# Raster calibration features
gf_cal <- gs_cal[,10:47] ## grid covariates

# Central place theory model <glm> -----------------------------------------
# select central place variables
gf_cpv <- gs_cal[c(13:22,41)] ## central-place covariates & slope

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats=5, classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
CP.gl <- train(gf_cpv, cp_cal, 
               method = "glmStepAIC",
               family = "binomial",
               preProc = c("center","scale"), 
               trControl = tc,
               metric ="ROC")

# model outputs & predictions
print(CP.gl) ## ROC's accross cross-validation
plot(varImp(CP.gl)) ## relative variable importance
confusionMatrix(CP.gl) ## cross-validation performance
cpgl.pred <- predict(grids, CP.gl, type = "prob") ## spatial predictions

stopCluster(mc)

# Random forest <randomForest> --------------------------------------------
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

# Model stacking setup ----------------------------------------------------
preds <- stack(1-cpgl.pred, 1-cprf.pred, 1-cpgb.pred, 1-cpnn.pred)
names(preds) <- c("gl","rf", "gb","nn")
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

# Receiver-operator characteristics <dismo> -------------------------------
# validation-set ROC
cp_pre <- predict(CP.st, gf_val, type="prob")
cp_val <- cbind(cp_val, cp_pre)
cpp <- subset(cp_val, cp_val=="Y", select=c(Y))
cpa <- subset(cp_val, cp_val=="N", select=c(Y))
cp_eval <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC on test set
cp_eval
plot(cp_eval, 'ROC') ## plot ROC curve

# complete-set ROC
# extract model predictions
coordinates(gsdat) <- ~x+y
projection(gsdat) <- projection(preds)
gspred <- extract(preds, gsdat)
gspred <- as.data.frame(cbind(gsdat, gspred))
write.csv(gsdat, "./Results/TZ_CP_pred.csv", row.names = F) ## write dataframe

# stacking model labels and features
cp_all <- gspred$CP ## subset validation labels
gf_all <- gspred[,46:49] ## subset validation features

# ROC calculation
cp_pre <- predict(CP.st, gf_all, type="prob")
cp_all <- cbind(cp_all, cp_pre)
cpp <- subset(cp_all, cp_all=="Y", select=c(Y))
cpa <- subset(cp_all, cp_all=="N", select=c(Y))
cp_eall <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC on complete set
cp_eall
plot(cp_eall, 'ROC') ## plot ROC curve

# Generate cropland mask --------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on validation-set ROC
r <- matrix(c(0, t[,2], 0, t[,2], 1, 1), ncol=3, byrow = T) ## set threshold value <spec_sens>
cp_all$mask <- ifelse(cp_all$Y > t[,2], "Y", "N") ## validation-set based classification
table(cp_all$cp_all, cp_all$mask) # confusion/classification table
mask <- reclassify(1-cpst.pred, r) ## classify stacked spatial predictions to produce mask
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
      addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
      addRasterImage(pred, colors = pal, opacity = 0.3) %>%
      addLegend(pal = pal, values = values(pred), title = "Cropland prob")
w ## plot widget 

# save widget
saveWidget(w, 'TZ_CP_prob.html', selfcontained = T)
