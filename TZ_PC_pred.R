# Stacked predictions of Tanzania presence/absence of croplands with buildings (priority croplands)
# M. Walsh, January 2018

# Required packages
# install.packages(c("devtools","caret","MASS","randomForest","gbm","nnet","glmnet","plyr","doParallel","dismo")), dependencies=T)
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
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$PC, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
cp_cal <- gs_cal$PC ## Cropland & buildings present? (Y/N)

# raster calibration features
gf_cal <- gs_cal[,11:48] ## grid covariates

# Central place theory model <glm> -----------------------------------------
# select central place variables
gf_cpv <- gs_cal[c(14:23,42)] ## central-place covariates & slope

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
PC1.gl <- train(gf_cpv, cp_cal, 
                method = "glmStepAIC",
                family = "binomial",
                preProc = c("center","scale"), 
                trControl = tc,
                metric ="ROC")

# model outputs & predictions
summary(PC1.gl)
print(PC1.gl) ## ROC's accross cross-validation
pcg1.pred <- predict(grids, PC1.gl, type = "prob") ## spatial predictions

stopCluster(mc)

# GLM with all covariates -------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
PC2.gl <- train(gf_cal, cp_cal, 
                method = "glmStepAIC",
                family = "binomial",
                preProc = c("center","scale"), 
                trControl = tc,
                metric ="ROC")

# model outputs & predictions
summary(PC2.gl)
print(PC2.gl) ## ROC's accross cross-validation
pcg2.pred <- predict(grids, PC2.gl, type = "prob") ## spatial predictions

stopCluster(mc)

# Random forest <randomForest> --------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
tg <- expand.grid(mtry=seq(1, 10, by=1))
PC.rf <- train(gf_cal, cp_cal,
               preProc = c("center","scale"),
               method = "rf",
               ntree = 501,
               metric = "ROC",
               tuneGrid = tg,
               trControl = tc)

# model outputs & predictions
print(PC.rf) ## ROC's accross tuning parameters
plot(varImp(PC.rf)) ## relative variable importance
pcrf.pred <- predict(grids, PC.rf, type = "prob") ## spatial predictions

stopCluster(mc)

# Generalized boosting <gbm> ----------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, summaryFunction = twoClassSummary,
                   allowParallel = T)

# model training
PC.gb <- train(gf_cal, cp_cal, 
               method = "gbm", 
               preProc = c("center", "scale"),
               trControl = tc,
               metric = "ROC")

# model outputs & predictions
print(PC.gb) ## ROC's accross tuning parameters
plot(varImp(PC.gb)) ## relative variable importance
pcgb.pred <- predict(grids, PC.gb, type = "prob") ## spatial predictions

stopCluster(mc)

# Neural network <nnet> ---------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
PC.nn <- train(gf_cal, cp_cal, 
               method = "nnet",
               preProc = c("center","scale"), 
               trControl = tc,
               metric ="ROC")

# model outputs & predictions
print(PC.nn) ## ROC's accross tuning parameters
plot(varImp(PC.nn)) ## relative variable importance
pcnn.pred <- predict(grids, PC.nn, type = "prob") ## spatial predictions

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
preds <- stack(1-pcg1.pred, 1-pcg2.pred, 1-pcrf.pred, 1-pcgb.pred, 1-pcnn.pred)
names(preds) <- c("gl1","gl2","rf", "gb","nn")
plot(preds, axes = F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(preds)
gspred <- extract(preds, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# stacking model validation labels and features
cp_val <- gspred$PC ## subset validation labels
gf_val <- gspred[,50:54] ## subset validation features

# Model stacking ----------------------------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
PC.st <- train(gf_val, cp_val,
               method = "glmnet",
               family = "binomial",
               metric = "ROC",
               trControl = tc)

# model outputs & predictions
print(PC.st)
plot(varImp(PC.st))
pcst.pred <- predict(preds, PC.st, type = "prob") ## spatial predictions
plot(1-pcst.pred, axes = F)

stopCluster(mc)

# Receiver-operator characteristics ---------------------------------------
cp_pre <- predict(PC.st, gf_val, type="prob")
cp_val <- cbind(cp_val, cp_pre)
cpp <- subset(cp_val, cp_val=="Y", select=c(Y))
cpa <- subset(cp_val, cp_val=="N", select=c(Y))
cp_eval <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC's on test set
plot(cp_eval, 'ROC') ## plot ROC curve

# complete-set ROC
# extract model predictions
coordinates(gsdat) <- ~x+y
projection(gsdat) <- projection(preds)
gspred <- extract(preds, gsdat)
gspred <- as.data.frame(cbind(gsdat, gspred))
write.csv(gsdat, "./Results/TZ_PC_pred.csv", row.names = F) ## write dataframe

# stacking model labels and features
cp_all <- gspred$PC ## subset validation labels
gf_all <- gspred[,50:54] ## subset model predictions

# ROC calculation
cp_pre <- predict(PC.st, gf_all, type="prob")
cp_all <- cbind(cp_all, cp_pre)
cpp <- subset(cp_all, cp_all=="Y", select=c(Y))
cpa <- subset(cp_all, cp_all=="N", select=c(Y))
cp_eall <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC on complete set
cp_eall
plot(cp_eall, 'ROC') ## plot ROC curve

# Generate mask -----------------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on ROC
r <- matrix(c(0, t[,2], 0, t[,2], 1, 1), ncol=3, byrow = T) ## set threshold value <spec_sens>
mask <- reclassify(1-pcst.pred, r) ## reclassify stacked predictions
plot(mask, axes=F)

# Write prediction files --------------------------------------------------
pcpreds <- stack(preds, 1-pcst.pred, mask)
names(pcpreds) <- c("gl1","gl2","rf","gb","nn","st","mk")
writeRaster(pcpreds, filename="./Results/TZ_pcpreds_2017.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Prediction map widget ---------------------------------------------------
require(leaflet)
require(htmlwidgets)

# ensemble prediction map 
pred <- 1-pcst.pred ## GeoSurvey ensemble probability

# set color pallet
pal <- colorBin("Greens", domain = 0:1) 

# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(pred, colors = pal, opacity = 0.5) %>%
  addLegend(pal = pal, values = values(pred), title = "Priority croplands")
w ## plot widget 

# save widget
saveWidget(w, 'TZ_PC_prob.html', selfcontained = T)
