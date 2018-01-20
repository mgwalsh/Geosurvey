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
gf_cpv <- gs_cal[c(14:23,42)] ## central-place covariates & slope
gf_cal <- gs_cal[,11:48] ## grid covariates

# Generalized linear models <MASS> ----------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Central-place variables only <gl1>
# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gl1 <- train(gf_cpv, cp_cal, 
             method = "glmStepAIC",
             family = "binomial",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model predictions
summary(gl1)
print(gl1) ## ROC's accross cross-validation
gl1.pred <- predict(grids, gl1, type = "prob") ## spatial predictions
stopCluster(mc)

# All covariates model <gl2>
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
gl2 <- train(gf_cal, cp_cal, 
             method = "glmStepAIC",
             family = "binomial",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model predictions
summary(gl2)
print(gl2) ## ROC's accross cross-validation
gl2.pred <- predict(grids, gl2, type = "prob") ## spatial predictions
stopCluster(mc)

# Random Forest models <randomForest> -------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Central-place variables only <rf1>
# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)
tg <- expand.grid(mtry=seq(1, 10, by=1)) ## tuning grid parameters

# model training
rf1 <- train(gf_cpv, cp_cal,
             preProc = c("center","scale"),
             method = "rf",
             ntree = 501,
             metric = "ROC",
             tuneGrid = tg,
             trControl = tc)

# model predictions
print(rf1) ## ROC's accross tuning parameters
plot(varImp(rf1)) ## relative variable importance
rf1.pred <- predict(grids, rf1, type = "prob") ## spatial predictions
stopCluster(mc)

# All covariates model <rf2>
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)
tg <- expand.grid(mtry=seq(1, 10, by=1)) ## tuning grid parameters

# model training
rf2 <- train(gf_cal, cp_cal,
             preProc = c("center","scale"),
             method = "rf",
             ntree = 501,
             metric = "ROC",
             tuneGrid = tg,
             trControl = tc)

# model predictions
print(rf2) ## ROC's accross tuning parameters
plot(varImp(rf2)) ## relative variable importance
rf2.pred <- predict(grids, rf2, type = "prob") ## spatial predictions
stopCluster(mc)

# Generalized boosting models <gbm> ---------------------------------------
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Central-place variables only <gb1>
# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, summaryFunction = twoClassSummary,
                   allowParallel = T)

# model training
gb1 <- train(gf_cpv, cp_cal, 
             method = "gbm", 
             preProc = c("center", "scale"),
             trControl = tc,
             metric = "ROC")

# model predictions
print(gb1) ## ROC's accross tuning parameters
plot(varImp(gb1)) ## relative variable importance
gb1.pred <- predict(grids, gb1, type = "prob") ## spatial predictions
stopCluster(mc)

# All covariates model <gb2>
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, summaryFunction = twoClassSummary,
                   allowParallel = T)

# model training
gb2 <- train(gf_cal, cp_cal, 
             method = "gbm", 
             preProc = c("center", "scale"),
             trControl = tc,
             metric = "ROC")

# model predictions
print(gb2) ## ROC's accross tuning parameters
plot(varImp(gb2)) ## relative variable importance
gb2.pred <- predict(grids, gb2, type = "prob") ## spatial predictions
stopCluster(mc)

# Neural network models <nnet> --------------------------------------------
# Central-place variables only <nn1>
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
nn1 <- train(gf_cpv, cp_cal, 
             method = "nnet",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model predictions
print(nn1) ## ROC's accross tuning parameters
plot(varImp(nn1)) ## relative variable importance
nn1.pred <- predict(grids, nn1, type = "prob") ## spatial predictions
stopCluster(mc)

# All covariates model <nn2>
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T,
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
nn2 <- train(gf_cal, cp_cal, 
             method = "nnet",
             preProc = c("center","scale"), 
             trControl = tc,
             metric ="ROC")

# model predictions
print(nn2) ## ROC's accross tuning parameters
plot(varImp(nn2)) ## relative variable importance
nn2.pred <- predict(grids, nn2, type = "prob") ## spatial predictions
stopCluster(mc)

# Model stacking setup ----------------------------------------------------
preds <- stack(1-gl1.pred, 1-gl2.pred, 1-rf1.pred, 1-rf2.pred, 1-gb1.pred, 1-gb2.pred, 1-nn1.pred, 1-nn2.pred)
names(preds) <- c("gl1", "gl2", "rf1", "rf2", "gb1", "gb2", "nn1", "nn2")
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
