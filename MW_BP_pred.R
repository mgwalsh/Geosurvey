# Initial stacked predictions of Malawi presence/absence of buildings
# M. Walsh, March 2018

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
# Run this first: https://github.com/mgwalsh/Geosurvey/blob/master/MW_GS100_data.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/blob/master/MW_GS100_data.R"
# source_url(SourceURL)
rm(list=setdiff(ls(), c("gsdat","grids","glist"))) ## scrub extraneous objects in memory

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(gsdat$BP, p = 4/5, list = F, times = 1)
gs_cal <- gsdat[ gsIndex,]
gs_val <- gsdat[-gsIndex,]

# GeoSurvey calibration labels
cp_cal <- gs_cal$BP ## Buildings present? (Y/N)

# raster calibration features
# gf_cal <- gs_cal[,11:48] ## grid covariates

# Central place theory model <glm> -----------------------------------------
# select central place variables
gf_cpv <- gs_cal[c(10:18)] ## central-place covariates

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

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

# model outputs & predictions
summary(gl1)
print(gl1) ## ROC's accross cross-validation
gl1.pred <- predict(grids, gl1, type = "prob") ## spatial predictions
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

# <rf1> model predictions
print(rf1) ## ROC's accross tuning parameters
plot(varImp(rf1)) ## relative variable importance
rf1.pred <- predict(grids, rf1, type = "prob") ## spatial predictions
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

# <gb1> model predictions
print(gb1) ## ROC's accross tuning parameters
plot(varImp(gb1)) ## relative variable importance
gb1.pred <- predict(grids, gb1, type = "prob") ## spatial predictions
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

# <nn1> model predictions
print(nn1) ## ROC's accross tuning parameters
plot(varImp(nn1)) ## relative variable importance
nn1.pred <- predict(grids, nn1, type = "prob") ## spatial predictions
stopCluster(mc)

# Model stacking setup ----------------------------------------------------
pred1 <- stack(1-gl1.pred, 1-rf1.pred, 1-gb1.pred, 1-nn1.pred) ## central place predictions
names(pred1) <- c("gl1","rf1","gb1","nn1")
plot(pred1, axes = F)
# pred2 <- stack(1-gl2.pred, 1-rf2.pred, 1-gb2.pred, 1-nn2.pred) ## predictions with all covariates
# names(pred2) <- c("gl2","rf2","gb2","nn2")
# plot(pred2, axes = F)

# extract model predictions
coordinates(gs_val) <- ~x+y
projection(gs_val) <- projection(pred1)
gspred <- extract(pred1, gs_val)
gspred <- as.data.frame(cbind(gs_val, gspred))

# Model stacking ----------------------------------------------------------
# stacking model validation labels and features
cp_val <- gspred$BP ## validation labels
gf_val <- gspred[,20:23] ## validation features

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = T, 
                   summaryFunction = twoClassSummary, allowParallel = T)

# model training
st <- train(gf_val, cp_val,
            method = "glmnet",
            family = "binomial",
            metric = "ROC",
            trControl = tc)

# <st> model predictions
print(st)
plot(varImp(st))
st.pred <- predict(pred1, st, type = "prob") ## spatial predictions
plot(1-st.pred, axes = F)
stopCluster(mc)

# Receiver-operator characteristics ---------------------------------------
cp_pre <- predict(st, gf_val, type="prob")
cp_val <- cbind(cp_val, cp_pre)
cpp <- subset(cp_val, cp_val=="Y", select=c(Y))
cpa <- subset(cp_val, cp_val=="N", select=c(Y))
cp_eval <- evaluate(p=cpp[,1], a=cpa[,1]) ## calculate ROC's on test set
plot(cp_eval, 'ROC') ## plot ROC curve

# Generate mask -----------------------------------------------------------
t <- threshold(cp_eval) ## calculate thresholds based on ROC
r <- matrix(c(0, t[,1], 0, t[,1], 1, 1), ncol=3, byrow = T) ## set threshold value <kappa>
mask <- reclassify(1-st.pred, r) ## reclassify stacked predictions
plot(mask, axes=F, legend=F)

# Write prediction files --------------------------------------------------
pcpreds <- stack(pred1, 1-st.pred, mask)
names(pcpreds) <- c("gl1","rf1","gb1","nn1","st","mk")
writeRaster(pcpreds, filename="./Results/MW_BP100_pred.tif", datatype="FLT4S", options="INTERLEAVE=BAND", overwrite=T)

# Prediction map widget ---------------------------------------------------
# ensemble prediction map 
pred <- 1-st.pred ## GeoSurvey ensemble probability

# set color palette
pal <- colorBin("Reds", domain = 0:1) 

# render map
w <- leaflet() %>% 
  addProviderTiles(providers$OpenStreetMap.Mapnik) %>%
  addRasterImage(pred, colors = pal, opacity = 0.5) %>%
  addLegend(pal = pal, values = values(pred), title = "BP")
w ## plot widget 

# save widget
saveWidget(w, 'MW_BP.html', selfcontained = T)
