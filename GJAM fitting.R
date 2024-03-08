####################################################################
###     ~  Fitting and Predicting Joint and Stacked SDMs   ~     ###
###  - An evaluation of conditional prediction, particularly...  ###
###      predicting fishery discards using retained catch.       ###
###  - James A. Smith Feb 2024                                   ###
####################################################################

## This script contains code to fit, evaluate, and predict with two
## joint species models using GJAM
## The first model is a presence-absence model, to evaluate performance using AUC
## The second is a continuous abundance models, to evaluate performance using RMSE

## Using 'gjam' version 2.6.2, for prediction I need to construct the model matrix 
## manually when using quadratics and a month factor. This may change as the 
## package is updated.

## Also, this code increases the example data to avoid dimension reduction by GJAM, 
## do not do this normally. Dimension reduction must be avoided for conditional
## prediction: so ensure there are numerous times more observations than species.

## Don't include in presence-absence model any taxa that are always present


library(gjam)
library(pROC)
library(corrplot)
library(gclus)

source("Performance functions.R")
source("Cross validation function.R")


########################
## Load catch data    ##
########################
## These catch and environment data are fake, used to illustrate the format
## The original observer data are confidential

## Load environmental variables and biomass (catch) data
enviro_data <- readRDS("example_enviro_data.rds")  #mean conditions per trip; ensure factors have been labelled as factors
biomass_data <- readRDS("example_catch_data.rds")  #catch data (kg per trip)

  enviro_data <- rbind(enviro_data, enviro_data)  #*** doubles example data to avoid dimension reduction - for example data only!
  biomass_data <- rbind(biomass_data, biomass_data)  #*** doubles example data to avoid dimension reduction - for example data only!

## Create presence hurdle component
biomass_data_pres <- biomass_data
biomass_data_pres[biomass_data_pres > 0] <- 1  #presence

## Species lists
spp_list_gjam <- names(biomass_data)
spp_list_ret <- readRDS("spp_list_ret.rds")  #retained taxa
spp_list_dis <- readRDS("spp_list_dis.rds")  #discarded taxa

## Prep
names(enviro_data) <- c("meanLat","meanLong","meanDepthftm","meanSST",
                        "meanMLD","lunar","areaSweptkm2","totalCatchkg",
                        "month","boatID") 
# ^ gjam does not allow underscores
ydata <- biomass_data[,c(spp_list_ret, spp_list_dis)]  #reorder taxa
#ydata <- ydata/10  # **this can help gjam convergence, but if used need to x10 the predictions 
ydata_pa <- ydata
ydata_pa[ydata_pa > 0] <- 1

## Model formula
formula <- as.formula( ~ meanLat + I(meanLat^2) +
                         meanDepthftm + I(meanDepthftm^2) +
                         meanSST + I(meanSST^2) +
                         areaSweptkm2 + meanMLD + 
                         lunar + month +
                         totalCatchkg )


###########################
## Fit full models       ##
###########################

## This model is fit using all the data

ml_pa  <- list(ng = 10000,  #gibbs steps
               burnin = 5000,  #burnin
               typeNames = "PA",    #presence-absence
               random = "boatID")

ml_ca  <- list(ng = 10000,  #gibbs steps
               burnin = 5000,  #burnin
               typeNames = "CA",    #continuous abundance
               random = "boatID")

myout_pa <- gjam(formula = formula,
                xdata = enviro_data,
                ydata = ydata_pa,
                modelList = ml_pa)

myout_ca <- gjam(formula = formula,
                 xdata = enviro_data,
                 ydata = ydata,
                 modelList = ml_ca)

summary(myout_pa)
summary(myout_ca)
saveRDS(myout_pa,"myout_pa_full.rds")
saveRDS(myout_ca,"myout_ca_full.rds")


###########################
## Evaluate full models  ##
###########################

## all result plots, saves in new plot folder
gp_pa <- gjamPlot(myout_pa, plotPars = list(GRIDPLOTS=T, SAVEPLOTS = T))

# this will overwrite the previous plots
gp_ca <- gjamPlot(myout_ca, plotPars = list(GRIDPLOTS=T, SAVEPLOTS = T))

## covariate importance for full response matrix
gjamSensitivity(myout_pa, nsim=100)
gjamSensitivity(myout_ca, nsim=100)

## residual correlations
pc <- myout_pa$parameters$corMu
corrplot(pc[order.single(pc), order.single(pc)],
         method = "color", type="lower", col = colorRampPalette(c("blue","white","red"))(200),
         title = paste0("GJAM presence, max=",round(max(pc[pc<1]),2),", min=",round(min(pc[pc<1]),2)),
         mar=c(0,0,1,0), tl.cex=0.6, tl.col="black")

pc <- myout_ca$parameters$corMu
corrplot(pc[order.single(pc), order.single(pc)],
         method = "color", type="lower", col = colorRampPalette(c("blue","white","red"))(200),
         title = paste0("GJAM presence, max=",round(max(pc[pc<1]),2),", min=",round(min(pc[pc<1]),2)),
         mar=c(0,0,1,0), tl.cex=0.6, tl.col="black")


###################################
##  Cross-Validation - k-folds   ##
###################################

## Set up folds
n_folds <- 5; n_repeats <- 3
my_folds <- CreateFolds(data = all_data,
                        n_folds = n_folds,
                        n_repeats = n_repeats,
                        seed = 111)  #same seed gives same folds
folds <- my_folds$folds
fold_rows <- my_folds$fold_rows

## Prep for saving
save_auc <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_auc),
                              ncol=length(spp_list_dis)))
colnames(dfspp) <- spp_list_dis
save_aucC <- cbind(save_auc, dfspp)
save_aucUn <- save_aucC

### Cross-validation for presence-absence model

## Prep taxa groups
cond_spp <- 1:22  #columns of ydata_pa to condition on (retained taxa)
pred_spp <- 23:46  #columns of ydata_pa to predict (discarded taxa)
num_dis_spp <- length(pred_spp)
ssx <- max(cond_spp)

## Cross-validation loop
for (nn in 1:(n_folds*n_repeats)) {
  
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  data_trainx <- enviro_data[-fold_rows[[nn]],]
  data_testx <- enviro_data[fold_rows[[nn]],]
  data_trainy <- ydata_pa[-fold_rows[[nn]],]
  data_testy <- ydata_pa[fold_rows[[nn]],]

  quads <- data.frame(lat = (data_testx$meanLat)^2,
                      depth = (data_testx$meanDepthftm)^2,
                      sst = (data_testx$meanSST)^2)
  names(quads) <- c("I(meanLat^2)", "I(meanDepthftm^2)","I(meanSST^2)")
  data_testx <- cbind(data_testx, quads)
  # ^ add the quadratics manually
  monthx <- as.data.frame(model.matrix(formula, data=data_testx))
  data_testx <- cbind(data_testx, monthx[,11:21])
  # ^ add months in model matrix format
  
  myoutx <- gjam(formula = formula,
                 xdata = data_trainx,
                 ydata = data_trainy,
                 modelList = ml_pa)
  
  newUn <- list(xdata = data_testx,  #unconditional
                nsim = 500)
  newC <- list(xdata = data_testx,
               ydataCond = data_testy[,cond_spp],  #conditional
               nsim = 1000)
  
  Pun <- gjamPredict(output = myoutx, newdata = newUn)$sdList$yMu[,pred_spp]  #unconditional predictions
  Pc <- gjamPredict(output = myoutx, newdata = newC)$sdList$yMu[,pred_spp]  #conditional predictions
  
  for (sppx in 1:num_dis_spp)  {
    save_aucC[nn,2+sppx] <- round(suppressMessages(auc(data_testy[,ssx+sppx],Pc[,sppx])),2)
    save_aucUn[nn,2+sppx] <- round(suppressMessages(auc(data_testy[,ssx+sppx],Pun[,sppx])),2)
  }
  
}

gjam_cv_pa_summary <- list(save_auc_unconditional = save_aucUn,
                             save_auc_conditional = save_aucC)
saveRDS(gjam_cv_pa_summary,"gjam_cv_pa_summary.rds")



### Cross-validation for continuous abundance model

save_r2 <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_r2),
                              ncol=length(spp_list_dis)))
colnames(dfspp) <- spp_list_dis
save_r2c <- cbind(save_r2, dfspp)
save_r2un <- save_r2c
save_rmsec <- save_r2c
save_rmseun <- save_r2c


save_all_test_predsC <- list()
save_all_test_predsUn <- list()

## Prep taxa groups
cond_spp <- 1:22  #columns of ydata to condition on (retained taxa)
pred_spp <- 23:46  #columns of ydata to predict (discarded taxa)
num_dis_spp <- length(pred_spp)
ssx <- max(cond_spp)

## Cross-validation loop
for (nn in 1:(n_folds*n_repeats)) {
  
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  data_trainx <- enviro_data[-fold_rows[[nn]],]
  data_testx <- enviro_data[fold_rows[[nn]],]
  data_trainy <- ydata[-fold_rows[[nn]],]
  data_testy <- ydata[fold_rows[[nn]],]
  
  quads <- data.frame(lat = (data_testx$meanLat)^2,
                      depth = (data_testx$meanDepthftm)^2,
                      sst = (data_testx$meanSST)^2)
  names(quads) <- c("I(meanLat^2)", "I(meanDepthftm^2)","I(meanSST^2)")
  data_testx <- cbind(data_testx, quads)
  # ^ add the quadratics manually
  monthx <- as.data.frame(model.matrix(formula, data=data_testx))
  data_testx <- cbind(data_testx, monthx[,11:21])
  # ^ add months in model matrix format
  
  myoutx <- gjam(formula = formula,
                 xdata = data_trainx,
                 ydata = data_trainy,
                 modelList = ml_ca)
  
  newUn <- list(xdata = data_testx,  #unconditional
                nsim = 500)
  newC <- list(xdata = data_testx,
               ydataCond = data_testy[,cond_spp],  #conditional
               nsim = 1000)
  
  Pun <- gjamPredict(output = myoutx, newdata = newUn)$sdList$yMu[,pred_spp]  #unconditional predictions
  Pc <- gjamPredict(output = myoutx, newdata = newC)$sdList$yMu[,pred_spp]  #conditional predictions
  
  save_all_test_predsUn[[nn]] <- Pun # save predictions
  save_all_test_predsC[[nn]] <- Pc # save predictions
  
  for (sppx in 1:num_dis_spp)  {
    save_r2un[nn,2+sppx] <- round(R2(Pun[,sppx], data_testy[,ssx+sppx]),2)
    save_r2c[nn,2+sppx] <- round(R2(Pc[,sppx], data_testy[,ssx+sppx]),2)
    save_rmseun[nn,2+sppx] <- round(RMSE(Pun[,sppx], data_testy[,ssx+sppx]),2)
    save_rmsec[nn,2+sppx] <- round(RMSE(Pc[,sppx], data_testy[,ssx+sppx]),2)
  }
  
}

gjam_cv_ca_summary <- list(save_r2_uncond = save_r2un,
                           save_r2_cond = save_r2c,
                           save_rmse_uncond = save_rmseun,
                           save_rmse_cond = save_rmsec,
                           save_cond_test_preds = save_all_test_predsC,
                           save_uncond_test_preds = save_all_test_predsUn)

saveRDS(gjam_cv_ca_summary,"gjam_cv_ca_summary.rds")


