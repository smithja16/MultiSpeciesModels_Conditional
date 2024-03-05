####################################################################
###     ~  Fitting and Predicting Joint and Stacked SDMs   ~     ###
###  - An evaluation of conditional prediction, particularly...  ###
###      predicting fishery discards using retained catch.       ###
###  - James A. Smith Feb 2024                                   ###
####################################################################

## This script contains code to fit, evaluate, and predict with a hurdle
## stacked single species model using random forests
## This code fits two models, one without other taxa as covariates (for marginal
## prediction), and one with other taxa as covariates (for conditional prediction). 
## We use k-folds cross-validation (CV) to measure predictive performance of
## these two types of predictions.


library(randomForest)
library(pROC)

source("RandomForest fitting function.R")
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

## Create presence hurdle component
biomass_data_pres <- biomass_data
biomass_data_pres[biomass_data_pres > 0] <- 1  #presence

all_data_pres <- cbind(enviro_data, biomass_data_pres)
all_data <- cbind(enviro_data, biomass_data)  #data is 

## Species lists
spp_list_rf <- names(biomass_data)
spp_list_ret <- readRDS("spp_list_ret.rds")  #retained taxa
spp_list_dis <- readRDS("spp_list_dis.rds")  #discarded taxa


############################
## Fit full hurdle models ##
############################

rf_marginal <- fit_rf_hurdle_cond(spp_list_response = spp_list_dis,  
                                  spp_list_covar = spp_list_ret,
                                  first_taxon = spp_list_rf[1],
                                  data_train_pres = all_data_pres,
                                  data_train = all_data,
                                  data_test_pres = NA,  #no test data for full model
                                  data_test = NA,
                                  pres_prevalence = T,
                                  conditional = F,  #for marginal = F
                                  save_models = T)  #save individual models for evaluation

rf_conditional <- fit_rf_hurdle_cond(spp_list_response = spp_list_dis,  
                                     spp_list_covar = spp_list_ret,
                                     first_taxon = spp_list_rf[1],
                                     data_train_pres = all_data_pres,
                                     data_train = all_data,
                                     data_test_pres = NA,  #no test data for full model
                                     data_test = NA,
                                     pres_prevalence = T,
                                     conditional = T,  #for conditional = T
                                     save_models = T)  #save individual models for evaluation

saveRDS(rf_marginal, paste0("Fit_rf_marg_summary_results.rds"))
saveRDS(rf_conditional, paste0("Fit_rf_cond_summary_results.rds"))


##########################
## Evaluate full models ##
##########################

## OOB performance of full models
quantile(rf_marginal$save_rf_perf$pres_oob_error_perc)
quantile(rf_marginal$save_rf_perf$posi_perc_var_expl)
quantile(rf_conditional$save_rf_perf$pres_oob_error_perc)
quantile(rf_conditional$save_rf_perf$posi_perc_var_expl)

## An example model
sppx_plot <- "Platycephalus.caeruleopunctatus_dis"
rf_presx <- readRDS(paste0("Mrf_pres_cond_spp_",sppx_plot,".rds"))

## Plot response to some covariates
partialPlot(rf_presx, pred.data=all_data,
            x.var="Platycephalus.caeruleopunctatus_ret",
            which.class="1")
partialPlot(rf_presx, pred.data=all_data,
            x.var="meanDepth_ftm",
            which.class="1")

## Variable importance
importance(rf_presx)


###################################
##  Cross-Validation - k-folds   ##
###################################

## Do "marginal" and "conditional"

cv_type <- "conditional"  #specify predictions to test: "conditional" or "marginal" (this just saves me writing almost identical code twice)

## Set up folds
n_folds <- 5; n_repeats <- 3
my_folds <- CreateFolds(data = all_data,
                        n_folds = n_folds,
                        n_repeats = n_repeats,
                        seed = 111)  #same seed gives same folds
folds <- my_folds$folds
fold_rows <- my_folds$fold_rows

## Objects for saving results
save_auc <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_auc), ncol=length(spp_list_dis)))
colnames(dfspp) <- spp_list_dis
save_auc <- cbind(save_auc, dfspp)
save_rmae <- save_auc
save_rmse <- save_auc
save_all_test_preds <- list()

## Cross validation loop
for (nn in 1:(n_folds*n_repeats)) {
  
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  data_train_presx <- all_data_pres[-fold_rows[[nn]],]
  data_trainx <- all_data[-fold_rows[[nn]],]
  data_test_presx <- all_data_pres[fold_rows[[nn]],]
  data_testx <- all_data[fold_rows[[nn]],]
  
  if (cv_type == "conditional") { conditional = T } else { conditional = F }
  
  rf_cond_cv <- fit_rf_hurdle_cond(spp_list_response = spp_list_dis,
                                   spp_list_covar = spp_list_ret,
                                   first_taxon = spp_list_rf[1],
                                   data_train_pres = data_train_presx,
                                   data_train = data_trainx,
                                   data_test_pres = data_test_presx,
                                   data_test = data_testx,
                                   pres_prevalence = T,
                                   conditional = conditional,
                                   save_models = F)  #no reason to save models during CV
  
  save_all_test_preds[[nn]] <- rf_cond_cv$test_pred_tot  #predicted discards for all discard taxa, test set
  save_auc[nn, 3:ncol(save_auc)] <- rf_cond_cv$test_cv_perf$auc_pres
  save_rmae[nn, 3:ncol(save_auc)] <- rf_cond_cv$test_cv_perf$rmae_total
  save_rmse[nn, 3:ncol(save_auc)] <- rf_cond_cv$test_cv_perf$rmse_total
  
}

saveRDS(list(rf_cv_auc_pres = save_auc,
             rf_cv_rmae_total = save_rmae,
             rf_cv_rmse_total = save_rmse,
             rf_cv_discard_test_preds = save_all_test_preds),
        paste0("RF_CV_results_",cv_type,".rds"))

AUCx <- apply(save_auc[3:ncol(save_auc)], 2, FUN=median, na.rm=T)
hist(AUCx, main=paste0("median AUC Cond=",conditional," mean=",round(mean(AUCx), 2)))

