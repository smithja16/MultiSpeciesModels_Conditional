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

cv_type <- "conditional"  #specify which predictions to test: "conditional" or "marginal" (this just saves me writing almost identical code twice)

## Set up folds
n_folds <- 5; n_repeats <- 3
my_folds <- CreateFolds(data=obs_dat_trip_split,
                        n_folds=n_folds,
                        n_repeats=n_repeats,
                        seed=111)
folds <- my_folds$folds
fold_rows <- my_folds$fold_rows

## Objects for saving results
save_auc <- folds
dfspp <- as.data.frame(matrix(data=0, nrow=nrow(save_auc), ncol=length(spp_list_dis)))
colnames(dfspp) <- spp_list_dis
save_auc <- cbind(save_auc, dfspp)
save_rmae <- save_auc
save_rmse <- save_auc
save_model_perf_oob_err <- save_auc
save_model_perf_var_expl <- save_auc
save_tot_dis <- folds
save_tot_dis$r2_tot_dis <- 0
save_tot_dis$rmae_tot_dis <- 0
save_tot_dis$sum_tot_dis_obs <- 0
save_tot_dis$sum_tot_dis_pred <- 0
save_all_test_preds <- list()

for (nn in 1:(n_folds*n_repeats)) {
  
  print(paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn]))
  
  data_train_presx <- obs_dat_trip_split_pres[-fold_rows[[nn]],]
  data_trainx <- obs_dat_trip_split[-fold_rows[[nn]],]
  data_test_presx <- obs_dat_trip_split_pres[fold_rows[[nn]],]
  data_testx <- obs_dat_trip_split[fold_rows[[nn]],]
  zeros <- colSums(data_testx[,colx:ncol(data_testx)])
  print(paste0("Taxa with zeros in test data = ", length(zeros[zeros==0])))
  
  rf_cond_cv <- fit_rf_hurdle_cond(spp_list_response=spp_list_dis,
                                   spp_list_covar=spp_list_ret,
                                   data_train_pres=data_train_presx,
                                   data_train=data_trainx,
                                   data_test_pres=data_test_presx,
                                   data_test=data_testx,
                                   pres_threshold=0,
                                   pres_prevalence=T,  #thresholding by prevalence
                                   conditional=F,  #T or F; T adds other spp as covariates; use F to estimate contribution of other spp to performance (an 'enviro only' scenario)
                                   save_models=F)  #don't save models during CV
  
  save_all_test_preds[[nn]] <- rf_cond_cv$test_pred_tot  #predicted discards for all discard taxa, test set (use this for comparison with simulated discard series)
  tot_disc_obs_pred <- rowSums(rf_cond_cv$test_pred_tot)  #predicted total discard biomass per trip in test set
  tot_disc_test_pred <- sum(tot_disc_obs_pred)  #predicted total discard biomass in test set (~84 observations)
  colcut <- which(names(data_testx) %in% spp_list_dis)  #only sum the predicted discarded taxa
  tot_disc_obs_obs <- rowSums(data_testx[,colcut])  #observed total discard biomass per trip in test set
  tot_disc_test_obs <- sum(tot_disc_obs_obs)  #observed total discard biomass in test set
  plot(tot_disc_obs_pred, tot_disc_obs_obs,
       main=paste0("Rep=",folds$rep[nn], ", Fold=",folds$fold[nn])); abline(a=0,b=1)  #generally pretty good
  to_save <- cbind(tot_disc_obs_pred,tot_disc_obs_obs)  #save all the obs vs pred total discard biomass for trips
  if (nn == 1) { save_tot_dis_all <- to_save }
  if (nn > 1) { save_tot_dis_all <- rbind(save_tot_dis_all, to_save)}
  
  save_tot_dis$r2_tot_dis[nn] <- R2(tot_disc_obs_pred, tot_disc_obs_obs)
  save_tot_dis$rmae_tot_dis[nn] <- RMAE(tot_disc_obs_pred, tot_disc_obs_obs)
  save_tot_dis$sum_tot_dis_obs[nn] <- tot_disc_test_obs
  save_tot_dis$sum_tot_dis_pred[nn] <- tot_disc_test_pred
  # ^ compare these to intercept only model
  
  save_auc[nn, 3:ncol(save_auc)] <- rf_cond_cv$test_cv_perf$auc_pres
  save_rmae[nn, 3:ncol(save_auc)] <- rf_cond_cv$test_cv_perf$rmae_total
  save_rmse[nn, 3:ncol(save_auc)] <- rf_cond_cv$test_cv_perf$rmse_total
  save_model_perf_oob_err[nn, 3:ncol(save_auc)] <- rf_cond_cv$save_rf_perf$pres_oob_error_perc
  save_model_perf_var_expl[nn, 3:ncol(save_auc)] <- rf_cond_cv$save_rf_perf$abund_perc_var_expl
  
}

saveRDS(list(rfH_cond_cv_auc_pres=save_auc,
             rfH_cond_cv_rmae_total=save_rmae,
             rfH_cond_cv_rmse_total=save_rmse,
             rfH_cond_cv_modelPerf_oobErr=save_model_perf_oob_err,
             rfH_cond_cv_modelPerf_varExpl=save_model_perf_var_expl,
             rfH_cond_cv_tot_dis=save_tot_dis,
             rfH_cond_cv_tot_dis_all=save_tot_dis_all,
             rfH_cond_cv_all_dis_test_preds=save_all_test_preds),
        paste0("C:/Users/smithj08/OneDrive - DPIE/3 Multispecies models and discard prediction/",
               "R code/RF_Spp46_CV_results_enviroOnly_presThreshold_noBoat_noLocationVars_280224.rds"))  #RF_Spp46_CV_results_presThreshold_210224.rds
