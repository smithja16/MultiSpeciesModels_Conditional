####################################################################
###     ~  Fitting and Predicting Joint and Stacked SDMs   ~     ###
###  - An evaluation of conditional prediction, particularly...  ###
###      predicting fishery discards using retained catch.       ###
###  - James A. Smith Feb 2024                                   ###
####################################################################

## This script contains a function to fit and evaluate a hurdle
## stacked random forest species model


fit_rf_hurdle_cond <- function(spp_list_response,  #the discard species are the response vraiables
                               spp_list_covar,  #the retained species are the covariates
                               first_taxon,  #first taxa in species list
                               data_train_pres,  #presence data (1/0) for fitting/training
                               data_train,  #biomass data for fitting/training
                               data_test_pres,  #NA or presence data (1/0) for testing (cross validation)
                               data_test,  #NA or biomass data for testing (cross validation)
                               pres_prevalence,  #T or F; T means creating zeros from < threshold predicted probabilities
                               conditional, #T fits models for conditional pred., F fits models for marginal pred.
                               save_models) {  #T or F; saves both hurdle components for every taxa
  
  save_rf_perf <- data.frame(species=spp_list_response,
                             pres_oob_error_perc=0,
                             abund_perc_var_expl=0)
  
  colx <- which(names(data_train_pres)==first_taxon)
  
  fitted_values <- data_train[,colx:ncol(data_train)]
  fitted_values[] <- 0
  test_cv_perf <- data.frame(species=spp_list_response,
                             auc_pres=NA,
                             r2_total=NA,
                             rmae_total=NA,
                             rmse_total=NA)
  
  if (is.na(data_test_pres)[1] == T) {
    save_test_pres <- 0
    save_test_total <- 0
  } else {
    save_test_pres <- as.data.frame(matrix(data=NA,
                                           nrow=nrow(data_test_pres),
                                           ncol=length(spp_list_response)))
    save_test_total <- save_test_pres
  }
  
  save_prev <- data.frame(species=spp_list, prevalence=0)  #if using 'prevalence' option
  
  pb <- txtProgressBar(min=0, max=length(spp_list_response))
  
  for (ss in 1:length(spp_list_response)) {
    
    setTxtProgressBar(pb,ss)
    
    sppx <- spp_list_response[ss]
    data_train_pres[,sppx] <- as.factor(data_train_pres[,sppx])
    
    ## hurdle model
    if (conditional == T) {
      Mrf_pres <- randomForest(as.formula(paste0(sppx," ~ meanLat + meanDepth_ftm + month + meanGlorys_sst +
                                         meanGlorys_mld + lunar_illum + area_swept_km2 + total_catch_kg + ",
                                                 paste(spp_list_covar, collapse=" + "))),
                               data=data_train_pres, importance=T,
                               keep.inbag=T, ntree=1201, mtry=3,
                               keep.forest=T)
    } else {
      Mrf_pres <- randomForest(as.formula(paste0(sppx," ~ meanLat + meanDepth_ftm + month + meanGlorys_sst +
                                         meanGlorys_mld + lunar_illum + area_swept_km2 + total_catch_kg")),
                               data=data_train_pres, importance=T,
                               keep.inbag=T, ntree=1201, mtry=2,
                               keep.forest=T)
    }
    
    obs_dat_trip_split_abundx <- data_train  #remove all rows that are zero for RESPONSE SPECIES ONLY
    rows0 <- which(obs_dat_trip_split_abundx[,sppx]==0)  #which rows are zero for this species
    if (length(rows0) > 0) {
      obs_dat_trip_split_abundx <- obs_dat_trip_split_abundx[-rows0,]
    }
    
    if (conditional == T) {
      Mrf_abund <- randomForest(as.formula(paste0(sppx," ~ meanLat + meanDepth_ftm + month + meanGlorys_sst +
                                         meanGlorys_mld + lunar_illum + area_swept_km2 + total_catch_kg + ",
                                                  paste(spp_list_covar, collapse=" + "))),
                                data=obs_dat_trip_split_abundx, importance=T,
                                keep.inbag=T, ntree=1201, mtry=3,
                                keep.forest=T)
    } else {
      Mrf_abund <- randomForest(as.formula(paste0(sppx," ~ meanLat + meanDepth_ftm + month + meanGlorys_sst +
                                         meanGlorys_mld + lunar_illum + area_swept_km2 + total_catch_kg")),
                                data=obs_dat_trip_split_abundx, importance=T,
                                keep.inbag=T, ntree=1201, mtry=2,
                                keep.forest=T)
    }
    
    if (save_models == T) {
      if (conditional == T) {
        saveRDS(Mrf_pres, paste0("Mrf_pres_cond_spp_",sppx,".rds"))
        saveRDS(Mrf_abund, paste0("Mrf_abund_cond_spp_",sppx,".rds"))
      } else {
        saveRDS(Mrf_pres, paste0("Mrf_pres_marg_spp_",sppx,".rds"))
        saveRDS(Mrf_abund, paste0("Mrf_abund_marg_spp_",sppx,".rds"))
      }
    }
    
    save_rf_perf$pres_oob_error_perc[ss] <- round((Mrf_pres$err.rate[nrow(Mrf_pres$err.rate),1]*100),2)
    save_rf_perf$abund_perc_var_expl[ss] <- round(Mrf_abund$rsq[nrow(Mrf_pres$err.rate)]*100,2)
    
    
    if (pres_prevalence == T) {  #calculate threshold for this species based on its prevalence
      colxx <- which(names(data_train_pres) == sppx)
      pres_thresholdx <- data_train_pres[,colxx]
      pres_threshold <- length(pres_thresholdx[pres_thresholdx==1])/nrow(data_train)
      if (pres_threshold > 0.95) { pres_threshold <- 0.95 }
      save_prev$prevalence[ss] <- pres_threshold
    }
    
    #calculate fitted values
    P_pres <- predict(Mrf_pres, newdata=data_train, type="prob")[,2]  #probability of presence
    P_pres[P_pres < pres_threshold] <- 0
    P_abund <- predict(Mrf_abund, newdata=data_train, type="response")
    P_total <- P_pres * P_abund
    
    fitted_values[,sppx] <- P_total
    
    #Calculate predictions on test data
    if (is.na(data_test_pres)[1] == T) { 
    } else {
      P_pres_test <- predict(Mrf_pres, data_test, type="prob")[,2]
      P_pres_test[P_pres_test < pres_threshold] <- 0
      P_abund_test <- predict(Mrf_abund, data_test, type="response")
      P_total_test <- P_pres_test*P_abund_test
      
      save_test_pres[,ss] <- P_pres_test
      save_test_total[,ss] <- P_total_test
      
      r2_total <- R2(P_total_test, data_test[,sppx])  #predicted then observed
      
      if (sum(data_test_pres[,sppx]) == 0 | mean(data_test_pres[,sppx]) == 1) {
        auc_pres <- NA  #AUC can't be calculated if test data are all zeros or all ones
      } else {
        auc_pres <- suppressMessages(auc(data_test_pres[,sppx],P_pres_test))
      }
      rmae_total <- RMAE(P_total_test, data_test[,sppx])  #predicted then observed
      rmse_total <- RMSE(P_total_test, data_test[,sppx])
      
      test_cv_perf$auc_pres[ss] <- auc_pres
      test_cv_perf$r2_total[ss] <- r2_total
      test_cv_perf$rmae_total[ss] <- rmae_total
      test_cv_perf$rmse_total[ss] <- rmse_total
      
    }
    
  }
  
  return(list(fitted_totals=fitted_values,
              save_rf_perf=save_rf_perf,
              test_pred_pres=save_test_pres,
              test_pred_tot=save_test_total,
              test_cv_perf=test_cv_perf))
}
