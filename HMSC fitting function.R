####################################################################
###     ~  Fitting and Predicting Joint and Stacked SDMs   ~     ###
###  - An evaluation of conditional prediction, particularly...  ###
###      predicting fishery discards using retained catch.       ###
###  - James A. Smith Feb 2024                                   ###
####################################################################

## This script contains a function to fit and evaluate a hurdle
## joint species model using Hmsc
## This uses spatial random effects for the positive hurdle component only
## due to issues fitting the abundance components in Hmsc


fit_hmsc_hurdle <- function(train_data_pres,  #presence data (1/0) for fitting/training
                            train_data_posi, #positive part for fitting/training (put NAs instead of zeros)
                            species_list,  #list of all discard and retained taxa
                            first_taxon,  #first taxa in species list
                            nChains,  #number of chains for MCMC
                            thin,  #thinning rate for MCMC
                            samples,  #number of final samples for MCMC
                            transient,  #burn-in for MCMC
                            nParallel,  #number of cores for model fitting
                            YScale,  #T or F; for scaling reponse of normal distribution (T recommended)
                            calc_fit,  #T or F; multiply hurdle components for total biomass?
                            save_model,  #T or F; save models?
                            save_model_suffix) {  #suffix for saved model objects
  #total number of iterations = thin * samples + transient
  
  # Save data
  model_fit <- data.frame(species=species_list, AUC_pres=0, R2_abund=0)
  
  # Set up models
  colx <- which(names(train_data_pres)==first_taxon)
  
  X_pres = train_data_pres[,1:(colx-1)]
  X_abund = train_data_posi[,1:(colx-1)]
  
  XFormula = ~ poly(meanLat, degree=2, raw=T) +
    poly(meanDepth_ftm, degree=2, raw=T) +
    poly(meanGlorys_sst, degree=2, raw=T) +
    meanGlorys_mld +
    lunar_illum +
    area_swept_km2 +
    total_catch_kg +
    month
  
  XFormula_abund = ~ poly(meanLat, degree=2, raw=T) +
    poly(meanDepth_ftm, degree=2, raw=T) +
    poly(meanGlorys_sst, degree=2, raw=T) +
    area_swept_km2 +
    total_catch_kg
  
  
  Y_pres <- as.matrix(train_data_pres[,colx:ncol(train_data_pres)])
  Y_abund <- as.matrix(train_data_posi[,colx:ncol(train_data_posi)])
  Y_abund_log <- log(Y_abund)
  
  xycoords <- train_data_pres[,c("meanLong", "meanLat")]
  colnames(xycoords) <- c("x-coordinate","y-coordinate")
  rownames(xycoords) <- 1:nrow(xycoords)
  while (anyDuplicated(xycoords) > 0) {  #coords must be unique for Hmsc
    dups <- anyDuplicated(xycoords)
    xycoords$`x-coordinate`[dups[1]] <- xycoords$`x-coordinate`[dups[1]] + 0.002  #slightly alter location if a duplicate
  }
  
  studyDesign = data.frame(sample=as.factor(1:nrow(X_pres)),
                           vessel=droplevels(train_data_pres$boatID))  #droplevels essential for data subsets
  
  Knots <- constructKnots(sData = xycoords, knotDist = 0.15, minKnotDist = 0.15)
  plot(xycoords[,1],xycoords[,2],pch=18, asp=1)
  points(Knots[,1],Knots[,2],col='red',pch=18)
  rL_spat <- HmscRandomLevel(sData = xycoords, sMethod = 'GPP', sKnot = Knots)  #spatial for presence part
  rL_v <- HmscRandomLevel(units=studyDesign$vessel)
  rL_s <- HmscRandomLevel(units=studyDesign$sample)  #non-spatial for positive part
  
  Mhm_pres <- Hmsc(Y=Y_pres,
                   XData=X_pres,
                   XFormula=XFormula,
                   studyDesign=studyDesign,
                   ranLevels=list(sample=rL_spat, vessel=rL_v),
                   distr="probit")
  
  Mhm_abund <- Hmsc(Y=Y_abund_log,
                    XData=X_abund,
                    XFormula=XFormula_abund,
                    studyDesign=studyDesign,
                    ranLevels=list(sample=rL_s, vessel=rL_v),
                    distr="normal",
                    YScale=YScale)
  # *^ note the logged response variable
  
  # Fit models
  print("Fitting Mpres spatial model")
  t1 <- Sys.time(); t1
  Mhm_pres <- sampleMcmc(Mhm_pres,
                         thin=thin,
                         samples=samples,
                         transient=transient,
                         nChains=nChains,
                         nParallel=nParallel, 
                         verbose=500)
  t2 <- Sys.time()
  print(t2-t1)
  
  print("Fitting Mposi model")
  t1 <- Sys.time(); t1
  Mhm_abund <- sampleMcmc(Mhm_abund,
                          thin=thin,
                          samples=samples,
                          transient=transient,
                          nChains=nChains,
                          nParallel=nParallel,
                          verbose=500)
  t2 <- Sys.time()
  print(t2-t1)
  
  if (save_model == T) {
    saveRDS(Mhm_pres, paste0("Mpres_",save_model_suffix,".rds"))
    saveRDS(Mhm_abund, paste0("Mposi_",save_model_suffix,".rds"))
  }
  
  if (calc_fit == T) {
    print("Calculating fits")
    preds_p <- computePredictedValues(Mhm_pres, expected=F)
    preds_a <- computePredictedValues(Mhm_abund, expected=T)
    
    Mhm_pres_fit <- evaluateModelFit(hM=Mhm_pres, predY=preds_p)
    Mhm_abund_fit <- evaluateModelFit(hM=Mhm_abund, predY=preds_a)

    model_fit$AUC_pres <- Mhm_pres_fit$AUC
    model_fit$R2_abund <- Mhm_abund_fit$R2
    
    # Calculate total biomass
    # probability of presence * abundance when present
    preds_t <- preds_a
    for (ff in 1:dim(preds_a)[3]) {
      preds_t[,,ff] <- preds_p[,,ff] * exp(preds_a[,,ff])  #undo the log
    }
    preds_t_mean <- apply(preds_t, c(1,2), mean)  #use mean then remove outliers if necessary
    
    saveRDS(Mhm_pres_fit, paste0("Mpres_fit_",save_model_suffix,".rds"))
    saveRDS(Mhm_abund_fit, paste0("Mposi_fit_",save_model_suffix,".rds"))
  }
  
  if (calc_fit == T) {
    return(list(model_fit=model_fit,
                fitted_totals=preds_t_mean))
  } else {
    return(list(model_fit=model_fit))
  }
  
}
