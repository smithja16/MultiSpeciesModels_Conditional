####################################################################
###     ~  Fitting and Predicting Joint and Stacked SDMs   ~     ###
###  - An evaluation of conditional prediction, particularly...  ###
###      predicting fishery discards using retained catch.       ###
###  - James A. Smith Feb 2024                                   ###
####################################################################

## This script contains code to fit, evaluate, and predict with a hurdle
## joint species model using Hmsc
## This code fits a full model, then does two stages of cross-validation (CV)
## the first stage is joint prediction, and the second stage is joint
## conditional prediction (i.e. assumes some species are known)
## By comparing the CV predictive performance of joint and joint conditional
## predictions, we can see whether there is value is using logbooks of
## retained species to help predict unobserved discards (i.e. conditional prediction)


library(Hmsc)
library(pROC)
library(corrplot)
library(gclus)

source("HMSC fitting function.R")


########################
## Load catch data    ##
########################
## These fishing trip-level data are fake, used to illustrate the format
## the original observer data are confidential

## Load environmental variables and biomass (catch) data
enviro_data <- readRDS("example_enviro_data.rds")  #mean conditions per trip; ensure factors have been labelled as factors; fake data
biomass_data <- readRDS("example_catch_data.rds")  #fake catch data (kg per trip)

## Create presence and positive hurdle components
biomass_data_pres <- biomass_data
biomass_data_pres[biomass_data_pres > 0] <- 1  #presence
biomass_data_posi <- biomass_data
biomass_data_posi[biomass_data_posi == 0] <- NA  #positive

all_data_pres <- cbind(enviro_data, biomass_data_pres)
all_data_posi <- cbind(enviro_data, biomass_data_posi)

## Species lists
spp_list_hmsc <- names(biomass_data)
spp_list_ret <- readRDS("spp_list_ret.rds")  #retained taxa
spp_list_dis <- readRDS("spp_list_dis.rds")  #discarded taxa

## Maximum biomass of each taxa (for later truncation)
spp_max <- apply(biomass_data, 2, FUN=max)


###########################
## Fit full hurdle model ##
###########################

save_model_suffix = "mymodel_date"

fit_hmsc <- fit_hmsc_hurdle(train_data_pres = obs_data_pres,
                            train_data_posi = obs_data_posi,
                            species_list = spp_list_hmsc,
                            first_taxon = spp_list_hmsc[1],
                            nChains = 3,
                            thin = 10,
                            samples = 1000,
                            transient = 10000,
                            nParallel = 1,  #recommend start with 1 to benchmark computation time
                            YScale = T,
                            calc_fit = T,
                            save_model = T,
                            save_model_suffix = save_model_suffix)

saveRDS(fit_hmsc, paste0("Fit_hmsc_summary_results.rds"))


#########################
## Evaluate full model ##
#########################

## Summary goodness-of-fits
fit_hmsc$model_fit  #any taxa always present will have an AUC = NA
mean(fit_hmsc$model_fit$AUC_pres, na.rm=T)  #presence component
mean(fit_hmsc$model_fit$R2_abund, na.rm=T)  #positive component

## Load saved models
Mpres <- readRDS(paste0("Mpres_",save_model_suffix,".rds"))
Mposi <- readRDS(paste0("Mposi_",save_model_suffix,".rds"))

## Evaluate coefficients
postBeta = getPostEstimate(Mpres, parName="Beta")
par(mar=c(5,5,3,3)); plotBeta(Mpres, post=postBeta)
postBeta = getPostEstimate(Mposi, parName="Beta")
plotBeta(Mposi, post=postBeta)

## Variance partitioning
VPpres = computeVariancePartitioning(Mpres)
plotVariancePartitioning(Mpres, VPpres)  #lat 23, depth 21, sst 10, mld 4, lunar 3, area 2, catch 1, month 10, sample 11, vessel 15
VPabund = computeVariancePartitioning(Mposi)
plotVariancePartitioning(Mposi, VPabund)  #lat 6, depth 4, sst 3, area 1, catch 1, sample 83, vessel 2

## Residual species correlations
OmegaCor_pres = computeAssociations(Mpres)
OmegaCor_abund = computeAssociations(Mposi)
supportLevel = 0.95  #try lowering this to see correlations with less support

LV_lvl <- 1  #1 or 2 (sample and vessel; Mpres$nr)
toPlot_pres = ((OmegaCor_pres[[LV_lvl]]$support>supportLevel)
               + (OmegaCor_pres[[LV_lvl]]$support<(1-supportLevel))>0)*OmegaCor_pres[[LV_lvl]]$mean  #only plots corrs with > 95% support
par(mfrow=c(1,1), mar=c(2,2,2,2))
corrplot(toPlot_pres[order.single(toPlot_pres), order.single(toPlot_pres)],  #reordered by similar correlations
         method = "color", type="lower", col = colorRampPalette(c("blue","white","red"))(200),
         title = paste("Hurdle-presence; random effect level:", Mpres$rLNames[LV_lvl]),
         mar=c(0,0,1,0), tl.cex=0.6)
toPlot_abund = ((OmegaCor_abund[[LV_lvl]]$support>supportLevel)
                + (OmegaCor_abund[[LV_lvl]]$support<(1-supportLevel))>0)*OmegaCor_abund[[LV_lvl]]$mean  #only plots corrs with > 95% support
par(mfrow=c(1,1), mar=c(2,2,2,2))
corrplot(toPlot_abund[order.single(toPlot_abund), order.single(toPlot_abund)],  #reordered by similar correlations
         method = "color", type="lower", col = colorRampPalette(c("blue","white","red"))(200),
         title = paste("Hurdle-positive; random effect level:", Mposi$rLNames[LV_lvl]), mar=c(0,0,1,0),
         tl.cex=0.6)


###################################
##  Cross-Validation - Joint     ##
###################################

## Load saved models
Mpres <- readRDS(paste0("Mpres_",save_model_suffix,".rds"))
Mposi <- readRDS(paste0("Mposi_",save_model_suffix,".rds"))

## For saving results
save_r2_rmse <- data.frame(species=spp_list_hmsc, auc_pres=NA, r2_abund=NA,
                           r2_total=NA, rmse_total=NA, rmae_total=NA)

## Create partition across observations; both hurdle components MUST use same partition
partitionP = createPartition(hM = Mpres, nfolds = 5, column = "sample")
saveRDS(partitionP, paste0("Partition_",save_model_suffix,".rds"))

## Presence component
t1 <- Sys.time()
predsP = pcomputePredictedValues(hM = Mpres,
                                 partition = partitionP,
                                 nParallel = 1,
                                 verbose = 500,
                                 expected = FALSE,
                                 updater = list(GammaEta=FALSE))
t2 <- Sys.time(); t2-t1 
saveRDS(predsP, paste0("PresPreds_CV_",save_model_suffix,".rds"))

## Positive component
t1 <- Sys.time()
predsA = pcomputePredictedValues(hM = Mposi,
                                 partition = partitionP,
                                 nParallel = 1,
                                 verbose = 500,
                                 expected = TRUE)
t2 <- Sys.time(); t2-t1
saveRDS(predsA, paste0("PosiPreds_CV_",save_model_suffix,".rds"))

## Predictive performance of hurdle components
cvP_res <- evaluateModelFit(hM=Mpres, predY=predsP)
mean(cvP_res$AUC, na.rm=T)
save_r2_rmse$auc_pres <- cvP_res$AUC

cvA_res <- evaluateModelFit(hM=Mposi, predY=predsA)
mean(cvA_res$R2, na.rm=T)
save_r2_rmse$r2_abund <- cvA_res$R2

## Multiply predictions together for total biomass
preds_t <- predsA  #calculate total predicted (fitted) values
num_outliers <- 0
for (ff in 1:dim(predsA)[3]) {
  preds_t[,,ff] <- predsP[,,ff] * exp(predsA[,,ff])  #undo the log

  for (cc in 1:ncol(preds_t[,,ff])) {  #remove outlier values [only if abundance component is predicting some spurious values]
    datx <- preds_t[,cc,ff]
    exceed_max <- length(datx[datx > 3*spp_max[cc]])  # > 3*maximum observed biomass = outlier
    num_outliers <- num_outliers+exceed_max
    datx[datx > 3*spp_max[cc]] <- NA  #delete the outlier predictions
    preds_t[,cc,ff] <- datx
  }
}
num_outliers   #number of outliers removed

## Take mean (or *median*) of posterior for point estimates of total biomass
preds_t_mean <- apply(preds_t, c(1,2), mean, na.rm=T)
saveRDS(preds_t_mean, paste0("TotalPreds_CV_",save_model_suffix,".rds"))

## Plot some observed vs predicted biomasses
obs_dat_hmsc <- biomass_data
par(mfrow=c(1,3))
plot(obs_dat_hmsc[,1], preds_t_mean[,1], main=spp_list_hmsc[1]); abline(a=0,b=1)
plot(obs_dat_hmsc[,2], preds_t_mean[,2], main=spp_list_hmsc[2]); abline(a=0,b=1)
plot(obs_dat_hmsc[,3], preds_t_mean[,3], main=spp_list_hmsc[3]); abline(a=0,b=1)

## Predicted performance for total biomass
for (ss in 1:length(spp_list_hmsc)) {
  sppx <- spp_list_hmsc[ss]
  save_r2_rmse$r2_total[ss] <- R2(preds_t_mean[,ss], obs_dat_hmsc[,sppx])  #predicted then observed
  save_r2_rmse$rmse_total[ss] <- RMSE(preds_t_mean[,ss], obs_dat_hmsc[,sppx])
  save_r2_rmse$rmae_total[ss] <- RMAE(preds_t_mean[,ss], obs_dat_hmsc[,sppx])
}
mean(save_r2_rmse$auc_pres, na.rm=T)
mean(save_r2_rmse$r2_total)
mean(save_r2_rmse$rmae_total)
mean(save_r2_rmse$rmse_total)

saveRDS(save_r2_rmse, paste0("SummaryFit_CV_",save_model_suffix,".rds"))



##############################################
##  Cross-Validation - Conditional Joint    ##
##############################################

## Load saved models
Mpres <- readRDS(paste0("Mpres_",save_model_suffix,".rds"))
Mposi <- readRDS(paste0("Mposi_",save_model_suffix,".rds"))

## For saving results
save_r2_rmse <- data.frame(species=spp_list_hmsc, auc_pres=NA, r2_abund=NA,
                           r2_total=NA, rmse_total=NA, rmae_total=NA)

## Create new partition or use same one as previous joint prediction CV
partition = createPartition(hM = Mpres, nfolds = 5, column = "sample")
saveRDS(partitionP, paste0("Partition_",save_model_suffix,".rds"))

## Create the presence and abundance response matrices of known (retained) taxa to input
retained_spp_hmsc <- spp_list_hmsc[spp_list_hmsc %in% spp_list_ret]
partition.sp.df <- data.frame(species = spp_list_hmsc, value = 2)
partition.sp.df$value[partition.sp.df$species %in% retained_spp_hmsc] <- 1
partition.sp <- as.vector(partition.sp.df$value)
# Hmsc will estimate all taxa in partition 1 given observed biomasses of 2, and vice versa
# we will only report the predictions of partition 2 - the discarded taxa

## Presence component, conditional 
t1 <- Sys.time()
predsPc = pcomputePredictedValues(hM = Mpres,
                                  partition = partitionP,
                                  partition.sp = partition.sp,
                                  expected=FALSE,
                                  updater=list(GammaEta=FALSE),
                                  mcmcStep = 100,  #try mcmcStep = 5,10,100 and look for any change in AUC (more is better though)
                                  verbose = 500,
                                  nParallel = 1)
t2 <- Sys.time(); t2- t1
saveRDS(predsPc, paste0("PresPredsCOND_CV_",save_model_suffix,".rds"))

## Positive component, conditional
t1 <- Sys.time()
predsAc = pcomputePredictedValues(hM = Mposi,
                                  partition=partitionP,
                                  partition.sp = partition.sp,
                                  expected = TRUE,
                                  mcmcStep = 100,
                                  verbose = 500,
                                  nParallel = 1) 
t2 <- Sys.time(); t2- t1
saveRDS(predsAc, paste0("PosiPredsCOND_CV_",save_model_suffix,".rds"))

## Predictive performance of hurdle components
cvP_res <- evaluateModelFit(hM=Mpres, predY=predsPc)
mean(cvP_res$AUC, na.rm=T)
save_r2_rmse$auc_pres <- cvP_res$AUC
cvA_res <- evaluateModelFit(hM=Mposi, predY=predsAc)
mean(cvA_res$R2, na.rm=T)
save_r2_rmse$r2_abund <- cvA_res$R2

## Multiply predictions together for total biomass
preds_t <- predsAc  #calculate total predicted (fitted) values
num_outliers <- 0
for (ff in 1:dim(predsAc)[3]) {
  preds_t[,,ff] <- predsPc[,,ff] * exp(predsAc[,,ff])  #undo the log
  
  for (cc in 1:ncol(preds_t[,,ff])) {  #remove outliers if necessary
    datx <- preds_t[,cc,ff]
    exceed_max <- length(datx[datx > 3*spp_max[cc]])
    num_outliers <- num_outliers+exceed_max
    datx[datx > 3*spp_max[cc]] <- NA
    preds_t[,cc,ff] <- datx
  }
}
num_outliers   #number of outliers removed 

## Take mean (or *median*) of posterior for point estimates of total biomass
preds_t_mean <- apply(preds_t, c(1,2), mean, na.rm=T)
saveRDS(preds_t_mean, paste0("TotalPredsCOND_CV_",save_model_suffix,".rds"))

## Plot some observed vs conditional predicted biomasses
obs_dat_hmsc <- biomass_data
par(mfrow=c(1,3))
plot(obs_dat_hmsc[,1], preds_t_mean[,1], main=spp_list_hmsc[1]); abline(a=0,b=1)
plot(obs_dat_hmsc[,2], preds_t_mean[,2], main=spp_list_hmsc[2]); abline(a=0,b=1)
plot(obs_dat_hmsc[,3], preds_t_mean[,3], main=spp_list_hmsc[3]); abline(a=0,b=1)

## Predicted performance for total biomass
for (ss in 1:length(spp_list_hmsc)) {
  sppx <- spp_list_hmsc[ss]
  save_r2_rmse$r2_total[ss] <- R2(preds_t_mean[,ss], obs_dat_hmsc[,sppx])  #predicted then observed
  save_r2_rmse$rmse_total[ss] <- RMSE(preds_t_mean[,ss], obs_dat_hmsc[,sppx])
  save_r2_rmse$rmae_total[ss] <- RMAE(preds_t_mean[,ss], obs_dat_hmsc[,sppx])
}
mean(save_r2_rmse$auc_pres, na.rm=T)
mean(save_r2_rmse$r2_total)
mean(save_r2_rmse$rmae_total)
mean(save_r2_rmse$rmse_total)

saveRDS(save_r2_rmse, paste0("SummaryFitCOND_CV_",save_model_suffix,".rds"))



