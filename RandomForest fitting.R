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


###########################
## Fit full hurdle model ##
###########################

rf_marginal <- fit_rf_hurdle_cond(spp_list_response = spp_list_dis,  
                                  spp_list_covar = spp_list_ret,
                                  first_taxon = spp_list_rf[1],
                                  data_train_pres = all_data_pres,
                                  data_train = all_data,
                                  data_test_pres = NA,  #no test data for full model
                                  data_test = NA,
                                  pres_threshold = 0,
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
                                     pres_threshold = 0,
                                     pres_prevalence = T,
                                     conditional = T,  #for conditional = T
                                     save_models = T)  #save individual models for evaluation
