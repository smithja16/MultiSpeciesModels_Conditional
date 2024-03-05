####################################################################
###     ~  Fitting and Predicting Joint and Stacked SDMs   ~     ###
###  - An evaluation of conditional prediction, particularly...  ###
###      predicting fishery discards using retained catch.       ###
###  - James A. Smith Feb 2024                                   ###
####################################################################

## This script contains some basic performance functions

## Performance functions
RMSE = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  sqrt(mean((p - o)^2))
}

RMAE = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  (sum(abs(p - o))/length(p)) / mean(o)
}

R2 = function(p, o) {
  if (length(o[is.na(o)]) > 0) {
    nas <- which(is.na(o))  #remove NAs
    p <- p[-nas]
    o <- o[-nas]
  }
  round(cor(p,o)^2,3)
}
