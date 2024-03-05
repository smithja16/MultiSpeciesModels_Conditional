##################################################################
###   Fitting, Predicting, Evaluating Joint and Stacked SDMs   ###  
###    - As used for the Ocean Prawn Trawl observer data       ###
###    - Publication: xxx                                      ###
###    - James A. Smith Feb 2024                               ###
##################################################################

## This script contains a function to create folds for k-folds cross-validation

# Every fold*repeat is a part of a list, which IDs the rows to be withheld
CreateFolds <- function(data, n_folds, n_repeats, seed) {
  n_folds <- n_folds
  n_repeats <- n_repeats
  n_per_fold <- floor(nrow(data)/n_folds)
  
  folds <- data.frame(fold=rep(1:n_folds, n_repeats), rep=rep(1:n_repeats, each=n_folds))
  fold_rows <- list()
  set.seed(seed)
  for (rr in 1:n_repeats) {
    rows <- data.frame(idx=1:nrow(data),
                       ran=0)
    rows$ran <- runif(nrow(data))
    rows <- rows[order(rows$ran),]
    for (kk in 1:n_folds) {
      idx1 <- kk*n_per_fold-(n_per_fold-1)
      idx2 <- kk*n_per_fold
      fold_rows[[kk+(rr*n_folds-n_folds)]] <- rows$idx[idx1:idx2]  #list is indexed by rows in 'folds'
    }
  }
  return(list(folds=folds, fold_rows=fold_rows))
}
