#' Prune the set of candidate change points using an information criterion (IC) (see paper for details)
#'
#' @param cp_set output of \code{getCandidateSet}
#' @param y vector of observations
#' @param Xsin matrix of sine functions X^(1) based on the grid
#' @param Xcos matrix of cosine functions X^(2) based on the grid
#' @param Ne number of effects in the SuSiE algorithm
#' @param sigma2 default is 0, meaning that the variance of the error term is estimated. It can be fixed to a specific value.
#' @param IC the information criterion to be used. "MDL" (minimum description length) and "mBIC" modified BIC are possible
#' @param option the type of pruning to perform (more details...)
#' @param Trace if 1, print progress of the pruning
#' @param pl if \code{TRUE}, plot the IC path

#'
#' @return  The estimated set of change points
#' @export
pruneSet = function(cp_set,y,Xsin,Xcos,Ne,sigma2=0,IC="MDL",option="best IC",Trace=0,pl=FALSE,re_test=3,fixCP=0,maxEval=8) {
  if (is.null(ncol(y))) {
    out = pruneSetU(cp_set=cp_set,y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,IC=IC,option=option,Trace=Trace,pl=pl,re_test=3,fixCP=0,maxEval=8)
  }
  
  if (!is.null(ncol(y))) {
    out = pruneSetMV(cp_set=cp_set,Y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,IC=IC,option=option,Trace=Trace,pl=pl,re_test=3,fixCP=0,maxEval=8)
  }
  
  out
}
