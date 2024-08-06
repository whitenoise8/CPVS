#' Segmentation for univariate or multiple time series sharing same change point locations using Optimistic Search and gain function based on Marginal Likelihood/ELBO (see paper for details)
#'
#' @param y vector of observations or matrix (\code{T x d}) if multivariate case
#' @param Xsin matrix (\code{T x p}) or array (\code{T x p x d}), if multivariate, of sine functions X^(1) based on the grid
#' @param Xcos matrix (\code{T x p}) or array (\code{T x p x d}), if multivariate, of cosine functions X^(2) based on the grid
#' @param Ne number of effects in the SuSiE algorithm
#' @param IC the information criterion to be used. "MDL" (minimum description length) and "mBIC" modified BIC are possible
#' @param sigma2 default is 0, meaning that the variance of the error term is estimated. It can be fixed to a specific value.
#' @param nu default is 1/2. \nu parameter in the optimistic search algorithm
#' @param thresh default is equal to 0. Threshold in the gain function to define a change point
#' @param rate_c defines the minimum distance from the boundaries for a change point rate_c*T
#' @param option the type of pruning to perform (more details...)

#'
#' @return  The estimated set of change points
#' @export
cpvs = function(y,Xsin,Xcos,Ne,IC,sigma2=0,nu=1/2,thresh=0,rate_c=0.1,option='best IC') {
  if (is.null(ncol(y))) {
    cp_set = getCandidateSetU(y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,nu=nu,thresh=thresh,rate_c=rate_c,Trace=0)
    out = pruneSetU(cp_set=cp_set,y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,IC=IC,option=option,Trace=0,pl=FALSE,re_test=3,fixCP=0,maxEval=8)
  }
  
  if (!is.null(ncol(y))) {
    cp_set = getCandidateSetMV(Y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,nu=nu,thresh=thresh,rate_c=rate_c,Trace=0)
    out = pruneSetMV(cp_set=cp_set,Y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,IC=IC,option=option,Trace=0,pl=FALSE,re_test=3,fixCP=0,maxEval=8)
  }
  
  out
}
