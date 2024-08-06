#' Get the candidate set of change points using Optimistic Search and gain function based on Marginal Likelihood/ELBO (see paper for details)
#'
#' @param y vector of observations or matrix (\code{T x d}) if multivariate case
#' @param Xsin matrix (\code{T x p}) or array (\code{T x p x d}), if multivariate, of sine functions X^(1) based on the grid
#' @param Xcos matrix (\code{T x p}) or array (\code{T x p x d}), if multivariate, of cosine functions X^(2) based on the grid
#' @param Ne number of effects in the SuSiE algorithm
#' @param sigma2 default is 0, meaning that the variance of the error term is estimated. It can be fixed to a specific value.
#' @param nu default is 1/2. \nu parameter in the optimistic search algorithm
#' @param thresh default is equal to 0. Threshold in the gain function to define a change point
#' @param rate_c defines the minimum distance from the boundaries for a change point rate_c*T
#' @param Trace if 1, print progress of the search

#'
#' @return  \code{cp_list}: list of the candidate change point locations
#' @return  \code{delbo}: list of the gain function values for each element in \code{cp_list}
#' @export
getCandidateSet = function(y,Xsin,Xcos,Ne,sigma2=0,nu=1/2,thresh=0,rate_c=0.1,Trace=0) {
  if (is.null(ncol(y))) {
    out = getCandidateSetU(y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,nu=nu,thresh=thresh,rate_c=rate_c,Trace=Trace)
  }
  
  if (!is.null(ncol(y))) {
    out = getCandidateSetMV(Y=y,Xsin=Xsin,Xcos=Xcos,Ne=Ne,sigma2=sigma2,nu=nu,thresh=thresh,rate_c=rate_c,Trace=Trace)
  }
  
  out
}
