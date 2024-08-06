#' Estimate the signal's parameters in each segment given a set of change points 
#'
#' @param y vector of observations
#' @param Xsin matrix of sine functions X^(1) based on the grid
#' @param Xcos matrix of cosine functions X^(2) based on the grid
#' @param Ne number of effects in the SuSiE algorithm (can be different from the one used for segmentation, if a better estimate of the signal is desired)
#' @param cp estimated change points from the proposed procedure
#' @param sigma2 default is 0, meaning that the variance of the error term is estimated. It can be fixed to a specific value.
#' @param pl if \code{TRUE}, plot the estimates over the observed values

#'
#' @return  \code{models}: list of models estimated in each segment (output of \code{SuSiE_group})
#' @return  \code{signal_estimate}: estimate of the signal, given the change points
#' @export
fit_segments = function(y,Xsin,Xcos,Ne,cp,sigma2=0,pl=TRUE){
  L = Ne
  
  out = list()
  yhat = NULL
  
  n = length(y)
  n_cp = length(cp)
  pts = c(0,cp,n)
  
  for (i in 1:(n_cp+1)) {
    segm = (pts[i]+1):pts[i+1]
    
    yW = y[segm]
    XsinW = Xsin[segm,] 
    XcosW = Xcos[segm,]
    
    elbos = rep(NA,L)
    for (l in 1:L) {
      fitW = SuSiE_group(yW,XsinW,XcosW,L=l,sigma2=sigma2)
      elbos[l] = fitW$elbo
    }
    maxlbf = max(elbos)
    w = exp(elbos - maxlbf);
    w_weighted = w;
    weighted_sum_w = sum(w_weighted);
    alpha = w_weighted/weighted_sum_w;
    lhat = which.max(alpha)
    
    fitW = SuSiE_group(yW,XsinW,XcosW,L=lhat,sigma2=sigma2)
    
    out[[i]] = fitW
    yhat = c(yhat,fitW$yhat)
  }
  
  if (pl) {
    plot(y)
    abline(v=cp,lty='dashed',lwd=2)
    lines(yhat,col=2,lwd=2)
  }
  
  list(models = out,
       signal_estimate = yhat)
}
