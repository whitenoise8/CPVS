#' Get the design matrices X^(1) and X^(2)
#'
#' @param x vector of times
#' @param om the grid of frequencies

#'
#' @return  \code{Xsin}: matrix X^(1) 
#' @return  \code{Xcos}: matrix X^(2) 
#' @export
getX = function(x,om){
  
  n = length(x)
  Xsin = matrix(0,nrow=n,ncol=length(om))
  Xcos = matrix(0,nrow=n,ncol=length(om))
  
  for (j in 1:length(om)){
    Xsin[,j] = sin(om[j]*2*pi*x)
    Xcos[,j] = cos(om[j]*2*pi*x)
  }
  
  list(Xsin=Xsin, Xcos=Xcos)
}