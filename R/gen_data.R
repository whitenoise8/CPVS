#' Generate data based on the functions used in the simulation study in the paper
#'
#' @param n length of the series
#' @param d number of series
#' @param ncp number of change points
#' @param minLen minimum number of observations between consecutive change points
#' @param sigma standard deviation of the errors: if scalar, the same value is applied for all segments
#' @param fix_seed seed can be fixed for reproducibility. If 0, no seed is fixed

#'
#' @return  \code{y}: the simulated data (a matrix \code{n x d} if multivariate case)
#' @return  \code{signal}: the true simulated signal (a matrix \code{n x d} if multivariate case)
#' @return  \code{cp}: the change point locations
#' @export
generateData = function(n,d,ncp,minLen,sigma=1,fix_seed=0) {
  
  #: Define the setting and the signal 
  signal = function(x,beta1,beta2,omega) {
    n = length(x)
    M = length(omega)
    S = matrix(NA,n,M)
    for (j in 1:M) S[,j] = beta1[j]*sin(2*pi*x*omega[j])+beta2[j]*cos(2*pi*x*omega[j])
    rowSums(S)
  }
  
  setting = list()
  setting[[1]] = matrix(c(2.0, 3.5, 1/30), 1, 3, byrow=TRUE)
  setting[[2]] = matrix(c(4.0, 3.0, 1/12), 1, 3, byrow=TRUE)
  
  setting[[3]] = matrix(c(1.5, 2.0, 1/35,
                          2.5, 3.0, 1/20), 2, 3, byrow=TRUE)
  setting[[4]] = matrix(c(2.5, 4.0, 1/22,
                          4.0, 2.0, 1/15), 2, 3, byrow=TRUE)
  
  setting[[5]] = matrix(c(3.0, 2.0, 1/40,
                          2.0, 4.0, 1/20,
                          1.0, 1.5, 1/10), 3, 3, byrow=TRUE)
  setting[[6]] = matrix(c(2.0, 3.0, 1/24,
                          4.0, 5.0, 1/15,
                          1.0, 2.5, 1/7), 3, 3, byrow=TRUE)
  
  colnames(setting[[1]]) = colnames(setting[[2]]) = colnames(setting[[3]]) = c("beta_1","beta_2","omega")
  colnames(setting[[4]]) = colnames(setting[[5]]) = colnames(setting[[6]]) = c("beta_1","beta_2","omega")
  
  if (d == 1) {
    #: Preliminaries
    if (fix_seed > 0) set.seed(fix_seed)
    nseg = ncp + 1
    if (length(sigma) == 1) sigma = rep(sigma,nseg)
    
    cp = rep(NA,ncp)
    id = rep(NA,nseg)
    sets = rep(NA,ncp+2)
    
    #: Get change-points
    id[1:nseg] = 1
    while (any(diff(id[1:nseg])==0)) id[1:nseg] = sample(1:6,nseg,replace=T)
    
    tau = rep(2,ncp)
    while (min(diff(c(1,tau,n))) < minLen) tau = sort(sample(1:n,ncp))
    cp[1:ncp] = tau
    sets[1:(nseg+1)] = c(0,tau,n)
    
    #: Sampling
    knots = c(0,tau,n)
    y = NULL
    s = NULL
    for (i in 1:nseg) {
      x = seq((knots[i]+1),knots[i+1])
      s_i = signal(x,setting[[id[i]]][,"beta_1"],setting[[id[i]]][,"beta_2"],setting[[id[i]]][,"omega"])
      s = c(s,s_i)
      y = c(y,s_i+rnorm(length(x),0,sigma[i]))
    }
    
    out = list(y=y,signal=s,cp=cp) 
  }
  
  
  if (d > 1) {
    #: Preliminaries
    if (fix_seed > 0) set.seed(fix_seed)
    nseg = ncp + 1
    if (length(sigma) == 1) sigma = rep(sigma,nseg)
    
    cp = rep(NA,ncp)
    sets = rep(NA,ncp+2)
    
    tau = rep(2,ncp)
    while (min(diff(c(1,tau,n))) < minLen) tau = sort(sample(1:n,ncp))
    cp[1:ncp] = tau
    sets[1:(nseg+1)] = c(0,tau,n)
    
    Y = S = matrix(NA,n,d)
    for (m in 1:d) {
      id = rep(NA,nseg)
      
      #: Get change-points
      id[1:nseg] = 1
      while (any(diff(id[1:nseg])==0)) id[1:nseg] = sample(1:6,nseg,replace=T)
      
      #: Sampling
      knots = c(0,tau,n)
      y = NULL
      s = NULL
      for (i in 1:nseg) {
        x = seq((knots[i]+1),knots[i+1])
        s_i = signal(x,setting[[id[i]]][,"beta_1"],setting[[id[i]]][,"beta_2"],setting[[id[i]]][,"omega"])
        s = c(s,s_i)
        y = c(y,s_i+rnorm(length(x),0,sigma[i]))
      }
      
      Y[,m] = y
      S[,m] = s
    }
    
    out = list(y=Y,signal=S,cp=cp) 
  }
  
  out
}