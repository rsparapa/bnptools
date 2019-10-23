
lpml <- function(y,    # data
                 mus,  # Posterior samples for the mean
                 sigs  # Posterior samples for sigma
                 ){
  S <- nrow(mus)
  n <- length(y)
  maxs <- NULL
  f <- NULL; finv <- NULL;CPO <- NULL
  for( i in 1:n){
    for(s in 1:S){
      #f[s] <- dnorm(y[i],mus[s,i],sigs[s]) 
      f[s] <-(1/(sqrt(2*pi*sigs[s]^2)))*exp(-(1/(2*sigs[s]^2))*(y[i]-mus[s,i])^2)
      finv[s] <- 1/f[s]  #check max (-log(f[s])) 
    }
    CPO[i] <- (mean(finv))^(-1)
    maxs[i] <-  max(-log(f)) 
  }
  LPML <- sum(log(CPO))
  return(list(maxs = maxs,
              CPOi = CPO,
              LPML = LPML))
}

