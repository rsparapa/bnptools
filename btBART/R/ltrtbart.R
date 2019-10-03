
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

ltrtbart=function(
  x.train1,trt.train,y.train,
  x.test1=matrix(0.0,0,0),trt.test=0,
  mu0=0,sig20=25,
  sparse=FALSE, a=0.5, b=1, augment=FALSE, rho=NULL,
  x1info=matrix(0.0,0,0),usequants=FALSE,
  cont=FALSE, rm.const=TRUE, tau.interval=0.95,
  k=2.0, power=2.0, base=.95,
  binaryOffset=NULL,
  ntree=200L, numcut1=100L, 
  ndpost=1000L, nskip=100L,
  keepevery=1L,
  nkeeptrain=ndpost, nkeeptest=ndpost,
  #nkeeptestmean=ndpost,
  nkeeptreedraws=ndpost,
  printevery=100, transposed=FALSE
  ##treesaslists=FALSE
)
{
  #--------------------------------------------------
  #data
  n = length(y.train)
  
  # if(binaryOffset!=0)
  #     stop('binaryOffset not supported by lbart')
  
  if(length(binaryOffset)==0) binaryOffset=qlogis(mean(y.train))
  
  if(!transposed) {
    temp1 = bartModelMatrix(x.train1, numcut1, usequants=usequants,
                            cont=cont, xinfo=x1info, rm.const=rm.const)
    x.train1 = t(temp1$X)
    numcut1 = temp1$numcut
    x1info = temp1$xinfo
     
    
    if(length(x.test1)>0) 
      x.test1 = t(bartModelMatrix(x.test1[ , temp1$rm.const]))
      rm.const <- temp1$rm.const  
    rm(temp1) 
  }
  
  if(n!=ncol(x.train1))
    stop('The length of y.train and the number of rows in x.train1 must be identical')
  
  p1 = nrow(x.train1) 
  np = ncol(x.test1)
  if(length(rho)==0) rho <- p1
  if(length(rm.const)==0) rm.const <- 1:p1
  if(tau.interval>0.5) tau.interval=1-tau.interval
  
  tau=qlogis(1-0.5*tau.interval)/(k*sqrt(ntree))
  
  #--------------------------------------------------
  #set  nkeeps for thinning
  if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
    nkeeptrain=ndpost
    cat('*****nkeeptrain set to ndpost\n')
  }
  if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
    nkeeptest=ndpost
    cat('*****nkeeptest set to ndpost\n')
  }
  ## if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
  ##    nkeeptestmean=ndpost
  ##    cat('*****nkeeptestmean set to ndpost\n')
  ## }
  if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
    nkeeptreedraws=ndpost
    cat('*****nkeeptreedraws set to ndpost\n')
  }
  #--------------------------------------------------
  #prior
  ## nu=sigdf
  ## if(is.na(lambda)) {
  ##    if(is.na(sigest)) {
  ##       if(p < n) {
  ##          df = data.frame(t(x.train),y.train)
  ##          lmf = lm(y.train~.,df)
  ##          rm(df)
  ##          sigest = summary(lmf)$sigma
  ##       } else {
  ##          sigest = sd(y.train)
  ##       }
  ##    }
  ##    qchi = qchisq(1.0-sigquant,nu)
  ##    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
  ## }
  
  ## if(is.na(sigmaf)) {
  ##    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree));
  ## } else {
  ##    tau = sigmaf/sqrt(ntree)
  ## }
  #--------------------------------------------------
  #call
  res = .Call("cltrtbart",
              n,  #number of observations in training data
              p1,  #dimension of x1 
              np, #number of observations in test data
              x.train1,   #p*n training data x
              trt.train,    
              as.integer(2*y.train-1),   #n*1 training data y
              x.test1,    #p*np test data x
              trt.test,    
              mu0,
              sig20,
              ntree,
              numcut1, 
              ndpost*keepevery,
              nskip,
              power,
              base,
              binaryOffset,
              tau,
              sparse,
              a,
              b,
              rho,
              augment,
              nkeeptrain,
              nkeeptest,
              nkeeptreedraws,
              printevery,
              x1info 
  )
  
  if(nkeeptrain>0) {
    ##res$yhat.train.mean <- NULL
    ##res$yhat.train.mean = res$yhat.train.mean+binaryOffset
    res$beta <- res$beta[(nskip+1):(ndpost+nskip)]
    res$yhat.train = res$yhat.train+binaryOffset
    res$prob.train = plogis(res$yhat.train)
    res$prob.train.mean <- apply(res$prob.train, 2, mean)
  } else {
    res$yhat.train <- NULL
    ##res$yhat.train.mean <- NULL
  }
  
  if(np>0) {
    ##res$yhat.test.mean <- NULL
    ##res$yhat.test.mean = res$yhat.test.mean+binaryOffset
    res$yhat.test = res$yhat.test+binaryOffset
    res$prob.test = plogis(res$yhat.test)
    res$prob.test.mean <- apply(res$prob.test, 2, mean)
  } else {
    res$yhat.test <- NULL
    ##res$yhat.test.mean <- NULL
  }
  
  if(nkeeptreedraws>0)
    names(res$treedraws1$cutpoints) = dimnames(x.train1)[[1]] 
  
  dimnames(res$varcount1)[[2]] = dimnames(x.train1)[[1]]
  dimnames(res$varprob1)[[2]] = dimnames(x.train1)[[1]] 
  res$varcount1.mean <- apply(res$varcount1, 2, mean)
  res$varprob1.mean <- apply(res$varprob1, 2, mean) 
  res$binaryOffset=binaryOffset
  res$rm.const <- rm.const
  attr(res, 'class') <- 'ltrtbart'
  return(res)
}
