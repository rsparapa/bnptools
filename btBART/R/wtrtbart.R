
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

wtrtbart=function(
  x.train1, trt.train,y.train, 
  x.test1=matrix(0.0,0,0),trt.test = 0, 
  mu0=0,sig20=25,
  sparse=FALSE, theta=0, omega=1,
  a=0.5, b=1, augment=FALSE, rho=NULL,
  x1info=matrix(0.0,0,0),usequants=FALSE,
  cont=FALSE, rm.const=TRUE,
  sigest=NA, sigdf=3, sigquant=.90,
  k=2.0, power=2.0, base=.95,
  sigmaf=NA, lambda=NA,
  fmean=mean(y.train),
  w=rep(1,length(y.train)),
  ntree=200L, numcut1=100L, 
  ndpost=1000L, nskip=100L, keepevery=1L,
  nkeeptrain=ndpost, nkeeptest=ndpost,
  nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
  printevery=100L, transposed=FALSE
)
{
  #--------------------------------------------------
  #data
  n = length(y.train)
  
  if(!transposed){
    temp1 = bartModelMatrix(x.train1, numcut1, usequants=usequants,
                            cont=cont, xinfo=x1info, rm.const=rm.const)
    x.train1 = t(temp1$X)
    numcut1 = temp1$numcut
    x1info = temp1$xinfo
    if(length(x.test1)>0)
      x.test1 = t(bartModelMatrix(x.test1[ , temp1$rm.const]))
    rm.const <- temp1$rm.const
    grp <- temp1$grp
    rm(temp1)
    
  }
  else {
    rm.const <- NULL
    grp <- NULL
  }
  
  if(n!=ncol(x.train1))
    stop('The length of y.train and the number of rows in x.train must be identical')
  
  x.train <- rbind(trt.train,x.train1)
  
  p1 = nrow(x.train1) 
  np = ncol(x.test1)
  p = p1+1  
  if(length(rho)==0) rho=(p1+1)
  if(length(rm.const)==0) rm.const <- 1:(p1+1)
  if(length(grp)==0) grp <- 1:(p1+1)
  
  ##if(p>1 & length(numcut)==1) numcut=rep(numcut, p)
  
  y.train = y.train-fmean
  #--------------------------------------------------
  #set nkeeps for thinning
  if((nkeeptrain!=0) & ((ndpost %% nkeeptrain) != 0)) {
    nkeeptrain=ndpost
    cat('*****nkeeptrain set to ndpost\n')
  }
  if((nkeeptest!=0) & ((ndpost %% nkeeptest) != 0)) {
    nkeeptest=ndpost
    cat('*****nkeeptest set to ndpost\n')
  }
  if((nkeeptestmean!=0) & ((ndpost %% nkeeptestmean) != 0)) {
    nkeeptestmean=ndpost
    cat('*****nkeeptestmean set to ndpost\n')
  }
  if((nkeeptreedraws!=0) & ((ndpost %% nkeeptreedraws) != 0)) {
    nkeeptreedraws=ndpost
    cat('*****nkeeptreedraws set to ndpost\n')
  }
  #--------------------------------------------------
  #prior
  nu=sigdf
  if(is.na(lambda)) {
    if(is.na(sigest)) {
      if(p < n) {
        df = data.frame(t(x.train),y.train)
        lmf = lm(y.train~.,df)
        sigest = summary(lmf)$sigma
      } else {
        sigest = sd(y.train)
      }
    }
    qchi = qchisq(1.0-sigquant,nu)
    lambda = (sigest*sigest*qchi)/nu #lambda parameter for sigma prior
  } else {
    sigest=sqrt(lambda)
  }
  
  if(is.na(sigmaf)) {
    tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
  } else {
    tau = sigmaf/sqrt(ntree)
  }
  #--------------------------------------------------
  #call
  res = .Call("cwtrtbart",
              n,  #number of observations in training data
              p1,  #dimension of x 
              np, #number of observations in test data
              x.train1,   #pxn training data x
              trt.train,
              y.train,   #pxn training data x
              x.test1,   #p*np test data x
              trt.test,
              mu0, 
              sig20,
              ntree,
              numcut1, 
              ndpost*keepevery,
              nskip,
              power,
              base,
              tau,
              nu,
              lambda,
              sigest,
              w,
              sparse,
              theta,
              omega,
              grp,
              a,
              b,
              rho,
              augment,
              nkeeptrain,
              nkeeptest,
              nkeeptestmean,
              nkeeptreedraws,
              printevery,
              x1info 
  )
  res$mu = fmean 
  res$yhat.train.mean = res$yhat.train.mean+fmean
  res$yhat.train = res$yhat.train+fmean
  res$yhat.test.mean = res$yhat.test.mean+fmean
  res$yhat.test = res$yhat.test+fmean
  res$beta <- res$beta[(nskip+1):(ndpost+nskip)]
  if(nkeeptreedraws>0)
    names(res$treedraws1$cutpoints) = dimnames(x.train1)[[1]] 
  
  dimnames(res$varcount1)[[2]] = as.list(dimnames(x.train1)[[1]])
  dimnames(res$varprob1)[[2]] = as.list(dimnames(x.train1)[[1]]) 
  ##res$nkeeptreedraws=nkeeptreedraws
  res$varcount1.mean <- apply(res$varcount1, 2, mean)
  res$varprob1.mean <- apply(res$varprob1, 2, mean) 
  res$rm.const <- rm.const
  attr(res, 'class') <- 'wtrtbart'
  return(res)
  
}
