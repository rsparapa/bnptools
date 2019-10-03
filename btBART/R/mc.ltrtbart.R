
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

mc.ltrtbart=function(
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
  printevery=100, keeptrainfits=TRUE,transposed=FALSE,
  mc.cores = 2L,nice=19L,seed=99L
  ##treesaslists=FALSE
)
{
  
  if(.Platform$OS.type!='unix')
    stop('parallel::mcparallel/mccollect do not exist on windows')
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  
  
  #--------------------------------------------------
  #data
 # n = length(y.train)
  
  # if(binaryOffset!=0)
  #     stop('binaryOffset not supported by lbart')
  
  #if(length(binaryOffset)==0) binaryOffset=qlogis(mean(y.train))
  
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
  
  mc.cores.detected <- detectCores()
  
  if(mc.cores>mc.cores.detected) mc.cores <- mc.cores.detected
  ## warning(paste0('The number of cores requested, mc.cores=', mc.cores,
  ##                ',\n exceeds the number of cores detected via detectCores() ',
  ##                'which yields ', mc.cores.detected, ' .'))
  
  mc.ndpost <- ceiling(ndpost/mc.cores)
  
  
  
  for(i in 1:mc.cores) {
    parallel::mcparallel({psnice(value=nice);
     ltrtbart(x.train1=x.train1,trt.train=trt.train,y.train=y.train,
        x.test1=x.test1,trt.test=trt.test,
        sparse=sparse,
        mu0=mu0,sig20=sig20,
        a=a, b=b, augment=augment, rho=rho,
        x1info=x1info, 
        k=k, power=power, base=base,
        binaryOffset=binaryOffset,
        ntree=ntree, numcut1=numcut1,
        ndpost=mc.ndpost, nskip=nskip, keepevery=keepevery, 
        nkeeptreedraws=ndpost,
        printevery=printevery, transposed=TRUE )}, 
      silent=(i!=1))
    
    
    
    ## to avoid duplication of output
    ## capture stdout from first posterior only
  }
  
  post.list <- parallel::mccollect()
  
  post <- post.list[[1]]
  
  if(mc.cores==1 | attr(post, 'class')!='ltrtbart') return(post)
  else {
    if(class(rm.const)!='logical') post$rm.const <- rm.const
    
    p1 <- nrow(x.train1[post$rm.const, ])
    ##p <- nrow(x.train[ , post$rm.const])
    
    ## if(length(rm.const)==0) rm.const <- 1:p
    ## post$rm.const <- rm.const
    
    old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree), ' ', as.character(p1))
    ##old.text <- paste0(as.character(mc.nkeep), ' ', as.character(ntree), ' ', as.character(p1))
    old.stop <- nchar(old.text)
    
    post$treedraws1$trees <- sub(old.text,
                                paste0(as.character(mc.cores*mc.ndpost), ' ', as.character(ntree), ' ',
                                       ##paste0(as.character(mc.cores*mc.nkeep), ' ', as.character(ntree), ' ',
                                       as.character(p1)),
                                post$treedraws1$trees)
    
    keeptestfits <- length(x.test1)>0
    
    for(i in 2:mc.cores) {
      if(keeptrainfits) {
        post$yhat.train <- rbind(post$yhat.train, post.list[[i]]$yhat.train)
        post$prob.train <- rbind(post$prob.train, post.list[[i]]$prob.train)
      }
      
      post$beta <- cbind(post$beta, post.list[[i]]$beta)
      
      if(keeptestfits) {
        post$yhat.test <- rbind(post$yhat.test, post.list[[i]]$yhat.test)
        post$prob.test <- rbind(post$prob.test, post.list[[i]]$prob.test)
      }
      
    
      post$varcount1 <- rbind(post$varcount1, post.list[[i]]$varcount)
      post$varprob1 <- rbind(post$varprob1, post.list[[i]]$varprob)
      
      post$treedraws1$trees <- paste0(post$treedraws1$trees,
                                     substr(post.list[[i]]$treedraws1$trees, old.stop+2,
                                            nchar(post.list[[i]]$treedraws1$trees)))
      
      ## if(treesaslists) post$treedraws$lists <-
      ##                      c(post$treedraws$lists, post.list[[i]]$treedraws$lists)
    }
    
    ## if(length(post$yhat.train.mean)>0)
    ##     post$yhat.train.mean <- apply(post$yhat.train, 2, mean)
    
    ## if(length(post$yhat.test.mean)>0)
    ##     post$yhat.test.mean <- apply(post$yhat.test, 2, mean)
    
    post$varcount1.mean <- apply(post$varcount1, 2, mean)
    post$varprob1.mean <- apply(post$varprob1, 2, mean)
    
    attr(post, 'class') <- 'ltrtbart'
    return(post)
  }
}

