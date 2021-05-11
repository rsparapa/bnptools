
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2019 Rodney Sparapani

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

mxbart=function(
               x.train, y.train,
               x.test=matrix(0,0,0), type='wbart',
               id.train=NULL,z.train=NULL,
               ntype=as.integer(
                   factor(type, levels=c('wbart', 'pbart'))),
                   ##factor(type, levels=c('wbart', 'pbart', 'lbart'))),
               sparse=FALSE, theta=0, omega=1,
               a=0.5, b=1, augment=FALSE, rho=NULL,
               xinfo=matrix(0,0,0), usequants=FALSE,
               rm.const=TRUE,
               sigest=NA, sigdf=3, sigquant=0.90,
               k=2, power=2, base=0.95,
               ##sigmaf=NA,
               lambda=NA, tau.num=c(NA, 3, 6)[ntype],
               mxps=NULL,
               offset=NULL, ##w=rep(1, length(y.train)),
               ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
               ndpost=1000L, nskip=100L,
               keepevery=c(1L, 10L, 10L)[ntype],
               printevery=100L, transposed=FALSE,
               hostname=FALSE
               )
{
    if(type=='wbart') ntype <- 1
    else if(type=='pbart') ntype <- 2

    n = length(y.train)
    if(length(id.train)==0)
        stop("the random effects indices must be provided")
    #### Process id.train; want i passed as an n x (l-1) matrix where 'l' is the number of level
    ## id.train is passed as list
    if(methods::is(id.train,'list')) {
        id.train.no <- length(id.train)
        id.train <- do.call(cbind,id.train)
    }
    ## id.train is passed as matrix
    else if(methods::is(id.train,'matrix')) {
        id.train.no <- NCOL(id.train)
        if(nrow(id.train)!=n) stop('id.train passed as a matrix must have the same number of rows as elements of y.train')
    }
    ## id.train is passed as vector (2-level model)
    else if(methods::is(id.train,'vector')) {
        id.train.no <- 1
        if(length(id.train)!=n) stop('id.train passed as a vector must have the same number of elements as y.train')
    }
    #### Process z.train;
    if(length(z.train)==0) {
        z.train <- vector(mode='list',length=id.train.no)
        for(i in 1:id.train.no) z.train[[i]] <- 1
    }
    ## 3+ level models
    if(id.train.no!=1){
        ## z.train isn't passed as a list
        if(!methods::is(z.train,'list')) stop("z.train must be passed as a list for 3+ level models")
        ## z.train is passed as a list
        else if(methods::is(z.train,'list')){
            z.cols <- sapply(z.train,NCOL)
            if(length(z.train)!=id.train.no) stop("the number of elements in z.train must be the same as the number of elements in id.train (if list) or as the number of columns in id.train (if matrix)")
            z.train <- do.call(cbind,z.train)
        }
    }
    ## 2-level model
    else{
        if(methods::is(z.train,'list')){
            if(length(z.train)!=1) stop("the number of elements in z.train must be the same as the number of elements in id.train (if list) or as the number of columns in id.train (if matrix)")
            z.train <- z.train[[1]]
        }
        z.cols <- NCOL(z.train)
    }
    c.train <- matrix(0,nrow=n,ncol=NCOL(id.train))
    r.train <- list()
    L <- integer(id.train.no)
    n.train <- matrix(0,nrow=n,ncol=NCOL(id.train))
    for(i in 1:NCOL(id.train)){
        c.train[,i] <- integer(n) ## changing from R/arbitrary indexing to C/0
        if(id.train.no!=1) {
            r.train[[i]] <- unique(id.train[,i])
            for(ii in 1:n) c.train[ii,i]=which(id.train[ii,i]==r.train[[i]])-1
        }
        else {
            r.train[[1]] <- unique(id.train)
            for(ii in 1:n) c.train[ii,1]=which(id.train[ii]==r.train[[1]])-1
        }
        r.train[[i]] <- unique(c.train[,i])
        L[i] <- length(r.train[[i]])
        tmp.n <- as.data.frame(table(c.train[,i]))$Freq
        n.train[,i] <- rep(tmp.n,times=tmp.n)
        rm(tmp.n)
        z.train <- matrix(as.numeric(z.train),nrow=n)      
    }
    if(length(mxps)==0) mxps <- vector('list',length=id.train.no)
    if(length(mxps$prior)!=0&id.train.no==1) mxps <- list(mxps)
    mixed.prior <- rep(0,id.train.no)
    mixed.prior.df <- rep(0,id.train.no)
    mixed.prior.scale <- list()
    for(i in 1:NCOL(id.train)){
        if(length(mxps[[i]]$prior)==0) mixed.prior[i] <- 1
        else mixed.prior[i] <- mxps[[i]]$prior
        if(length(mxps[[i]]$df)==0) mixed.prior.df[i] <- 3
        else mixed.prior.df[i] <- mxps[[i]]$df
        ## cmxbart.cpp uses an Inverse-Gamma distribution instead of a Scaled-Inverse ChiSquare so the parameters need to be converted.
    }
    mixed.prior.scale <- unlist(mixed.prior.scale)
    if(!transposed) {
        temp = BART3::bartModelMatrix(x.train, numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        ## if(length(x.test)>0)
        ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if(length(x.test)>0) {
            x.test = BART3::bartModelMatrix(x.test)
            x.test = t(x.test[ , temp$rm.const])
        }
        rm.const <- temp$rm.const
        grp <- temp$grp
        rm(temp)
    }
    else {
        rm.const <- NULL
        grp <- NULL
    }
    if(n!=ncol(x.train))
        stop('The length of y.train and the number of rows in x.train must be identical')

    p = nrow(x.train)
    np = ncol(x.test)
    if(length(rho)==0) rho=p
    if(length(rm.const)==0) rm.const <- 1:p
    if(length(grp)==0) grp <- 1:p

    check <- unique(sort(y.train))
    if(length(check)==2) {
        if(!all(check==0:1))
            stop('Binary y.train must be coded as 0 and 1')
        if(type=='wbart')
            stop("The outcome is binary so set type to 'pbart'")
            ##stop("The outcome is binary so set type to 'pbart' or 'lbart'")
    }

    ## check <- c('wbart', 'pbart', 'lbart')

    ## if(!(type %in% check))
    ##     stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    if(length(offset)==0) {
        if(type=='wbart') offset=mean(y.train)
        if(type=='pbart') offset=stats::qnorm(mean(y.train))
        ##else if(type=='lbart') offset=qlogis(offset)
    }
    if(type=='wbart') {
        if(is.na(tau.num)) {
            tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
        }
        else {
            tau=tau.num/(k*sqrt(ntree))
        }
        y.train = y.train-offset
        if(is.na(lambda)) {
            if(is.na(sigest)) {
                if(p < n) {
                    tx.train <- t(x.train)
                    tmp.lmer <- lme4::lmer(y.train~tx.train+(1|id.train))
                    sdu.est <- sqrt(lme4::VarCorr(tmp.lmer)[1]$id[1])
                    sigest <- summary(tmp.lmer)$sigma
                } else sigest = stats::sd(y.train)
            }
            qchi = stats::qchisq(sigquant, sigdf)
            lambda = (sigest^2)*qchi/sigdf #lambda parameter for sigma prior
        }
        else {
            sigest=sqrt(lambda)
        }
    }
    else {
         tau <- 3/(k*sqrt(ntree))
         if(p<n) {
             tx.train <- t(x.train)
             tmp.glmer <- lme4::glmer(y.train~tx.train+(1|id.train),family=binomial("probit"))
             sdu.est <- sqrt(lme4::VarCorr(tmp.glmer)[1]$id[1])
         }
    }
    for(i in 1:NCOL(id.train)) {
        if(mixed.prior[i]==1){
            if(length(mxps[[i]]$scale)==0)
                mixed.prior.scale[[i]] <- sdu.est^2/stats::qt(.5*1.95,df=mixed.prior.df[i])
            else
                mixed.prior.scale[[i]] <- mxps[[i]]$scale
            if(length(mxps[[i]]$scale!=z.cols[i]))
               mixed.prior.scale[[i]] <- rep(mixed.prior.scale[[i]][1],z.cols[i])
        }
        else {
            if(length(mxps[[i]]$scale)==0)
                mixed.prior.scale[[i]] <- sdu.est^2*stats::qchisq(0.95,mixed.prior.df[i])/mixed.prior.df[i]
            else
                mixed.prior.scale[[i]] <- mxps[[i]]$scale
            if(length(mxps[[i]]$scale)!=z.cols[i])
                mixed.prior.scale[[i]] <- rep(mixed.prior.scale[[i]][1],z.cols[i])
            ## Convert from scaled-inverse Chi-square parameterization to inverse-gamma parameterization
            mixed.prior.df[i] <- .5*mixed.prior.df[i]
            mixed.prior.scale[[i]] <- mixed.prior.df[i]*mixed.prior.scale[[i]]
        }
    }
    if(type=='pbart')
        tau=qnorm(1-0.5*tau)/(k*sqrt(ntree))
    
    #else if(type=='lbart')
    ##     tau=qlogis(1-0.5*tau)/(k*sqrt(ntree))

    ## hot deck missing imputation
    ## must be conducted here since it would
    ## cause trouble with multi-threading on the C++ side

    check=(np>0 && np==n)

    for(i in 1:n)
        for(j in 1:p) {
            if(check) check=((is.na(x.train[j, i]) && is.na(x.test[j, i])) ||
                             (!is.na(x.train[j, i]) && !is.na(x.test[j, i]) &&
                              x.train[j, i]==x.test[j, i]))

            while(is.na(x.train[j, i])) {
                h=sample.int(n, 1)
                x.train[j, i]=x.train[j, h]
            }
        }

    if(check) x.test=x.train
    else if(np>0) {
        for(i in 1:np)
            for(j in 1:p)
                while(is.na(x.test[j, i])) {
                    h=sample.int(np, 1)
                    x.test[j, i]=x.test[j, h]
                }
    }
    
    if(.Platform$OS.type!='unix') hostname <- FALSE
    else if(hostname)
        hostname <- system('hostname', intern=TRUE)

    ptm <- proc.time()
    
    res = .Call("cmxbart",
                ntype, ##as.integer(factor(type, levels=check))-1,
                n,  #number of observations in training data
                p,  #dimension of x
                np, #number of observations in test data
                x.train,   #pxn training data x
                y.train,   #pxn training data x
                x.test,    #p*np test data x
                c.train,
                z.train,
                z.cols,
                n.train,
                id.train.no,
                L,
                mixed.prior,
                mixed.prior.df,
                mixed.prior.scale,
                ntree,
                numcut,
                ndpost*keepevery,
                nskip,
                keepevery,
                power,
                base,
                offset,
                tau,
                sigdf,
                lambda,
                sigest,
                ##w,
                sparse,
                theta,
                omega,
                grp,
                a,
                b,
                rho,
                augment,
                printevery,
                xinfo
                )

    res$proc.time <- proc.time()-ptm
    res$hostname <- hostname

    tmp.varcov <- list()
    tmp.varcov.mean <- list()
    tmp.corr <- list()
    tmp.corr.mean <- list()
    for(i in 1:id.train.no){
        tmp.varcov[[i]] <- aperm(array(unlist(res$re.varcov[[i]]),dim=c(z.cols[i],z.cols[i],ndpost)),c(3,1,2))
        tmp.varcov.mean[[i]] <- apply(tmp.varcov[[i]],2:3,mean)
        tmp.corr[[i]] <- aperm(array(apply(tmp.varcov[[i]],1,stats::cov2cor),dim=c(z.cols[i],z.cols[i],ndpost)),c(3,1,2))
        tmp.corr.mean[[i]] <- apply(tmp.corr[[i]],2:3,mean)
    }
    res$re.varcov <- tmp.varcov
    res$re.varcov.mean <- tmp.varcov.mean
    res$re.corr <- tmp.corr
    res$re.corr.mean <- tmp.corr.mean
    
    tmp.re <- list()
    tmp.re.mean <- list()
    for(i in 1:id.train.no){
        tmp.re[[i]] <- aperm(array(unlist(res$re.train[[i]]),dim=c(L[i],z.cols[i],ndpost)),c(3,1,2))
        tmp.re.mean[[i]] <- apply(tmp.re[[i]],2:3,mean)
    }
    res$re.train <- tmp.re
    res$re.train.mean <- tmp.re.mean

    
    res$fhat.train.mean <- apply(res$fhat.train, 2, mean)
    if(type=='pbart') {
        res$prob.train = stats::pnorm(res$fhat.train)
        ##else if(type=='lbart') res$prob.train = plogis(res$yhat.train)
        res$prob.train.mean <- apply(res$prob.train, 2, mean)
    }
    if(np>0) {
        res$fhat.test.mean <- apply(res$fhat.test, 2, mean)
        if(type=='pbart') {
            res$prob.test = stats::pnorm(res$fhat.test)
            res$prob.test.mean <- apply(res$prob.test, 2, mean)
        }
    }
    ## L.mu.post <- matrix(0,nrow=ndpost,ncol=n)
    ## cpo <- double(n)
    ## if(z.cols>1){
    ##     u.train.expanded <- res$u.train[,rep(1:L,times=n.train),]
    ##     for(i in 1:n){
    ##         for(s in 1:ndpost){
    ##             L.mu.post[s,i] <- res$yhat.train[s,i]+sum(z.train[i,1:z.cols]*u.train.expanded[s,i,1:z.cols])
    ##         }
    ##         cpo[i] <- ndpost/sum(1/dnorm(y.train[i],mean=L.mu.post[,i],sd=res$sigma))
    ##     }
    ## }
    ## else{
    ##     u.train.expanded <- res$u.train[,rep(1:L,times=n.train),1]
    ##     for(i in 1:n){
    ##         for(s in 1:ndpost){
    ##             L.mu.post[s,i] <- res$yhat.train[s,i]+z.train[i,1]*u.train.expanded[s,i]
    ##         }
    ##         cpo[i] <- ndpost/sum(1/dnorm(y.train[i],mean=L.mu.post[,i],sd=res$sigma))
    ##     }
    ## }
    ## res$lpml <- mean(log(cpo))
         
    res$offset = offset
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    res$varcount.mean <- apply(res$varcount, 2, mean)
    res$varprob.mean <- apply(res$varprob, 2, mean)
    res$rm.const <- rm.const
    res$sigest <- sigest
    res$z.cols <- z.cols
    if(type=='wbart'){
        if(!sparse) 
            res <- res[c('fhat.train','fhat.train.mean','fhat.test','fhat.test.mean','sigma',
                         're.train','re.train.mean','re.varcov','re.varcov.mean',
                         're.corr','re.corr.mean','z.cols','offset','sigest',
                         'varprob','varprob.mean','varcount','varcount.mean',
                         'rm.const','hostname','proc.time')]
        else
            res <- res[c('fhat.train','fhat.train.mean','fhat.test','fhat.test.mean','sigma',
                         're.train','re.train.mean','re.varcov','re.varcov.mean',
                         're.corr','re.corr.mean','z.cols','offset','sigest',
                         'varprob','varprob.mean','varcount','varcount.mean','theta.train',
                         'rm.const','hostname','proc.time')]
    }
    else {
        if(!sparse)
            res <- res[c('fhat.train','fhat.train.mean','fhat.test','fhat.test.mean',
                         'prob.train','prob.train.mean','prob.test','prob.test.mean',
                         're.train','re.train.mean','re.varcov','re.varcov.mean',
                         're.corr','re.corr.mean','z.cols','offset','sigest',
                         'varprob','varprob.mean','varcount','varcount.mean',
                         'rm.const','hostname','proc.time')]
        else
            res <- res[c('fhat.train','fhat.train.mean','fhat.test','fhat.test.mean',
                         'prob.train','prob.train.mean','prob.test','prob.test.mean',
                         're.train','re.train.mean','re.varcov','re.varcov.mean',
                         're.corr','re.corr.mean','z.cols','offset','sigest',
                         'varprob','varprob.mean','varcount','varcount.mean','theta.train',
                         'rm.const','hostname','proc.time')]
        }
        attr(res, 'class') <- 'mxbart'
        return(res)
}
