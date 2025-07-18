
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2018-2025 Robert McCulloch and Rodney Sparapani
## gbart.R

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

gbart=function(
               x.train, y.train, x.test=matrix(0,0,0),
               ##z.train=NULL,
               type='wbart',
               ntype=as.integer(
                   factor(type, levels=c('wbart', 'pbart', 'lbart'))),
               ##treeinit=FALSE, trees=NULL,
               sparse=FALSE, theta=0, omega=1,
               a=0.5, b=1, augment=FALSE, rho=0, grp=NULL,
               varprob=NULL,
               xinfo=matrix(0,0,0), usequants=FALSE,
               rm.const=TRUE,
               sigest=NA, sigdf=3, sigquant=0.90,
               k=2, power=2, base=0.95,
               impute.mult=NULL, impute.prob=NULL,
               impute.miss=NULL,
               lambda=NA, tau.num=c(NA, 3, 6)[ntype],
               offset=NULL, w=rep(1, length(y.train)),
               ntree=c(200L, 50L, 50L)[ntype], numcut=100L,
               ndpost=1000L, nskip=100L,
               keepevery=c(1L, 10L, 10L)[ntype],
               printevery=100L, transposed=FALSE,
               probs=c(0.025, 0.975),
               mc.cores = 1L, ## mc.gbart only
               nice = 19L,    ## mc.gbart only
               seed = 99L,    ## mc.gbart only
               meta = FALSE,  ## mc.gbart only
               TSVS=FALSE,    ## gbart only
               verbose = 1L,
               shards = 1L, weight=rep(NA, shards)
               )
{

   if(is.na(ntype))
        stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    n = length(y.train)

    if(!transposed) {
        temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                               xinfo=xinfo, rm.const=rm.const)
        x.train = t(temp$X)
        numcut = temp$numcut
        xinfo = temp$xinfo
        ## if(length(x.test)>0)
        ##     x.test = t(bartModelMatrix(x.test[ , temp$rm.const]))
        if(length(x.test)>0) {
            x.test = t(bartModelMatrix(x.test))
            if(class(rm.const)[1]=='logical' && rm.const)
                x.test = rbind(x.test[temp$rm.const, ])
        }
        rm.const <- temp$rm.const
        if(length(grp)==0) grp <- temp$grp
        rm(temp)
    } else {
        rm.const <- NULL
        ##grp <- NULL
    }

    if(n!=ncol(x.train))
        stop('The length of y.train and the number of rows in x.train must be identical')

    p = nrow(x.train)
    np = ncol(x.test)

    if(np>0 && p!=nrow(x.test))
        stop('The number of columns in x.train and x.test must be identical')

    if(length(grp)==0) grp=rep(1, p)
    if(rho==0) rho=sum(1/grp)
    ##if(length(rho)==0) rho=p
    if(length(rm.const)==0) rm.const <- 1:p
    ##if(length(grp)==0) grp <- 1:p
    if(length(varprob)==0) varprob=rep(1, p)
    varprob=varprob/sum(varprob)

    check = length(impute.mult)
    if(check==1)
        stop("The number of multinomial columns must be greater than 1\nConvert a binary into two columns")
    if(check>1) {
        impute.flag=TRUE
        if(length(impute.miss)==0) {
            impute.miss=integer(n)
            for(j in 1:check) {
                i=impute.mult[j]
                impute.miss = pmax(impute.miss, is.na(x.train[i, ]))
            }
        }
        if(length(impute.prob)==0) {
            impute.prob=double(check)
            for(j in 1:check) {
                i=impute.mult[j]
                impute.prob[j]=sum(x.train[i, ]==1, na.rm = TRUE)
            }
            impute.prob=impute.prob/sum(impute.prob)
            impute.prob=matrix(impute.prob,
                               nrow=n, ncol=check, byrow=TRUE)
        }
        impute.mult=as.integer(impute.mult-1) ## convert from R index to C/C++
    } else {
        impute.flag=FALSE
        impute.mult=integer(0) ## integer vector of column indicators for missing covariates
        impute.miss=integer(0) ## integer vector of row indicators for missing values
        impute.prob=matrix(nrow=0, ncol=0)
    }

    check <- unique(sort(y.train))

    if(length(check)==2) {
        if(!all(check==0:1))
            stop('Binary y.train must be coded as 0 and 1')
        if(type=='wbart')
            stop("The outcome is binary so set type to 'pbart' or 'lbart'")
    }

    ## check <- c('wbart', 'pbart', 'lbart')

    ## if(!(type %in% check))
    ##     stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    if(length(offset)==0) {
        offset=mean(y.train)
        if(type=='pbart') offset=qnorm(offset)
        else if(type=='lbart') offset=qlogis(offset)
    }

    ##if(length(z.train)==0) z.train=y.train
    ## trees=''
    ## treeinit = FALSE
    ## treeinit = (treeinit && length(trees)>0)
    ## if(treeinit) {
    ##    if(trees=="") stop('Tree string is empty')
    ##    ## more likely a previous BART fit than a Random Forest
    ##    ## rf=randomForest(t(x.train), y.train, ntree=ntree,
    ##    ##                 maxnodes=4, forest=TRUE)
    ##    ## sigest=sd(rf$predicted-y.train)
    ##    ## trees=read.forest(rf, xinfo=xinfo)
    ## }

    if(type=='wbart') {
        y.train = y.train-offset
        ##z.train = z.train-offset
        if(!is.na(sigest) && !is.na(lambda) && lambda==0) {
            ##no op: sigma is fixed and known at given sigest value
        }
        else if(is.na(lambda)) {
            if(is.na(sigest)) {
                if(p < n) 
                    sigest =
                        summary(lm(y.train~.,
                                   data.frame(t(x.train),y.train)))$sigma
                if(is.na(sigest)) sigest = sd(y.train)
            }
            qchi = qchisq(1-sigquant, sigdf)
            lambda = (sigest^2)*qchi/sigdf #lambda parameter for sigma prior
        } else {
            sigest=sqrt(lambda)
        }
        if(is.na(tau.num)) {
            tau=(max(y.train)-min(y.train))/(2*k*sqrt(ntree))
        } else {
            tau=tau.num/(k*sqrt(ntree))
        }
    } else {
        lambda=1
        sigest=1
        tau=tau.num/(k*sqrt(ntree))
        ## tau=1-tau.interval
        ## if(type=='pbart')
        ##     tau=qnorm(1-0.5*tau)/(k*sqrt(ntree))
        ## else if(type=='lbart')
        ##     tau=qlogis(1-0.5*tau)/(k*sqrt(ntree))
    }
    ## k can now be re-used

    ## hot deck missing imputation
    ## must be conducted here since it would
    ## cause trouble with multi-threading on the C++ side
    ## ignore columns for multinomial imputation in training

    check=(np>0 && np==n)

    if(check)
        for(i in 1:n) {
            for(j in 1:p) {
                check=((is.na(x.train[j, i]) &&
                        is.na(x.test[j, i])) ||
                       (!is.na(x.train[j, i]) &&
                        !is.na(x.test[j, i]) &&
                        x.train[j, i]==x.test[j, i]))
                if(!check) break
            }
            if(!check) break
        }

    for(i in 1:n)
        for(j in 1:p)
            if(impute.flag && !(j %in% impute.mult)) {
                k = is.na(x.train[ , i])
                if(impute.flag) k[impute.mult]=FALSE
                while(is.na(x.train[j, i])) {
                    h=sample.int(n, 1)
                    x.train[which(k), i]=x.train[which(k), h]
                }
            }

    if(check && !impute.flag) { x.test=x.train
    } else if(np>0) {
        for(i in 1:np)
            for(j in 1:p) {
                k = is.na(x.test[ , i])
                while(is.na(x.test[j, i])) {
                    h=sample.int(np, 1)
                    x.test[which(k), i]=x.test[which(k), h]
                }
            }
    }

    ## if(hotdeck) ## warnings are suppressed with mc.gbart anyways
    ##     warning('missing elements of x imputed with hot decking')

    ## if(.Platform$OS.type!='unix') hostname <- FALSE
    ## else if(hostname)
    ##     hostname <- system('hostname', intern=TRUE)

    ptm <- proc.time()

    ## for debugging
    ## return(list(ntype=ntype, ##as.integer(factor(type, levels=check))-1,
    ##             n=n,  #number of observations in training data
    ##             p=p,  #dimension of x
    ##             np=np, #number of observations in test data
    ##             x.train=x.train,   #pxn training data x
    ##             y.train=y.train,   #nx1 training data y
    ##             x.test=x.test,    #pxnp test data x
    ##             ##z.train,
    ##             ntree=ntree,
    ##             numcut=numcut,
    ##             ndpost=ndpost*keepevery,
    ##             nskip=nskip,
    ##             keepevery=keepevery,
    ##             power=power,
    ##             base=base,
    ##             offset=offset,
    ##             tau=tau,
    ##             sigdf=sigdf,
    ##             lambda=lambda,
    ##             sigest=sigest,
    ##             w=w,
    ##             ##as.integer(treeinit),
    ##             ##trees,
    ##             sparse=sparse,
    ##             theta=theta,
    ##             omega=omega,
    ##             grp=grp,
    ##             a=a,
    ##             b=b,
    ##             rho=rho,
    ##             augment=augment,
    ##             varprob=varprob,
    ##             printevery=printevery,
    ##             xinfo=xinfo,
    ##             shards=shards,
    ##             verbose=verbose,
    ##             impute.mult=impute.mult,
    ##             impute.miss=impute.miss,
    ##             impute.prob=impute.prob
    ##             ))

    res = .Call("cgbart",
                ntype, ##as.integer(factor(type, levels=check))-1,
                n,  #number of observations in training data
                p,  #dimension of x
                np, #number of observations in test data
                x.train,   #pxn training data x
                y.train,   #nx1 training data y
                x.test,    #pxnp test data x
                ##z.train,
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
                w,
                ##as.integer(treeinit),
                ##trees,
                sparse,
                theta,
                omega,
                grp,
                a,
                b,
                rho,
                augment,
                varprob,
                printevery,
                xinfo,
                shards,
                verbose,
                impute.mult,
                impute.miss,
                impute.prob
                )

    res$proc.time <- proc.time()-ptm
    dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    res$varcount.mean <- apply(res$varcount, 2, mean)
    dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    res$varprob.mean <- apply(res$varprob, 2, mean)
if(TSVS) return(res)
##    res$hostname <- hostname

    type1.sigest=(type=='wbart' && nskip>0)
    if(type1.sigest && !is.na(sigest) && !is.na(lambda) && lambda==0)
        type1.sigest=FALSE

    Y=t(matrix(y.train, nrow=n, ncol=ndpost))

    if(type=='wbart') {
        res$yhat.train.mean <- apply(res$yhat.train, 2, mean)
        res$yhat.train.lower <- apply(res$yhat.train, 2, quantile,
                                      probs=min(probs))
        res$yhat.train.upper <- apply(res$yhat.train, 2, quantile,
                                      probs=max(probs))
        res$sigest <- sigest
        if(type1.sigest) {
            res$sigma.=res$sigma[-(1:nskip)]
            SD=matrix(res$sigma., nrow=ndpost, ncol=n)
            ##CPO=1/apply(1/dnorm(Y, res$yhat.train, SD), 2, mean)
            log.pdf=dnorm(Y, res$yhat.train, SD, TRUE)
            res$sigma.mean=mean(SD[ , 1])
        } else {
            res$sigma.mean = sigest
        }
    }
    else {
        if(type=='pbart') res$prob.train = pnorm(res$yhat.train)
        else if(type=='lbart') res$prob.train = plogis(res$yhat.train)

        ##CPO=1/apply(1/dbinom(Y, 1, res$prob.train), 2, mean)
        log.pdf=dbinom(Y, 1, res$prob.train, TRUE)

        res$prob.train.mean <- apply(res$prob.train, 2, mean)
    }

    if(type1.sigest | type!='wbart') {
        min.log.pdf=t(matrix(apply(log.pdf, 2, min), nrow=n, ncol=ndpost))
        log.CPO=log(ndpost)+min.log.pdf[1, ]-
            log(apply(exp(min.log.pdf-log.pdf), 2, sum))
        res$LPML=sum(log.CPO)
        ##res$CPO=exp(log.CPO)
        ##res$LPML=sum(log(CPO))
    }

    res$offset = offset
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
    res$x.train = t(x.train)

    keeptestfits <- (np>0)
    ##if(length(keeptestfits)==0) keeptestfits <- (np>0)

    if(keeptestfits) {
        res$x.test = t(x.test)
        if(type=='wbart') {
            res$yhat.test.mean <- apply(res$yhat.test, 2, mean)
            res$yhat.test.lower <- apply(res$yhat.test, 2, quantile,
                                         probs=min(probs))
            res$yhat.test.upper <- apply(res$yhat.test, 2, quantile,
                                         probs=max(probs))
        }
        else {
            if(type=='pbart') res$prob.test = pnorm(res$yhat.test)
            else if(type=='lbart') res$prob.test = plogis(res$yhat.test)

            res$prob.test.mean <- apply(res$prob.test, 2, mean)
            res$prob.test.lower <- apply(res$prob.test, 2, quantile,
                                         probs=min(probs))
            res$prob.test.upper <- apply(res$prob.test, 2, quantile,
                                         probs=max(probs))
        }
    }

    ## dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
    ## res$varcount.mean <- apply(res$varcount, 2, mean)
    ## dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
    ## res$varprob.mean <- apply(res$varprob, 2, mean)
    if(impute.flag) res$impute.miss = impute.miss
    res$rm.const <- rm.const
    res$ndpost = ndpost
##    res$trees = trees
    res$grp = grp
    res$rho = rho
    res$chains = mc.cores
    ##if(treeinit) res$trees = trees
    attr(res, 'class') <- type

    return(res)
}
