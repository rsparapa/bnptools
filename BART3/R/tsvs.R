
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2023 Robert McCulloch and Rodney Sparapani
## tsvs.R

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

tsvs <- function(
                 x.train=matrix(0,0,0), y.train=NULL,
                 ##tsvs args
                 T=20, a.=1, b.=0.5, C=0.5,
                 rds.file='tsvs.rds',
                 pdf.file='tsvs.pdf',
                 type='wbart',
                 ntype=as.integer(
                     factor(type,
                            levels=c('wbart', 'pbart', 'lbart'))),
                 sparse=TRUE, theta=0, omega=1,
                 a=0.5, b=1, augment=FALSE, rho=0, grp=NULL,
                 varprob=NULL,
                 xinfo=matrix(0,0,0), usequants=FALSE,
                 rm.const=TRUE,
                 sigest=NA, sigdf=3, sigquant=0.90,
                 k=2, power=2, base=0.95,
                 impute.mult=NULL, impute.prob=NULL, impute.miss=NULL,
                 lambda=NA, tau.num=c(NA, 3, 6)[ntype],
                 offset=NULL, w=rep(1, length(y.train)),
                 ntree=10L, numcut=100L,
                 ndpost=1000L, nskip=100L,
                 keepevery=c(1L, 10L, 10L)[ntype],
                 printevery=100L, transposed=FALSE,
                 probs=c(0.025, 0.975),
                 verbose = 1L
                 ## shards=1L, weight=rep(NA, shards),
                 ## meta = FALSE
                 )
{

    if(T==0) return(T)

    if(transposed) 
        stop('tsvs must be run with x.train untransposed')

    if(is.na(ntype))
        stop("type argument must be set to either 'wbart', 'pbart' or 'lbart'")

    if(length(y.train)==0)
        stop('Supply a non-zero length y.train vector')
    if(length(x.train)==0)
        stop('Supply a non-zero length x.train matrix')

    check <- unique(sort(y.train))

    if(length(check)==2) {
        if(!all(check==0:1))
            stop('Binary y.train must be coded as 0 and 1')
        if(type=='wbart')
            stop("The outcome is binary so set type to 'pbart' or 'lbart'")
    }

    ## if(.Platform$OS.type!='unix')
    ##     stop('parallel::mcparallel/mccollect do not exist on windows')

    ## RNGkind("L'Ecuyer-CMRG")
    ## parallel::mc.reset.stream()

    if(length(impute.mult)==1)
        stop("The number of multinomial columns must be greater than 1\nConvert a binary into two columns")

    temp = bartModelMatrix(x.train, numcut, usequants=usequants,
                           xinfo=xinfo, rm.const=rm.const)
    x.train = temp$X
    numcut = temp$numcut
    xinfo = temp$xinfo
    if(length(grp)==0) grp <- temp$grp
    rm(temp)
    
    P=ncol(x.train)
    if(length(grp)==0) grp <- rep(1, P)

    grp.=0
    h=1
    j=0
    for(i in grp) {
        grp.[h]=i
        if(j==0) {
            h=h+1
            j=i-1
        } else { j=j-1 }
        ##return(list(grp=grp, grp.=grp.))
    }
    dummy=matrix(0, nrow=2, ncol=P)
    h=1
    l=1
    for(i in 1:length(grp.)) {
        for(j in 1:grp.[i]) {
            dummy[1, h]=l
            dummy[2, h]=l+grp.[i]-1
            h=h+1
        }
        l=h
    }
    dimnames(dummy)[[2]]=dimnames(x.train)[[2]]
    ##return(dummy)

    Names=dimnames(x.train)[[2]] 
    A=matrix(a., nrow=T, ncol=P)
    B=matrix(b., nrow=T, ncol=P)
    S=matrix(0,  nrow=T, ncol=P)
    dimnames(S)[[2]]=Names
    prob=matrix(nrow=T, ncol=P)
    dimnames(prob)[[2]]=Names
    reward=matrix(0, nrow=T, ncol=P)
    dimnames(reward)[[2]]=Names
    vimp=matrix(nrow=T, ncol=P)
    dimnames(vimp)[[2]]=Names
    varcount=matrix(0, nrow=T, ncol=P)
    dimnames(varcount)[[2]]=Names

    for(i in 1:T) {
        set.seed(i)
        print(paste('Step:', i))
        if(i>1) {
            for(j in 1:P) {
                A[i, j]=A[i-1, j]
                B[i, j]=B[i-1, j]
            }
        }
        prob[i, ]=rbeta(P, A[i, ], B[i, ])
        S[i, which(prob[i, ]>=C)]=1
        
        j=sum(S[i, ])
        if(j==0) S[i, sample.int(P, 2)]=1
        else if(j==1) S[i, sample(which(S[i, ]==0), 1)]=1
        
        for(j in 1:P)
            if(S[i, j]==1) S[i, dummy[1, j]:dummy[2, j] ]=1

        set.seed(T+i*T)

        pick=(S[i, ]==1)
        h=0
        numcut.=0
        grp.=0
        for(j in 1:P)
            if(pick[j]) {
                h=h+1
                if(h==1) xinfo.=rbind(xinfo[j, ])
                else xinfo.=rbind(xinfo., xinfo[j, ])
                numcut.[h]=numcut[j]
                grp.[h]=grp[j]
            }
        ##dimnames(xinfo.)[[2]]=dimnames(xinfo)[[2]]
        
        x.train.=cbind(x.train[ , pick])
        dimnames(x.train.)[[2]]=Names[pick]
        
        post=gbart(x.train=t(x.train.), y.train=y.train,
                   type=type, ntype=ntype,
                   sparse=sparse, theta=theta, omega=omega,
                   a=a, b=b, augment=augment, rho=rho, grp=grp.,
                   varprob=varprob,
                   xinfo=xinfo., usequants=usequants,
                   rm.const=rm.const,
                   sigest=sigest, sigdf=sigdf, sigquant=sigquant,
                   k=k, power=power, base=base,
                   impute.mult=impute.mult, impute.prob=impute.prob,
                   impute.miss=impute.miss,
                   lambda=lambda, tau.num=tau.num,
                   offset=offset, w=w,
                   ntree=ntree, numcut=numcut.,
                   ndpost=ndpost, nskip=nskip,
                   keepevery=keepevery, printevery=printevery,
                   probs=probs,
                   ##mc.cores = mc.cores,
                   ##nice = nice,
                   TSVS=TRUE, verbose = verbose,
                   ##shards=shards, weight=weight, ##meta = meta,
                   transposed=TRUE
                   )
        ##return(post)

        names.=dimnames(post$varcount)[[2]]
        M=nrow(post$varcount)
        for(j in 1:P) {
            if(S[i, j]==1) {
                h=which(Names[j]==names.)
                l=post$varcount[M, h]
                if(l>0) {
                    varcount[i, j]=l
                    reward[i, j]=1
                }
                A[i, j]=A[i, j]+reward[i, j]
                B[i, j]=B[i, j]+1-reward[i, j]
            } else {
                B[i, j]=B[i, j]+1
            }
            vimp[i, j]=A[i, j]/(A[i, j]+B[i, j])
        }
        if(length(warnings())>0) print(warnings())
        res=list(step=i, vimp=vimp, S=S, a=A, b=B, reward=reward,
                 prob=prob, varcount=varcount)
        saveRDS(res, rds.file)

        pdf(file=pdf.file)
        par(mfrow=c(2, 1))
        plot(1:i, vimp[1:i, 1], type='n', ylim=c(0, 1), xlim=c(0, T),
             xlab='Steps', ylab='Inclusion Probability')
        abline(h=0:1, v=c(0, T))
        abline(h=0.5, col=8, lty=3)
        for(j in 1:P) 
            if(vimp[i, j]>0.5) {
                if(i==1) points(i, prob[i, j], col=j)
                else lines(1:i, prob[1:i, j], col=j)
                h=sample(1:i, 1)
                text(h, prob[h, j], Names[j], col=j, pos=1)
            }
        plot(1:i, vimp[1:i, 1], type='n', ylim=c(0, 1), xlim=c(0, T),
             xlab='Steps', ylab='VIMP Probability')
        abline(h=0:1, v=c(0, T))
        abline(h=0.5, col=8, lty=3)
        for(j in 1:P) 
            if(vimp[i, j]>0.5) {
                if(i==1) points(i, vimp [i, j], col=j)
                else lines(1:i, vimp [1:i, j], col=j)
                h=sample(1:i, 1)
                text(h, vimp [h, j], Names[j], col=j, pos=1)
            }
        par(mfrow=c(1, 1))
        dev.off()
    }

    return(res)
}
