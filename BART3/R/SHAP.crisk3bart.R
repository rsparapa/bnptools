
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani

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

## Shapley additive explanation (SHAP) partial dependence function
SHAP.crisk3bart=function(object,       ## object returned from BART
                         x.train,      ## x.train to estimate coverage
                         x.test,       ## settings of x.test: only x.test[ , S]
                                       ## are used but they must all be given
                         S,            ## indices of subset
                         type='pbart', ## type of probability model
                         probs=c(0.025, 0.975))
{
    for(v in S)
        if(any(is.na(x.test[ , v])))
            stop(paste0('x.test column with missing values:', v))

    P = ncol(x.train)

    if(!all(S %in% 1:P))
        stop('some elements of S are not in x.train')

    if(P!=ncol(x.test))
        stop('the number of columns in x.train and x.test are not the same')

    if(P!=length(object$treedraws$cutpoints))
        stop('the number of columns in x.train and length of cutpoints are not the same')

    if(P!=length(object$treedraws2$cutpoints))
        stop('the number of columns in x.train and length of cutpoints2 are not the same')

    if(P!=length(object$treedraws3$cutpoints))
        stop('the number of columns in x.train and length of cutpoints3 are not the same')

    Trees =read.trees(object$treedraws,  x.train)
    Trees2=read.trees(object$treedraws2, x.train)
    Trees3=read.trees(object$treedraws3, x.train)

    EXPVALUE = function(trees)
    {
        H = nrow(x.test)
        M = dim(trees)[1]
        T = dim(trees)[2]
        G = function(n) ## node
        {
            if(trees[i, j, n, 1]==2) return(trees[i, j, n, 4]) ## a leaf
            else { ## a branch
                v=trees[i, j, n, 2]
                c=trees[i, j, n, 3]
                n=2*n
                m=n+1
                if(v %in% S) {
                    if(x.test[h, v]<c) return(G(n))
                    else return(G(m))
                } else {
                    a=trees[i, j, n, 5]
                    b=trees[i, j, m, 5]
                    return((a*G(n)+b*G(m))/(a+b))
                }
            }
        }
        A = matrix(nrow=M, ncol=H)
        B = matrix(nrow=M, ncol=T)
        for(h in 1:H) { ## settings
            for(i in 1:M) ## samples
                for(j in 1:T) ## trees
                    B[i, j]=G(1)
            A[ , h]=apply(B, 1, sum)
        }
        return(A)
    }

    pred=list()

    pred$times <- object$times
    K <- object$K
    pred$K = K

    pred$offset  <- object$offset
    pred$offset2 <- object$offset2
    pred$offset3 <- object$offset3

    pred$yhat.test  <- object$offset +EXPVALUE(Trees)
    pred$yhat.test2 <- object$offset2+EXPVALUE(Trees2)
    pred$yhat.test3 <- object$offset3+EXPVALUE(Trees3)

    if(type=='pbart') {
        pred$prob.test  <- pnorm(pred$yhat.test)
        pred$prob.test2 <- pnorm(pred$yhat.test2)
        pred$prob.test3 <- pnorm(pred$yhat.test3)
    }
    else if(type=='lbart') {
        pred$prob.test  <- plogis(pred$yhat.test)
        pred$prob.test2 <- plogis(pred$yhat.test2)
        pred$prob.test3 <- plogis(pred$yhat.test3)
    }

    pred$surv.test <- (1-pred$prob.test)*(1-pred$prob.test2)*
        (1-pred$prob.test3)
    pred$cif.test  <- pred$prob.test
    pred$cif.test2 <- (1-pred$prob.test)*pred$prob.test2
    pred$cif.test3 <- (1-pred$prob.test)*(1-pred$prob.test2)*pred$prob.test3

    H <- nrow(x.test)/K ## the number of different settings

    for(h in 1:H)
        for(j in 2:K) {
            l <- K*(h-1)+j

            pred$cif.test[ , l] <- pred$cif.test[ , l-1]+
                pred$surv.test[ , l-1]*pred$cif.test[ , l]
            pred$cif.test2[ , l] <- pred$cif.test2[ , l-1]+
                pred$surv.test[ , l-1]*pred$cif.test2[ , l]
            pred$cif.test3[ , l] <- pred$cif.test3[ , l-1]+
                pred$surv.test[ , l-1]*pred$cif.test3[ , l]
            pred$surv.test[ , l] <- pred$surv.test[ , l-1]*
                pred$surv.test[ , l]
        }

    pred$cif.test.mean   <- apply(pred$cif.test, 2, mean)
    pred$cif.test.lower  <- apply(pred$cif.test, 2, quantile, probs=probs[1])
    pred$cif.test.upper  <- apply(pred$cif.test, 2, quantile, probs=probs[2])
    pred$cif.test2.mean  <- apply(pred$cif.test2, 2, mean)
    pred$cif.test2.lower <- apply(pred$cif.test2, 2, quantile, probs=probs[1])
    pred$cif.test2.upper <- apply(pred$cif.test2, 2, quantile, probs=probs[2])
    pred$cif.test3.mean  <- apply(pred$cif.test3, 2, mean)
    pred$cif.test3.lower <- apply(pred$cif.test3, 2, quantile, probs=probs[1])
    pred$cif.test3.upper <- apply(pred$cif.test3, 2, quantile, probs=probs[2])
    pred$surv.test.mean  <- apply(pred$surv.test, 2, mean)
    pred$surv.test.lower <- apply(pred$surv.test, 2, quantile, probs=probs[1])
    pred$surv.test.upper <- apply(pred$surv.test, 2, quantile, probs=probs[2])
    pred$prob.test.mean  <- apply(pred$prob.test, 2, mean)
    pred$prob.test2.mean <- apply(pred$prob.test2, 2, mean)
    pred$prob.test3.mean <- apply(pred$prob.test3, 2, mean)

    attr(pred, 'class') <- 'crisk3bart'

    return(pred)
}
