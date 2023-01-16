
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

## Hot deck (HD) SHAP additive explanation function
HDSHAP.wbart=function(object,  ## object returned from BART
                      x.train, ## x.train to estimate coverage
                      x.test,  ## settings of x.test: only x.test[ , S]
                      ## are used but they must all be given
                      S,       ## indices of subset
                      seed=99L,
                      mult.impute=1L,
                      comb.draw=0L,
                      hotd.var=FALSE, ## hot-deck variance adjustment
                      alpha=0.05, ## hot-deck symmetric credible interval
                      probs=c(0.025, 0.975),
                      ## hot-deck asymmetric credible interval
                      mc.cores=1L,
                      nice=19L)
{
    if(mc.cores>1L) {
        if(.Platform$OS.type!='unix')
            stop('parallel::mcparallel/mccollect do not exist on windows')
        else {
            RNGkind("L'Ecuyer-CMRG")
            set.seed(seed)
            parallel::mc.reset.stream()
        }
    } else { set.seed(seed) }

    for(v in S)
        if(any(is.na(x.test[ , v])))
            stop(paste0('x.test column with missing values:', v))

    P = ncol(x.train)

    if(!all(S %in% 1:P))
        stop('some elements of S are not in x.train')

    if(P!=ncol(x.test))
        stop('the number of columns in x.train and x.test are not the same')

    if(P!=length(object$treedraws$cutpoints))
        stop(paste0('the number of columns in x.train and\n',
                    'the length of cutpoints are not the same'))

    Q=nrow(x.test)
    for(i in 1:(Q-1))
        for(j in (i+1):Q)
            if(all(x.test[i, S]==x.test[j, S]))
                warning(paste0('Row ', i, ' and ', j,
                               ' of x.test are equal with respect to S'))
    M=P-length(S)

    if(mc.cores>1L) call.=mc.hotdeck
    else call.=hotdeck

    shap=call.(x.train, x.test, S, object$treedraws, mult.impute=mult.impute,
              hotd.var=hotd.var, alpha=alpha, probs=probs,
              mc.cores=mc.cores, nice=nice)

    shap.=list()

    if(hotd.var) {
        shap.$yhat.test.=shap$yhat.test.
        shap.$yhat.test.var.=shap$yhat.test.var.
    } else {
        shap.$yhat.test=shap$yhat.test
    }

    ## weighted difference
    if(M>0) {
        varprob=object$varprob.mean
        varprob[S]=0
        for(k in 1:M) {
            C=comb(M, k, (1:P)[-S])
            R=nrow(C)
            if(comb.draw>0 && comb.draw<R) {
                R=comb.draw
                C=matrix(0, nrow=R, ncol=k)
                for(i in 1:R) {
                    varprob.=varprob
                    for(j in 1:k) {
                        h=(t(rmultinom(1, 1, varprob.))%*%(1:P))[1, 1]
                        if(h %in% S)
                            stop('comb.draw picked a column in S')
                        else if(h %in% C[i, 1:j])
                            stop('comb.draw picked a column twice')
                        C[i, j]=h
                        varprob.[h]=0
                    }
                }
            }
            for(i in 1:R) {
                shap.in=call.(x.train, x.test, c(C[i, ], S),
                             object$treedraws, mult.impute=mult.impute,
                             hotd.var=hotd.var, alpha=alpha, probs=probs,
                             mc.cores=mc.cores, nice=nice)
                shap.ex=call.(x.train, x.test, C[i, ],
                             object$treedraws, mult.impute=mult.impute,
                             hotd.var=hotd.var, alpha=alpha, probs=probs,
                             mc.cores=mc.cores, nice=nice)
                if(hotd.var) {
                    shap.$yhat.test.=shap.$yhat.test.+
                        (shap.in$yhat.test.-shap.ex$yhat.test.)/R
                        ##(shap.in$yhat.test.-shap.ex$yhat.test.)/choose(M, k)
                    ## shap.var=shap.in$yhat.test.var.+shap.ex$yhat.test.var.
                    ## for(i in 1:Q) {
                    ##     C=2*var(shap.in$yhat.test.[ , i],
                    ##             shap.ex$yhat.test.[ , i])
                    ##     if(shap.var[i]>C)
                    ##         shap.var[i]=shap.var[i]-C
                    ## }
                    ## shap.$yhat.test.var.=shap.$yhat.test.var.+
                    ##     shap.var/(choose(M, k)^2)
                } else {
                    shap.$yhat.test=shap.$yhat.test+
                        (shap.in$yhat.test-shap.ex$yhat.test)/R
                        ##(shap.in$yhat.test-shap.ex$yhat.test)/choose(M, k)
                }
            }
        }
    }

    if(hotd.var) {
        shap.$yhat.test.=shap.$yhat.test./P
        shap.$yhat.test.mean =apply(shap.$yhat.test., 2, mean)
        shap.$yhat.test.lower=apply(shap.$yhat.test., 2, quantile,
                                    probs=probs[1])
        shap.$yhat.test.upper=apply(shap.$yhat.test., 2, quantile,
                                    probs=probs[2])
        ##shap.$yhat.test.var.=shap.$yhat.test.var./(P^2)
    } else {
        shap.$yhat.test=shap.$yhat.test/P
        shap.$yhat.test.mean =apply(shap.$yhat.test, 2, mean)
        shap.$yhat.test.lower=apply(shap.$yhat.test, 2, quantile,
                                    probs=probs[1])
        shap.$yhat.test.upper=apply(shap.$yhat.test, 2, quantile,
                                    probs=probs[2])
    }

    return(shap.)
}

## if(mc.cores>1L)
##     D=mc.hotdeck(x.train, x.test, S, object$treedraws,
##                  mc.cores=mc.cores, nice=nice)
## else D=hotdeck(x.train, x.test, S, object$treedraws)

## ## weighted difference
## if(M>0) {
##     for(k in 1:M) {
##         C=comb(M, k, (1:P)[-S])
##         R=nrow(C)
##         for(i in 1:R) {
##             if(mc.cores>1L)
##                 D=D+(mc.hotdeck(x.train, x.test, c(C[i, ], S),
##                                 object$treedraws,
##                                 mc.cores=mc.cores, nice=nice)-
##                      mc.hotdeck(x.train, x.test, C[i, ],
##                                 object$treedraws,
##                                 mc.cores=mc.cores, nice=nice))/
##                     choose(M, k)
##             else
##                 D=D+(hotdeck(x.train, x.test, c(C[i, ], S),
##                              object$treedraws)-
##                      hotdeck(x.train, x.test, C[i, ],
##                              object$treedraws))/choose(M, k)
##         }
##     }
## }
## return(D/P)

