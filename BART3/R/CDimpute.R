
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020-2022 Robert McCulloch and Rodney Sparapani
## CDimpute.R

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

CDimpute=function(x.train,
                  x.test=matrix(0, 0, 0),
                  impute.mult=NULL
                  )
{   ## NOT TRANSPOSED
    P = ncol(x.train)
    N = nrow(x.train)
    Q = nrow(x.test)

    if(Q>0 && P!=ncol(x.test))
        stop('The number of columns in x.train and x.test must be identical')

    check = length(impute.mult)
    if(check==1)
        stop("The number of multinomial columns must be greater than 1\nConvert a binary into two columns")

    impute.flag=(check>1)

    ## cold deck missing imputation
    ## ignore columns for multinomial imputation in training

    same = FALSE
    miss.train=apply(is.na(x.train), 2, sum)
    names(miss.train)=dimnames(x.train)[[2]]
    if(impute.flag) miss.train[impute.mult]=0
    else impute.mult=0
    miss.train.=(sum(miss.train)>0)

    miss.test.=0
    miss.test=NULL
    if(Q>0) {
        miss.test=apply(is.na(x.test), 2, sum)
        names(miss.test)=dimnames(x.test)[[2]]
        if(impute.flag) miss.test[impute.mult]=0
        miss.test.=(sum(miss.test)>0)
    }

    if(miss.train.>0 || miss.test.>0) {
        same=(Q>0 && Q==N) ## are x.train and x.test the same?

        if(same)
            for(i in 1:N) {
                for(j in 1:P) {
                    same=((is.na(x.train[i, j]) &&
                            is.na(x.test[i, j])) ||
                           (!is.na(x.train[i, j]) &&
                            !is.na(x.test[i, j]) &&
                            x.train[i, j]==x.test[i, j]))
                    if(!same) break
                }
                if(!same) break
            }

        if(!impute.flag) {
            for(i in 1:N) {
                k = is.na(x.train[i, ])
                while(any(k)) {
                    h=sample.int(N, 1)
                    x.train[i, which(k)]=x.train[h, which(k)]
                    k = is.na(x.train[i, ])
                }
            }
        } else {
            for(i in 1:N) {
                k = is.na(x.train[i, ])
                k[impute.mult]=FALSE
                while(any(k)) {
                    h=sample.int(N, 1)
                    x.train[i, which(k)]=x.train[h, which(k)]
                    k = is.na(x.train[i, ])
                    k[impute.mult]=FALSE
                }
            }
        }
        ## for(i in 1:N)
        ##     for(j in 1:P)
        ##         if(!(j %in% impute.mult)) {
        ##             k = is.na(x.train[i, ])
        ##             if(impute.flag) k[impute.mult]=FALSE
        ##             while(is.na(x.train[i, j])) {
        ##                 h=sample.int(N, 1)
        ##                 x.train[i, which(k)]=x.train[h, which(k)]
        ##             }
        ##         }

        if(same && !impute.flag) x.test=x.train
        else if(Q>0) {
            if(same) x.test=x.train ## to cold-deck impute.mult columns only
            for(i in 1:Q) {
                k = is.na(x.test[i, ])
                while(any(k)) {
                    h=sample.int(Q, 1)
                    x.test[i, which(k)]=x.test[h, which(k)]
                    k = is.na(x.test[i, ])
                }
            }
        }
    }
            ## for(i in 1:Q)
            ##     for(j in 1:P) {
            ##         k = is.na(x.test[i, ])
            ##         while(is.na(x.test[i, j])) {
            ##             h=sample.int(Q, 1)
            ##             x.test[i, which(k)]=x.test[h, which(k)]
            ##         }
            ##     }

    return(list(x.train=x.train, x.test=x.test,
                miss.train=miss.train, miss.test=miss.test,
                impute.flag=impute.flag, same=same))
}
