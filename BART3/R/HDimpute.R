
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani
## HDimpute.R

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

HDimpute=function(x.train,
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

    ## hot deck missing imputation
    ## ignore columns for multinomial imputation in training

    check=(Q>0 && Q==N) ## are x.train and x.test the same?

    if(check)
        for(i in 1:N) {
            for(j in 1:P) {
                check=((is.na(x.train[i, j]) &&
                        is.na(x.test[i, j])) ||
                       (!is.na(x.train[i, j]) &&
                        !is.na(x.test[i, j]) &&
                        x.train[i, j]==x.test[i, j]))
                if(!check) break
            }
            if(!check) break
        }

    for(i in 1:N)
        for(j in 1:P)
            if(impute.flag && !(j %in% impute.mult)) {
                k = is.na(x.train[i, ])
                if(impute.flag) k[impute.mult]=FALSE
                while(is.na(x.train[i, j])) {
                    h=sample.int(N, 1)
                    x.train[i, which(k)]=x.train[h, which(k)]
                }
            }

    if(check && !impute.flag) x.test=x.train
    else if(Q>0) {
        if(check) x.test=x.train ## to hot-deck impute.mult columns only
        for(i in 1:Q)
            for(j in 1:P) {
                k = is.na(x.train[i, ])
                while(is.na(x.test[i, j])) {
                    h=sample.int(Q, 1)
                    x.test[i, which(k)]=x.test[h, which(k)]
                }
            }
    }

    return(list(x.train=x.train, x.test=x.test))
}
