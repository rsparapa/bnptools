
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2020 Rodney Sparapani
## stratrs

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

stratrs <- function(y, C, seed=NULL)
{
    N <- length(y)
    ## if(C*(N%/%C)!=N)
    ##     stop("length of y is not divisible by C")
    if(length(seed)>0) set.seed(seed)
    old.order = 1:N
    new.order = order(y)
    old.order = old.order[new.order]
    strat=0
    for(i in 1:((N%/%C)+1)) {
        j=(i-1)*C
        strat[j+(1:C)]=sample.int(C)
    }
##    strat=strat[1:N]
    return(strat[order(old.order)])
}

## set.seed(12)
## y = (c(rnorm(100), rnorm(100, 5)))
## strat=stratrs(y, 2)
## print(table(strat))
## print(summary(y[strat==1]))
## print(summary(y[strat==2]))

##     set.seed(12)
##     x <- rbinom(25000, 1, 0.1)
##     a <- stratrs(x, 5)
##     table(a, x)
##     z <- pmin(rpois(25000, 0.8), 5)
##     b <- stratrs(z, 5)
##     table(b, z)

## old/dumb version
    ## stratrs <- function(y, C=5, P=0)
    ## {
    ##     N <- length(y)
    ##     rs <- integer(N)

    ##     if(P>0) {
    ##         ## for continuous variables
    ##         Q=quantile(y, probs=(1:(P-1))/P)

    ##         for(i in 1:P) {
    ##             if(i<P) {
    ##                 M <- sum(y<=Q[i])
    ##                 rs[y<=Q[i]] <- sample(1:M, M)%%C
    ##             } else {
    ##                 M <- sum(y>Q[P-1])
    ##                 rs[y>Q[P-1] ] <- sample(1:M, M)%%C
    ##             }
    ##         }
    ##     } else {
    ##         ## for categorical variables
    ##         strat <- unique(sort(y))

    ##         for(i in strat) {
    ##             M <- sum(y==i)
    ##             rs[y==i] <- sample(1:M, M)%%C
    ##         }
    ##     }

    ##     return(rs+1)
    ## }
