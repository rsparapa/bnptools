
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

## Hot deck (HD) SHAP additive explanation
HDSHAP=function(object,  ## object returned from BART
            x.train, ## x.train to estimate coverage
            x.test,  ## settings of x.test: only x.test[ , S]
                     ## are used but they must all be given
            S,       ## indices of subset
            seed=99L,
            mult.impute=1L,
            hotd.var=FALSE, ## hot-deck variance adjustment
            alpha=0.05, ## hot-deck symmetric credible interval
            probs=c(0.025, 0.975),
                     ## hot-deck asymmetric credible interval
            mc.cores=1L,
            nice=19L)
{
    UseMethod('HDSHAP')
}
