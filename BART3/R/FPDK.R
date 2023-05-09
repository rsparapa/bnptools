
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020-2023 Robert McCulloch and Rodney Sparapani

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

## kernel sampling Friedman's partial dependence (FPD) function
FPDK=function(object, ## object returned from BART
             x.test,  ## settings of x.test: only x.test[ , S]
                      ## are used but they must all be given
             S,       ## indices of subset
             x.train, 
             mult.impute=4L,
             kern.var=TRUE, ## kernel sampling variance adjustment
             alpha=0.05, ## kernel sampling symmetric credible interval
             probs=c(0.025, 0.975),
                      ## kernel sampling asymmetric credible interval
             mc.cores=getOption('mc.cores', 1L),
             seed=99L,
             nice=19L)
{
    UseMethod('FPDK')
}
