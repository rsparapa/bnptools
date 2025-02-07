
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2025 Robert McCulloch and Rodney Sparapani
## SHNN2.R

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

## Shapley additive explanation for 2-way interactions by Nearest Neighbors
SHNN2=function(object,  ## object returned from BART
              x.test,  ## settings of x.test
              S,       ## indices of two variables
              x.train,
              seed,
              mult.impute=5L
)
{
    UseMethod('SHNN2')
}
