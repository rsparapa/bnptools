
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2025 Robert McCulloch and Rodney Sparapani
## SHNN.R

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

## Shapley additive explanation with Nearest Neighbors
SHNN=function(object,  ## object returned from BART
              x.test,  ## settings of x.test
              S,       ## indices of subset
              x.train,
              type,
              probs,
              ##call,
              seed,
              mult.impute
              )
{
    UseMethod('SHNN')
}
