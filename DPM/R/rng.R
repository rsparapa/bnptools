
 ## DPM: Dirichlet Process Mixtures With Low Information Omnibus Priors
 ## Copyright (C) 2019 Prakash Laud and Rodney Sparapani
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2 of the License, or
 ## (at your option) any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, a copy is available at
 ## https://www.R-project.org/Licenses/GPL-2

rcat=function(prob) .Call("call_rcat", prob, PACKAGE="DPM")+1

.rgamma=function(n, shape, rate=1)
    .Call("call_rgamma", n, shape, rate, PACKAGE="DPM")

rlgamma=function(n, shape) .Call("call_rlgamma", n, shape, PACKAGE="DPM")

rtnorm=function(n, tau, mean=0, sd=1)
    .Call("call_rtnorm", n, tau, mean, sd, PACKAGE="DPM")
