##     mbrt.R: R script wrapper to call single tree homoscedastic BART model.
##     Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
##
##     This file is part of LISA.
##
##     LISA is free software: you can redistribute it and/or modify
##     it under the terms of the GNU Affero General Public License as published by
##     the Free Software Foundation, either version 3 of the License, or
##     (at your option) any later version.
##
##     LISA is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU Affero General Public License for more details.
##
##     You should have received a copy of the GNU Affero General Public License
##     along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##     Author contact information
##     Matthew T. Pratola: mpratola@gmail.com
##     Robert E. McCulloch: robert.e.mculloch@gmail.com
##     Hugh A. Chipman: hughchipman@gmail.com


mbrt = function(
x.train,
y.train,
x.test=matrix(0.0,0,0),
ndpost=1000, nskip=100,
power=1.0, base=.95,
tc=1,
sigmav=rep(1,length(y.train)),
fmean=mean(y.train),
chv = cor(x.train,method="spearman"),
pbd=.5,
pb=.5,
stepwpert=.01,
probchv=.1,
minnumbot=5,
printevery=100
)
{
#require(Rcpp)

#--------------------------------------------------
nd = ndpost
burn = nskip
#--------------------------------------------------
#data
n = length(y.train)
p = ncol(x.train)
np = nrow(x.test)
x = t(x.train)
xp = t(x.test)
y.train = y.train-fmean
#--------------------------------------------------
rgy = range(y.train)
tau =  (rgy[2]-rgy[1])/2.0;

#--------------------------------------------------
#call
res=.Call("cmbrt",
   x,
   y.train,
   xp,
   nd,
   burn,
   tau,
   base,
   power,
   tc,
   sigmav,
   chv,
   pbd,
   pb,
   stepwpert,
   probchv,
   minnumbot,
   printevery,PACKAGE="rbart"
)
res$yhat.train.mean = res$yhat.train.mean+fmean
res$yhat.train = res$yhat.train+fmean
res$yhat.test = res$yhat.test+fmean
return(res)
}
