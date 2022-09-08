## Copyright (C) 2021-2022 Rodney A. Sparapani

## This file is part of nftbart.
## nft.R

## nftbart is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.

## nftbart is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Author contact information
## Rodney A. Sparapani: rsparapa@mcw.edu

## for code maintainability, nft is just a wrapper of nft2

nft = function(## data
               x.train, times, delta=NULL, 
               x.test=matrix(nrow=0, ncol=0),
               rm.const=TRUE, rm.dupe=TRUE,
               impute.bin=NULL, impute.prob=NULL,
               ## multi-threading
               tc=getOption("mc.cores", 1), ##OpenMP thread count
               ##MCMC
               nskip=1000, ndpost=2000, 
               nadapt=1000, adaptevery=100,
               chv = NULL,
               method="spearman", use="pairwise.complete.obs",
               pbd=c(0.7, 0.7), pb=c(0.5, 0.5),
               stepwpert=c(0.1, 0.1), probchv=c(0.1, 0.1),
               minnumbot=c(5, 5),
               ## BART and HBART prior parameters
               ntree=c(50, 10), numcut=100, xicuts=NULL,
               power=c(2, 2), base=c(0.95, 0.95),
               ## f function
               k=5, sigmaf=NA,
               dist='weibull', 
               ## s function
               sigmav=NULL, total.lambda=NA, total.nu=10, mask=NULL,
               ## survival analysis 
               K=100, events=NULL, TSVS=FALSE,
               ## DPM LIO
               drawDPM=1L, 
               alpha=1, alpha.a=1, alpha.b=0.1, alpha.draw=1,
               neal.m=2, constrain=1, 
               m0=0, k0.a=1.5, k0.b=7.5, k0=1, k0.draw=1,
               a0=3, b0.a=2, b0.b=1, b0=1, b0.draw=1,
               ## misc
               na.rm=FALSE, probs=c(0.025, 0.975), printevery=100,
               transposed=FALSE ##, type='nft'
               )
{
res=nft2(## data
               xftrain=x.train, xstrain=x.train,
               times=times, delta=delta, 
               xftest=x.test,
               xstest=x.test,
               rm.const=rm.const, rm.dupe=rm.dupe,
               ## multi-threading
               tc=tc, ##OpenMP thread count
               ##MCMC
               nskip=nskip, ndpost=ndpost, 
               nadapt=nadapt, adaptevery=adaptevery, 
               chvf = chv, chvs = chv,
               method=method, use=use,
               pbd=pbd, pb=pb,
               stepwpert=stepwpert, probchv=probchv,
               minnumbot=minnumbot,
               ## BART and HBART prior parameters
               ntree=ntree, numcut=numcut,
               xifcuts=xicuts, xiscuts=xicuts,
               power=power, base=base,
               ## f function
               k=k, sigmaf=sigmaf, dist=dist,
               ## s function
               sigmav=sigmav, total.lambda=total.lambda,
               total.nu=total.nu, mask=mask,
               ## survival analysis 
               K=K, events=events, TSVS=TSVS,
               ## DPM LIO
               drawDPM=drawDPM, 
               alpha=alpha, alpha.a=alpha.a, alpha.b=alpha.b,
               alpha.draw=alpha.draw, neal.m=neal.m, constrain=constrain, 
               m0=m0, k0.a=k0.a, k0.b=k0.b, k0=k0, k0.draw=k0.draw,
               a0=a0, b0.a=b0.a, b0.b=b0.b, b0=b0, b0.draw=b0.draw,
               ## misc
               na.rm=na.rm, probs=probs, printevery=printevery,
               transposed=transposed
               )
attr(res, 'class') <- 'nft'
res$x.train=res$xftrain
res$xftrain=NULL
res$xstrain=NULL
res$xicuts=res$xifcuts
res$xifcuts=NULL
res$xiscuts=NULL
np=length(res$xftest)
if(np>0) {
    res$x.test=res$xftest
    res$xftest=NULL
    res$xstest=NULL
}
return(res)
}
