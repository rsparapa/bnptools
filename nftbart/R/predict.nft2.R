## Copyright (C) 2021-2022 Rodney A. Sparapani

## This file is part of nftbart.
## predict.nft2.R

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

predict.nft2 = function(
                       ## data
                       object,
                       xftest=object$xftrain,
                       xstest=object$xstrain,
                       ## multi-threading
                       tc=getOption("mc.cores", 1), ##OpenMP thread count
                       ## current process fit vs. previous process fit
                       XPtr=TRUE,
                       ## predictions
                       K=0,
                       events=object$events,
                       FPD=FALSE,
                       hazard=FALSE,
                       density=FALSE,
                       probs=c(0.025, 0.975),
                       take.logs=TRUE,
                       na.rm=FALSE,
                       seed=NULL,
                       ## default settings for NFT:BART/HBART/DPM
                       fmu=object$fmu,
                       soffset=object$soffset,
                       drawDPM=object$drawDPM,
                       ##mask=FALSE,
                       ## etc.
                       ...)
{
    if(is.null(object)) stop("No fitted model specified!\n")
    np = nrow(xftest)
    if(np!=nrow(xstest))
        stop('The number of rows in xftest and xstest must be the same')
    ptm <- proc.time()
    xftest = t(xftest)
    xstest = t(xstest)
    ndpost=object$ndpost
    m=object$ntree[1]
    if(length(object$ntree)==2) mh=object$ntree[2]
    else mh=object$ntreeh
    xif=object$xifcuts
    xis=object$xiscuts
    n = nrow(object$xftrain)
    pf = ncol(object$xftrain)
    ps = ncol(object$xstrain)
    events.matrix=FALSE
    
    if(length(K)==0) {
        K=0
        take.logs=FALSE
    } else if(K>0) {
        if(length(events)==0) {
            events = unique(quantile(object$z.train.mean,
                                      probs=(1:K)/(K+1)))
            attr(events, 'names') = NULL
            take.logs=FALSE
            K = length(events)
        } else if(length(events)!=K) {
            stop("K and the length of events don't match")
        }
    } else if(K==0 && length(events)>0) {
        events.matrix=(class(events)[1]=='matrix')
        if(events.matrix) {
            if(FPD)
                stop("Friedman's partial dependence function: can't be used with a matrix of events")
            K=ncol(events)
        } else K = length(events)
    }
    if(K>0 && take.logs) events=log(events)
    
    if(length(drawDPM)==0) drawDPM=0

    draw.logt=(length(seed)>0)
    nd=length(object$s.train.mask)
    mask=(nd>0)
    if(!mask) nd=ndpost
    
    q.lower=min(probs)
    q.upper=max(probs)

    if(XPtr) {
        res=.Call("cpsambrt_predict2",
                  xftest, xstest,
                  m,
                  mh,
                  ndpost,
                  xif, xis,
                  tc,
                  object,
                  PACKAGE="nftbart"
                  )
        if(mask) {
            res$f.test.=res$f.test.[object$s.train.mask, ]
            res$s.test.=res$s.test.[object$s.train.mask, ]
        }
        res$f.test.=res$f.test.+fmu
        res$f.test.mean.=apply(res$f.test.,2,mean)
        res$f.test.lower.=
            apply(res$f.test.,2,quantile,probs=q.lower,na.rm=na.rm)
        res$f.test.upper.=
            apply(res$f.test.,2,quantile,probs=q.upper,na.rm=na.rm)
        res$s.test.mean.=apply(res$s.test.,2,mean)
        res$s.test.lower.=
            apply(res$s.test.,2,quantile,probs=q.lower,na.rm=na.rm)
        res$s.test.upper.=
            apply(res$s.test.,2,quantile,probs=q.upper,na.rm=na.rm)
    } else {
        res=list()
    }

    if(np>0) {
        res.=.Call("cprnft2",
                   object,
                   xftest, xstest,
                   xif, xis,
                   tc,
                   PACKAGE="nftbart"
                   )
        if(mask) {
            res.$f.test=res.$f.test[object$s.train.mask, ]
            res.$s.test=res.$s.test[object$s.train.mask, ]
        }
        res$f.test=res.$f.test+fmu
        res$fmu=fmu
        
        m=length(soffset)
        if(m==0) soffset=0
        else if(m>1) 
            soffset=sqrt(mean(soffset^2, na.rm=TRUE)) ## Inf set to NA
        
        res$s.test=exp(res.$s.test-soffset)
        res$s.test.mean =apply(res$s.test, 2, mean)
        res$s.test.lower=apply(res$s.test, 2,
                               quantile, probs=q.lower, na.rm=na.rm)
        res$s.test.upper=apply(res$s.test, 2,
                               quantile, probs=q.upper, na.rm=na.rm)
        res$soffset=soffset

        res$f.test.mean =apply(res$f.test, 2, mean)
        res$f.test.lower=apply(res$f.test, 2,
                               quantile,probs=q.lower,na.rm=na.rm)
        res$f.test.upper=apply(res$f.test, 2,
                               quantile,probs=q.upper,na.rm=na.rm)
        
        if(K>0) {
            ## if(length(events)==0) {
            ##     events <- unique(quantile(object$z.train.mean,
            ##                               probs=(1:K)/(K+1)))
            ##     attr(events, 'names') <- NULL
            ## } else if(take.logs) events=log(events)
            ## events.matrix=(class(events)[1]=='matrix')
            ## if(events.matrix) K=ncol(events)
            ## else K <- length(events)

            if(FPD) {
                H=np/n
                if(drawDPM>0) {
                    for(h in 1:H) {
                        if(h==1) {
                            mu. = object$dpmu
                            sd. = object$dpsd
                        } else {
                            mu. = cbind(mu., object$dpmu)
                            sd. = cbind(sd., object$dpsd)
                        }
                    }
                    mu. = mu.*res$s.test
                    sd. = sd.*res$s.test
                    mu. = mu.+res$f.test
                } else {
                    mu. = res$f.test
                    sd. = res$s.test
                }

                ## these FPD predict objects can become too large
                ## to be loaded in typical interactive sessions
                ## therefore, we drop the intermediate results
                if(K>1) {
                    res$f.test = NULL
                    if(length(res$f.test.)>0) res$f.test. = NULL
                    res$s.test = NULL
                    if(length(res$s.test.)>0) res$s.test. = NULL
                }
                
                surv.fpd=list()
                surv.test=list()
                if(hazard) {
                    haz.fpd =list()
                    haz.test =list()
                }
                if(density) {
                    pdf.fpd =list()
                    pdf.test =list()
                }
                for(i in 1:H) {
                    h=(i-1)*n+1:n
                    for(j in 1:K) {
                        if(j==1) {
                            surv.fpd[[i]]=list()
                            surv.test[[i]]=list()
                            if(hazard) {
                                haz.fpd[[i]] =list()
                                haz.test[[i]] =list()
                            }
                            if(density) {
                                pdf.fpd[[i]] =list()
                                pdf.test[[i]] =list()
                            }
                        }
                        z=events[j]
                        t=exp(z)
                        surv.fpd[[i]][[j]]=
                            apply(matrix(pnorm(z, mu.[ , h], sd.[ , h], lower.tail=FALSE),
                                         nrow=nd, ncol=n), 1, mean)
                        if(hazard)
                            haz.fpd[[i]][[j]]=
                            apply(matrix((dnorm(z, mu.[ , h], sd.[ , h])/(t*sd.[ , h]))
                                         /pnorm(z, mu.[ , h], sd.[ , h], lower.tail=FALSE),
                                         nrow=nd, ncol=n), 1, mean)
                        if(density)
                            pdf.fpd[[i]][[j]]=
                            apply(matrix(dnorm(z, mu.[ , h], sd.[ , h]),
                                         nrow=nd, ncol=n), 1, mean)
                        if(i==1 && j==1) {
                            res$surv.fpd=cbind(surv.fpd[[1]][[1]])
                            if(hazard) res$haz.fpd =cbind(haz.fpd[[1]][[1]])
                            if(density) res$pdf.fpd =cbind(pdf.fpd[[1]][[1]])
                        } else {
                            res$surv.fpd=cbind(res$surv.fpd, surv.fpd[[i]][[j]])
                            if(hazard) res$haz.fpd =cbind(res$haz.fpd, haz.fpd[[i]][[j]])
                            if(density) res$pdf.fpd =cbind(res$pdf.fpd, pdf.fpd[[i]][[j]])
                        }
                        if(K==1) {
                            surv.test[[i]][[j]]=matrix(pnorm(z, mu.[ , h], sd.[ , h], lower.tail=FALSE),
                                                       nrow=nd, ncol=n)
                            if(hazard)
                                haz.test[[i]][[j]] =matrix((dnorm(z, mu.[ , h], sd.[ , h])/(t*sd.[ , h]))/
                                                          pnorm(z, mu.[ , h], sd.[ , h], lower.tail=FALSE),
                                                       nrow=nd, ncol=n)
                            if(density)
                                pdf.test[[i]][[j]] =matrix(dnorm(z, mu.[ , h], sd.[ , h]),
                                                       nrow=nd, ncol=n)
                            if(i==1 && j==1) {
                                res$surv.test=cbind(surv.test[[1]][[1]])
                                if(hazard) res$haz.test =cbind(haz.test[[1]][[1]])
                                if(density) res$pdf.test =cbind(pdf.test[[1]][[1]])
                            } else {
                                res$surv.test=cbind(res$surv.test, surv.test[[i]][[j]])
                                if(hazard) res$haz.test =cbind(res$haz.test, haz.test[[i]][[j]])
                                if(density) res$pdf.test =cbind(res$pdf.test, pdf.test[[i]][[j]])
                            }
                        }
                    }
                }
                
                res$surv.fpd.mean=apply(cbind(res$surv.fpd), 2, mean)
                res$surv.fpd.lower=
                    apply(cbind(res$surv.fpd), 2, quantile, probs=q.lower)
                res$surv.fpd.upper=
                    apply(cbind(res$surv.fpd), 2, quantile, probs=q.upper)
                if(hazard) {
                    res$haz.fpd.mean =apply(cbind(res$haz.fpd), 2, mean)
                    res$haz.fpd.lower=
                        apply(cbind(res$haz.fpd), 2, quantile, probs=q.lower)
                    res$haz.fpd.upper=
                        apply(cbind(res$haz.fpd), 2, quantile, probs=q.upper)
                }
                if(density) 
                    res$pdf.fpd.mean =apply(cbind(res$pdf.fpd), 2, mean)
                if(K==1) {
                    res$surv.test.mean=apply(res$surv.test, 2, mean)
                    if(hazard)
                        res$haz.test.mean =apply(res$haz.test,  2, mean)
                    if(density)
                        res$pdf.test.mean =apply(res$pdf.test,  2, mean)
                }
            } else if(drawDPM>0) {
                res$surv.test=matrix(0, nrow=nd, ncol=np*K)
                if(hazard) res$haz.test=matrix(0, nrow=nd, ncol=np*K)
                ##res$aft.test=matrix(0, nrow=nd, ncol=np*K)
                H=ncol(object$dpwt.)
                ##H=max(c(object$dpn.))
                events.=events
                ## draw.time=(length(seed)>0)
                ## if(draw.time) {
                ##     set.seed(seed)
                ##     res$time.test=matrix(0, nrow=nd, ncol=np)
                ## }
                ## sd0=rep(0, nd)
                ## sd.=apply(log(res$s.test), 1, mean)
                ## for(h in 1:H) 
                ##    sd0=sd0+object$dpwt.[ , h]*sd.*log(object$dpsd.[ , h])
                ## sd0=exp(sd0)
                
                for(i in 1:np) {
                    if(events.matrix) events.=events[i, ]
                    mu.=res$f.test[ , i]
                    sd.=res$s.test[ , i]
                    ## if(draw.time)
                    ##     for(h in 1:H) {
                    ##         res$time.test[ , i]=res$time.test[ , i]+
                    ##             object$dpwt.[ , h]*
                    ##             rnorm(nd, 
                    ##                   mu.+sd.*object$dpmu.[ , h],
                    ##                   sd.*object$dpsd.[ , h])
                    ##     }
                    for(j in 1:K) {
                        k=(i-1)*K+j
                        z=events.[j]
                        t=exp(z)
                        for(h in 1:H) {
                            ##if(draws)
                            res$surv.test[ , k]=res$surv.test[ , k]+
                                object$dpwt.[ , h]*
                                pnorm(z, mu.+sd.*object$dpmu.[ , h],
                                      sd.*object$dpsd.[ , h], FALSE)
                            if(hazard)
                                res$haz.test[ , k]=res$haz.test[ , k]+
                                object$dpwt.[ , h]*
                                    dnorm(z, mu.+sd.*object$dpmu.[ , h],
                                      sd.*object$dpsd.[ , h])/
                                    (t*sd.*object$dpsd.[ , h]*
                                pnorm(z, mu.+sd.*object$dpmu.[ , h],
                                      sd.*object$dpsd.[ , h], FALSE))
                            ## res$aft.test[ , k]=res$aft.test[ , k]+
                            ##     object$dpwt.[ , h]*(
                            ##         dnorm(events.[j], mu.+sd0*object$dpmu.[ , h], sd0)/
                            ##         (exp(events.[j])*sd0))/
                            ##     pnorm(events.[j], mu.+sd0*object$dpmu.[ , h], sd0, FALSE)
                            ## else
                            ##     res$surv.test[ , k]=res$surv.test[ , k]+
                            ##         object$dpwt.[ , h]*
                            ##         pnorm(events.[j],
                            ##               mu.+object$dpmu.[ , h],
                            ##               object$sigma, FALSE)
                        }
                        ## if(hazard)
                        ##     for(h in 1:H)
                        ##         res$haz.test[ , k]=res$haz.test[ , k]+
                        ##         (object$dpwt.[ , h]*dnorm(z, mu.+sd.*object$dpmu.[ , h],
                        ##                                   sd.*object$dpsd.[ , h])/
                        ##             (t*sd.*object$dpsd.[ , h]*res$surv.test[ , k]))
                    }
                }
                res$surv.test.mean=apply(res$surv.test, 2, mean)
                if(hazard) res$haz.test.mean =apply(res$haz.test, 2, mean)
                ## res$aft.test.mean =apply(res$aft.test, 2, mean)
                ## res$aft.sd = sd0
                ## if(draw.time)
                ##     res$time.test.mean=apply(res$time.test, 2, mean)
            }

            if(take.logs) res$events=exp(events)
            else res$events=events
        }
        res$K=K
        
        if(draw.logt) {
            set.seed(seed)
            ##res$logt.test=matrix(0, nrow=nd, ncol=np)

            if(drawDPM>0) {
                res$logt.test=matrix(0, nrow=nd, ncol=np)
                H=ncol(object$dpwt.)
                for(i in 1:np) {
                    mu. = res$f.test[ , i]
                    sd. = res$s.test[ , i]
                    for(h in 1:H) {
                        res$logt.test[ , i]=res$logt.test[ , i]+
                            object$dpwt.[ , h]*
                            rnorm(nd, mu.+sd.*object$dpmu.[ , h],
                                  sd.*object$dpsd.[ , h])
                    }
                }
            } else {
                res$logt.test=matrix(nrow=nd, ncol=np)
                for(i in 1:np) {
                    mu. = res$f.test[ , i]
                    sd. = res$s.test[ , i]
                    res$logt.test[ , i]=rnorm(nd, mu., sd.)
                }
            }
            res$logt.test.mean=apply(res$logt.test, 2, mean)
        }
    } else {
        res$f.test=res$f.test.
        res$s.test=res$s.test.
    }

    res$probs=probs
    res$elapsed <- (proc.time()-ptm)['elapsed']
    attr(res$elapsed, 'names')=NULL

    return(res)
}



