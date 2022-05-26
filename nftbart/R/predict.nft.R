## Copyright (C) 2021 Rodney A. Sparapani

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

predict.nft = function(
                       ## data
                       object,
                       x.test=object$x.train,
                       ## multi-threading
                       tc=1, ##OpenMP thread count
                       ## current process fit vs. previous process fit
                       XPtr=TRUE,
                       ## predictions
                       K=0,
                       events=object$events,
                       FPD=FALSE,
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
    ptm <- proc.time()
    ndpost=object$ndpost
    m=object$ntree[1]
    if(length(object$ntree)==2) mh=object$ntree[2]
    else mh=object$ntreeh
    xi=object$xicuts
    n = nrow(object$x.train)
    p = ncol(object$x.train)
    np = nrow(x.test)
    xp = t(x.test)
    if(is.null(object)) stop("No fitted model specified!\n")
    if(length(K)==0) {
        K=0
        take.logs=FALSE
    }
    if(length(drawDPM)==0) drawDPM=0

    draw.logt=(length(seed)>0)
    nd=length(object$s.train.mask)
    mask=(nd>0)
    if(!mask) nd=ndpost
    
    q.lower=min(probs)
    q.upper=max(probs)

    if(XPtr) {
        res=.Call("cpsambrt_predict",
                  xp,
                  m,
                  mh,
                  ndpost,
                  xi,
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
        res.=.Call("cprnft",
                   object,
                   xp,
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
            if(length(events)==0) {
                events <- unique(quantile(object$z.train.mean,
                                          probs=(1:K)/(K+1)))
                attr(events, 'names') <- NULL
            } else if(take.logs) events=log(events)
            events.matrix=(class(events)[1]=='matrix')
            if(events.matrix) K=ncol(events)
            else K <- length(events)

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
                pdf.fpd =list()
                surv.test=list()
                pdf.test =list()
                for(i in 1:H) {
                    h=(i-1)*n+1:n
                    for(j in 1:K) {
                        if(j==1) {
                            surv.fpd[[i]]=list()
                            pdf.fpd[[i]] =list()
                            surv.test[[i]]=list()
                            pdf.test[[i]] =list()
                        }
                        surv.fpd[[i]][[j]]=
                            apply(matrix(pnorm(events[j],
                                               mu.[ , h], sd.[ , h], lower.tail=FALSE),
                                         nrow=nd, ncol=n), 1, mean)
                        pdf.fpd[[i]][[j]]=
                            apply(matrix(dnorm(events[j], mu.[ , h], sd.[ , h]),
                                         nrow=nd, ncol=n), 1, mean)
                        if(i==1 && j==1) {
                            res$surv.fpd=cbind(surv.fpd[[1]][[1]])
                            res$pdf.fpd =cbind(pdf.fpd[[1]][[1]])
                        } else {
                            res$surv.fpd=cbind(res$surv.fpd, surv.fpd[[i]][[j]])
                            res$pdf.fpd =cbind(res$pdf.fpd, pdf.fpd[[i]][[j]])
                        }
                        if(K==1) {
                            surv.test[[i]][[j]]=matrix(pnorm(events[j],
                                                             mu.[ , h], sd.[ , h], lower.tail=FALSE),
                                                       nrow=nd, ncol=n)
                            pdf.test[[i]][[j]] =matrix(dnorm(events[j], mu.[ , h], sd.[ , h]),
                                                       nrow=nd, ncol=n)
                            if(i==1 && j==1) {
                                res$surv.test=cbind(surv.test[[1]][[1]])
                                res$pdf.test =cbind(pdf.test[[1]][[1]])
                            } else {
                                res$surv.test=cbind(res$surv.test, surv.test[[i]][[j]])
                                res$pdf.test =cbind(res$pdf.test, pdf.test[[i]][[j]])
                            }
                        }
                    }
                }
                
                res$surv.fpd.mean=apply(cbind(res$surv.fpd), 2, mean)
                res$surv.fpd.lower=
                    apply(cbind(res$surv.fpd), 2, quantile, probs=q.lower)
                res$surv.fpd.upper=
                    apply(cbind(res$surv.fpd), 2, quantile, probs=q.upper)
                res$pdf.fpd.mean =apply(cbind(res$pdf.fpd), 2, mean)
                if(K==1) {
                    res$surv.test.mean=apply(res$surv.test, 2, mean)
                    res$pdf.test.mean =apply(res$pdf.test,  2, mean)
                }
            } else if(drawDPM>0) {
                res$surv.test=matrix(0, nrow=nd, ncol=np*K)
                H=ncol(object$dpwt.)
                ##H=max(c(object$dpn.))
                events.=events
                ## draw.time=(length(seed)>0)
                ## if(draw.time) {
                ##     set.seed(seed)
                ##     res$time.test=matrix(0, nrow=nd, ncol=np)
                ## }
                    
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
                        for(h in 1:H) {
                            ##if(draws)
                            res$surv.test[ , k]=res$surv.test[ , k]+
                                object$dpwt.[ , h]*
                                pnorm(events.[j],
                                      mu.+sd.*object$dpmu.[ , h],
                                      sd.*object$dpsd.[ , h], FALSE)
                            ## else
                            ##     res$surv.test[ , k]=res$surv.test[ , k]+
                            ##         object$dpwt.[ , h]*
                            ##         pnorm(events.[j],
                            ##               mu.+object$dpmu.[ , h],
                            ##               object$sigma, FALSE)
                        }
                    }
                }
                res$surv.test.mean=apply(res$surv.test, 2, mean)
                ## if(draw.time)
                ##     res$time.test.mean=apply(res$time.test, 2, mean)
            }

            if(take.logs) res$events=exp(events)
            else res$events=events
        }
        res$K=K
        
        if(draw.logt) {
            set.seed(seed)
            res$logt.test=matrix(0, nrow=nd, ncol=np)

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



