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
                         object,
                         x.test=object$x.train,
                         mc.cores=1,
                         tc=mc.cores,
                         fmu=object$fmu,
                         probs=c(0.025, 0.975),
                         XPtr=TRUE,
                         draws=object$draws,
                         soffset=object$soffset,
                         drawMuTau=object$drawMuTau,
                         K=0,
                         events=object$times,
                         take.logs=TRUE,
                         FPD=FALSE,
                         AFT=FALSE,
                         na.rm=FALSE,
                         ...)
{
    ptm <- proc.time()
    nd=object$ndpost
    m=object$ntree[1]
    if(!draws) mh=0
    else if(length(object$ntree)==2) mh=object$ntree[2]
    else mh=object$ntreeh
    xi=object$xicuts
    n = nrow(object$x.train)
    p = ncol(object$x.train)
    np = nrow(x.test)
    x = t(object$x.train)
    xp = t(x.test)
    if(is.null(object)) stop("No fitted model specified!\n")
    if(length(draws)==0) draws=0
    if(length(K)==0) {
        K=0
        take.logs=FALSE
    }
    if(length(drawMuTau)==0) drawMuTau=0
    
    q.lower=min(probs)
    q.upper=max(probs)

    if(XPtr) {
        res=.Call("cpsambrt_predict",
                  x,
                  xp,
                  m,
                  mh,
                  nd,
                  xi,
                  tc,
                  object,
                  PACKAGE="nftbart"
                  )
        res$mdraws=res$mdraws+fmu
        res$mmean=apply(res$mdraws,2,mean)
        res$msd=apply(res$mdraws,2,sd)
        ##res$m.q2=apply(res$mdraws,2,quantile,0.5)
        res$m.lower=apply(res$mdraws,2,quantile,probs=q.lower,na.rm=na.rm)
        res$m.upper=apply(res$mdraws,2,quantile,probs=q.upper,na.rm=na.rm)
        if(draws) {
            res$smean=apply(res$sdraws,2,mean)
            res$ssd=apply(res$sdraws,2,sd)
            ##res$s.q2=apply(res$sdraws,2,quantile,0.5)
            res$s.lower=apply(res$sdraws,2,quantile,probs=q.lower,na.rm=na.rm)
            res$s.upper=apply(res$sdraws,2,quantile,probs=q.upper,na.rm=na.rm)
        }
    } else {
        res=list()
    }

    if(np>0) {
        res.=.Call("cprnft",
                   object,
                   xp,
                   tc,
                   as.integer(draws),
                   PACKAGE="nftbart"
                   )
        res$f.test=res.$f.test+fmu
        res$f.test.mask=is.finite(apply(abs(res$f.test), 1, max))
        res$fmu=fmu
        
        if(draws) {
            ## if(XPtr) {
            ##     res.$s.test.mean=apply(res.$s.test,2,mean)
            ##     fit=lm(log(res$smean)~res.$s.test.mean)
            ##     if(length(soffset)==0) soffset=-coef(fit)[1]
            ##     res$coef=coef(fit)
            ## } else {
            ##     if(length(soffset)==0) soffset=mh
            ## }

            m=length(soffset)
            if(m==0) soffset=0
            else if(m>1) 
                soffset=sqrt(mean(soffset^2, na.rm=TRUE)) ## Inf set to NA
            
            res$s.test=exp(res.$s.test-soffset)
            res$s.test.mask=is.finite(apply(res$s.test, 1, max))
            res$f.test.mask=(res$f.test.mask&res$s.test.mask)
            res$s.test.mean =apply(res$s.test[res$f.test.mask, ], 2, mean)
            res$s.test.lower=apply(res$s.test[res$f.test.mask, ], 2,
                                   quantile, probs=q.lower, na.rm=na.rm)
            res$s.test.upper=apply(res$s.test[res$f.test.mask, ], 2,
                                   quantile, probs=q.upper, na.rm=na.rm)
            res$soffset=soffset
        } 

        res$f.test.mean =apply(res$f.test[res$f.test.mask, ], 2, mean)
        res$f.test.lower=apply(res$f.test[res$f.test.mask, ], 2,
                               quantile,probs=q.lower,na.rm=na.rm)
        res$f.test.upper=apply(res$f.test[res$f.test.mask, ], 2,
                               quantile,probs=q.upper,na.rm=na.rm)
        
        if(K>0) {
            if(length(events)==0) {
                events <- unique(quantile(object$z.train.mean,
                                          probs=(1:K)/(K+1)))
                                          ##probs=(0:(K-1))/(K-1)))
                attr(events, 'names') <- NULL
            } else if(take.logs) events=log(events)
            events.matrix=(class(events)[1]=='matrix')
            if(events.matrix) K=ncol(events)
            else K <- length(events)

            if(FPD) {
                H=np/n
                if(drawMuTau>0) {
                    for(h in 1:H) {
                        if(h==1) {
                            mu. = object$dpmu
                            sd. = object$dpsd
                        } else {
                            mu. = cbind(mu., object$dpmu)
                            sd. = cbind(sd., object$dpsd)
                        }
                    }
                    if(draws) {
                        mu. = mu.*res$s.test
                        sd. = sd.*res$s.test
                    } 
                    mu. = mu.+res$f.test
                } else {
                    mu. = res$f.test
                    if(draws) sd. = res$s.test
                    else sd. = matrix(object$sigma, nrow=nd, ncol=np)
                }

                ## these FPD predict objects can become too large
                ## to be loaded in typical interactive sessions
                ## therefore, we drop the intermediate results
                if(K>1) {
                    res$f.test = NULL
                    res$s.test = NULL
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
                
                mask=res$f.test.mask
                
                res$surv.fpd.mean=apply(cbind(res$surv.fpd[mask, ]), 2, mean)
                res$surv.fpd.lower=
                    apply(cbind(res$surv.fpd[mask, ]), 2, quantile, probs=q.lower)
                res$surv.fpd.upper=
                    apply(cbind(res$surv.fpd[mask, ]), 2, quantile, probs=q.upper)
                res$pdf.fpd.mean =apply(cbind(res$pdf.fpd[mask, ]), 2, mean)
                if(K==1) {
                    res$surv.test.mean=apply(res$surv.test[mask, ], 2, mean)
                    res$pdf.test.mean =apply(res$pdf.test[mask, ],  2, mean)
                }
            } else if(drawMuTau>0) {
                res$surv.test=matrix(0, nrow=nd, ncol=np*K)
                ##H=ncol(object$dpmu.)
                H=max(c(object$dpn.))
                events.=events
                for(i in 1:np) {
                    if(events.matrix) events.=events[i, ]
                    mu.=res$f.test[ , i]
                    if(draws) sd.=res$s.test[ , i]
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
            }

            if(take.logs) res$times=exp(events)
            else res$times=events
        }
        res$K=K
    }

    res$probs=probs
    res$elapsed <- (proc.time()-ptm)['elapsed']
    attr(res$elapsed, 'names')=NULL

    return(res)
}



