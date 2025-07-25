predict.hbart = function(
                         object,
                         x.test=object$x.train,
                         tc=getOption('mc.cores', 1L),
                         fmean=mean(object$y.train),
                         probs=c(0.025, 0.975),
                         XPtr=TRUE,
                         soffset=object$soffset,
                         ...)
{
    nd=object$ndpost
    m=object$ntree[1]
    mh=object$ntreeh
    xi=object$xicuts
    p = ncol(object$x.train)
    np = nrow(x.test)
    x = t(object$x.train)
    xp = t(x.test)
    if(is.null(object)) stop("No fitted model specified!\n")

    q.lower=probs[1]
    q.upper=probs[2]

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
                  PACKAGE="hbart"
                  )
        res$f.mean=fmean
        res$f.train=res$f.train+fmean
        res$f.train.mean=apply(res$f.train,2,mean)
        ##res$s.mean=sqrt(mean(apply(res$s.train^2, 2, mean)))
        ##res$s.train=res$s.train/res$s.mean
        res$s.train.mean=apply(res$s.train,2,mean)
        res$s.mean=sqrt(mean(res$s.train.mean^2))
        if(np>0) {
            res$f.test=res$f.test+fmean
            res$f.test.mean=apply(res$f.test,2,mean)
            res$f.test.sd=apply(res$f.test,2,sd)
            ##res$m.q2=apply(res$f.test,2,quantile,0.5)
            res$f.test.lower=apply(res$f.test,2,quantile,q.lower)
            res$f.test.upper=apply(res$f.test,2,quantile,q.upper)
            res$s.test.mean=apply(res$s.test,2,mean)
            res$s.test.sd=apply(res$s.test,2,sd)
            ##res$s.q2=apply(res$s.test,2,quantile,0.5)
            res$s.test.lower=apply(res$s.test,2,quantile,q.lower)
            res$s.test.upper=apply(res$s.test,2,quantile,q.upper)
        }
    } else {
        res=list()
    }
    
    if(np>0) {
        res.=.Call("cphbart",
                   object,
                   xp,
                   tc,
                   PACKAGE="hbart"
                   )
        if(XPtr) {
            res.$s.test.mean=apply(res.$s.test,2,mean)
            fit=lm(log(res$s.test.mean)~res.$s.test.mean)
            if(length(soffset)==0) soffset=-coef(fit)[1]
            ## fit=lm(log(res$s.mean)~0+I(res.$s.test.mean-mh))
            ## if(length(soffset)==0) soffset=coef(fit)[1]
        } else {
            if(length(soffset)==0) soffset=mh
        }
        
        ##res$s.test=exp((res.$s.test-mh)*soffset)
        res$s.test=exp(res.$s.test-soffset)
        res$s.test.mean=apply(res$s.test,2,mean)
        res$s.test.lower=apply(res$s.test,2,quantile,q.lower)
        res$s.test.upper=apply(res$s.test,2,quantile,q.upper)
        res$soffset=soffset
        if(XPtr) res$coef=coef(fit)
        res$f.test=res.$f.test+fmean
        res$f.test.mean=apply(res$f.test,2,mean)
        res$f.test.lower=apply(res$f.test,2,quantile,q.lower)
        res$f.test.upper=apply(res$f.test,2,quantile,q.upper)
    }
    res$probs=probs

    return(res)
}



