predict.hbart = function(
                         object,
                         x.test=object$x.train,
                         tc=1,
                         fmean=mean(object$y.train),
                         probs=c(0.025, 0.975),
                         XPtr=TRUE,
                         b0=NULL,
                         ...)
{
    nd=object$ndpost
    m=object$ntree[1]
    mh=object$ntreeh
    xi=object$xicuts
    p = ncol(object$x.train)
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
        res$mdraws=res$mdraws+fmean
        res$mmean=apply(res$mdraws,2,mean)
        res$msd=apply(res$mdraws,2,sd)
        ##res$m.q2=apply(res$mdraws,2,quantile,0.5)
        res$m.lower=apply(res$mdraws,2,quantile,q.lower)
        res$m.upper=apply(res$mdraws,2,quantile,q.upper)
        res$smean=apply(res$sdraws,2,mean)
        res$ssd=apply(res$sdraws,2,sd)
        ##res$s.q2=apply(res$sdraws,2,quantile,0.5)
        res$s.lower=apply(res$sdraws,2,quantile,q.lower)
        res$s.upper=apply(res$sdraws,2,quantile,q.upper)
    } else {
        res=list()
    }

    res.=.Call("cphbart",
               object,
               xp,
               tc,
               PACKAGE="hbart"
               )
    if(XPtr) {
        res.$s.test.mean=apply(res.$s.test,2,mean)
        fit=lm(log(res$smean)~res.$s.test.mean)
        if(length(b0)==0) b0=-coef(fit)[1]
    } else {
        if(length(b0)==0) b0=mh
    }
    
    res$s.test=exp(res.$s.test-b0)
    res$s.test.mean=apply(res$s.test,2,mean)
    res$s.test.lower=apply(res$s.test,2,quantile,q.lower)
    res$s.test.upper=apply(res$s.test,2,quantile,q.upper)
    res$b0=b0
    if(XPtr) res$coef=coef(fit)
    res$f.test=res.$f.test+fmean
    res$f.test.mean=apply(res$f.test,2,mean)
    res$f.test.lower=apply(res$f.test,2,quantile,q.lower)
    res$f.test.upper=apply(res$f.test,2,quantile,q.upper)
    res$probs=probs

    return(res)
}



