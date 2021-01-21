predict.rbart = function(
object,
x.test=object$x.train,
tc=1,
fmean=mean(object$y.train),
q.lower=0.025,
q.upper=0.975,...)
{
#require(Rcpp)

#--------------------------------------------------
# params
nd=object$ndpost
m=object$ntree
mh=object$ntreeh
#numcut=object$numcut
xi=object$xicuts
p = ncol(object$x.train)
x = t(object$x.train)
xp = t(x.test)
if(is.null(object)) stop("No fitted model specified!\n")

#--------------------------------------------------
#call
res=.Call("cpsambrt_predict",
   x,
   xp,
   m,
   mh,
   nd,
#   numcut,
   xi,
   tc,
   object,
   PACKAGE="rbart"
)
res$mdraws=res$mdraws+fmean
res$mmean=apply(res$mdraws,2,mean)
res$smean=apply(res$sdraws,2,mean)
res$msd=apply(res$mdraws,2,sd)
res$ssd=apply(res$sdraws,2,sd)
res$m.5=apply(res$mdraws,2,quantile,0.5)
res$m.lower=apply(res$mdraws,2,quantile,q.lower)
res$m.upper=apply(res$mdraws,2,quantile,q.upper)
res$s.5=apply(res$sdraws,2,quantile,0.5)
res$s.lower=apply(res$sdraws,2,quantile,q.lower)
res$s.upper=apply(res$sdraws,2,quantile,q.upper)
res$q.lower=q.lower
res$q.upper=q.upper

return(res)
}



