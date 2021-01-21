
hbartqqplot=function(y,rbmod,nunif=10000,linecolor="red",linewd=3,...) {

   ##helper function
   qinsamp = function(y,ysamp) { ###get quantile of y in sample ysamp
   n=length(ysamp)
   return(which.min(abs(y-sort(ysamp)))/n)
   }

   ##helper function
   qsamp = function(y,yd) { ###get quantile of yi in ith column of yd
   nd=nrow(yd)
   n=ncol(yd)
   qvec=rep(0,n)
   for(i in 1:n) {
      qvec[i]=qinsamp(y[i],yd[,i])
   }
   return(qvec)
   }

   ##dims
   nd = nrow(rbmod$mdraws)
   np = ncol(rbmod$mdraws)

   ##draw from predictive
   pdraw = rbmod$mdraws + rbmod$sdraws * matrix(rnorm(nd*np),nrow=nd)

   ## predictive qqplot
   yquantile = qsamp(y,pdraw)
   unifdr = runif(nunif)
   qqplot(yquantile,unifdr,...)
   abline(0,1,col=linecolor,lwd=linewd)

   ##return yquantiles
   return(yquantile)
}
