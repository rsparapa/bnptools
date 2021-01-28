plotFunctionDraws=function(fd,complevel=mean(fd),probs=c(.025,.975),
      xlab="posterior mean of function",ylab="posterior intervals",
      intervalcol="green",linecol="red",
      pts=NA,ptscol="blue",ptspch=1, ptscex=1, ...) {
   nd = nrow(fd); nx = ncol(fd)
   mvec = apply(fd,2,mean)
   oo = order(mvec)
   qm = apply(fd,2,quantile,probs=probs)
   plot(c(mvec[oo[1]],mvec[oo[nx]]),range(qm),type="n",xlab=xlab, ylab=ylab,...)
   for(i in 1:nx) lines(rep(mvec[oo[i]],2),qm[,oo[i]],col=intervalcol)
   abline(h=complevel,col=linecol,lwd=3)
   if(length(pts) == nx) points(mvec[oo],pts[oo],col=ptscol,pch=ptspch,cex=ptscex)
}
