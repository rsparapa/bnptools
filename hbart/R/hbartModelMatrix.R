hbartModelMatrix=function(xdf) {
#note: for BART, if a factor has k levels, use all k dummies (not k-1 as in linear regression).
p=ncol(xdf)
Xnum = NULL #this will be all the numeric columns
Xfac = NULL #this will be all the factors expanded into dummies
for(i in 1:p) { # for each xi add to Xnum if numeric else add dummies to Xfac
   xnm = names(xdf)[i]
   if(is.factor(xdf[[i]])) {
      Xtemp = nnet::class.ind(xdf[[i]])
      colnames(Xtemp) = paste(xnm,1:ncol(Xtemp),sep='')
      Xfac = cbind(Xfac,Xtemp)
   } else {
      Xnum=cbind(Xnum,xdf[[i]])
      colnames(Xnum)[ncol(Xnum)]=xnm
   }
}
return(cbind(Xnum,Xfac))
}
