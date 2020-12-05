getaprior = function(n,Imin,Imax,psi,gs) {### prior on alpha, grid of size gs
   amin = amode(n,Imin)
   amax = amode(n,Imax)
   av =  seq(from=amin,to=amax,length.out=gs)
   temp = 1 - (av-amin)/(amax-amin)
   temp = temp^psi
   temp = temp/sum(temp)
   return(list(a=av,p=temp,amin=amin,amax=amax))
}
