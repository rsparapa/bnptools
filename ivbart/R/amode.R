amode = function(n,Imode) {### alpha value which gives p(k|..) with mode at Imode
   egamm = 0.5772156649
   temp = digamma(Imode) - log(egamm+log(n))
   return(exp(temp))
}
