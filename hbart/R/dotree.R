dotree=function(x,tmat,check,tc=2) {
   res = .Call("cdotree",
            t(x),
            tmat,
            check,
            tc,
            PACKAGE="hbart")
   return(res)
}
