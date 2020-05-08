
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2020 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

comb=function(M, K, S=1:M) {
    if(!all(sort(S)==unique(sort(S))))
        stop('S must consist of unique integers')
    if(K>M)
        stop('K cannot be greater than M')
    L=choose(M, K)
    C=matrix(nrow=L, ncol=K)
    if(K==0) {
        return(C)
    } else if(K==1) {
        C[ , 1]=sort(S[1:M])
        return(C)
    } else if(K==M) {
        C[1, ]=sort(S[1:M])
        return(C)
    } else if(K>1) {
        ## if(K>26)
        ##     stop('K only supported up to 26 at this time')
        a='0'
        for(N in 1:(K-1)) {
            ## b=letters[N]
            b=paste0('i', N)
            if(N==1)
                cmd=paste0('for(', b, ' in (', a, '+1):(M-K+', N, '))\n')
            else cmd=paste0(cmd, 'for(', b, ' in (', a, '+1):(M-K+', N, '))\n')
            a=b
        }
        ## b=letters[K]
        b=paste0('i', K)
        cmd=paste0(cmd, 'for(', b, ' in (', a, '+1):M) {\n')
        cmd=paste0(cmd, 'C[N, ]=sort(c(')
        ##for(N in 1:(K-1)) cmd=paste0(cmd, 'S[', letters[N], '],')
        for(N in 1:(K-1)) cmd=paste0(cmd, 'S[', paste0('i', N), '],')
        cmd=paste0(cmd, 'S[', b, ']))\n')
        cmd=paste0(cmd, 'N=N+1\n }')
        N=1
        eval(parse(text=cmd))
        for(i in 1:L) {
            if(!all(C[i, ]==sort(C[i, ])))
                warning(paste0('Row ', i, ' is not sorted'))
            if(i<L)
                for(j in (i+1):L)
                    if(all(C[i, ]==C[j, ]))
                        warning(paste0('Rows ', i, ' and ', j, ' are the same'))
        }
        return(C)
    }
}


comb.=function(M, K, S=1:M) {
    if(!all(sort(S)==unique(sort(S))))
        stop('S must consist of unique integers')
    if(K>M)
        stop('K cannot be greater than M')
    L=choose(M, K)
    C=matrix(nrow=L, ncol=K)
    if(K==0) {
        return(C)
    } else if(K==1) {
        C[ , 1]=sort(S[1:M])
        return(C)
    } else if(K==M) {
        C[1, ]=sort(S[1:M])
        return(C)
    } else if(K>1) {
        if(K>14)
            stop('K only supported up to 14 at this time')
        N=1
        for(a in 1:(M-K+1)) {
            for(b in (a+1):(M-K+2)) {
                if(K==2) {
                    C[N, ]=sort(c(S[a], S[b]))
                    N=N+1
                } else {
                    for(c in (b+1):(M-K+3)) {
                        if(K==3) {
                            C[N, ]=sort(c(S[a], S[b], S[c]))
                            N=N+1
                        } else {
                            for(d in (c+1):(M-K+4)) {
                                if(K==4) {
                                    C[N, ]=sort(c(S[a], S[b], S[c], S[d]))
                                    N=N+1
                                } else {
                                    for(e in (d+1):(M-K+5)) {
                                        if(K==5) {
                                            C[N, ]=sort(c(S[a], S[b], S[c], S[d], S[e]))
                                            N=N+1
                                        } else {
                                            for(f in (e+1):(M-K+6)) {
                                                if(K==6) {
                                                    C[N, ]=sort(c(S[a], S[b], S[c], S[d], S[e], S[f]))
                                                    N=N+1
                                                } else {
                                                    for(g in (f+1):(M-K+7)) {
                                                        if(K==7) {
                                                            C[N, ]=sort(c(S[a], S[b], S[c], S[d], S[e], S[f], S[g]))
                                                            N=N+1
                                                        } else {
                                                            for(h in (g+1):(M-K+8)) {
                                                                if(K==8) {
                                                                    C[N, ]=sort(c(S[a], S[b], S[c], S[d], S[e], S[f], S[g], S[h]))
                                                                    N=N+1
                                                                } else {
                                                                    for(i in (h+1):(M-K+9)) {
                                                                        if(K==9) {
                                                                            C[N, ]=sort(c(S[a], S[b], S[c], S[d], S[e],
                                                                                          S[f], S[g], S[h], S[i]))
                                                                            N=N+1
                                                                        } else {
                                                                            for(j in (i+1):(M-K+10)) {
                                                                                if(K==10) {
                                                                                    C[N, ]=sort(c(S[a], S[b], S[c], S[d],
                                                                                                  S[e], S[f], S[g], S[h],
                                                                                                  S[i], S[j]))
                                                                                    N=N+1
                                                                                } else {
                                                                                    for(k in (j+1):(M-K+11)) {
                                                                                        if(K==11) {
                                                                                            C[N, ]=sort(c(S[a], S[b], S[c],
                                                                                                          S[d], S[e], S[f],
                                                                                                          S[g], S[h], S[i],
                                                                                                          S[j], S[k]))
                                                                                            N=N+1
                                                                                        } else {
                                                                                            for(l in (k+1):(M-K+12)) {
                                                                                                if(K==12) {
                                                                                                    C[N, ]=sort(c(S[a], S[b], S[c], S[d], S[e], S[f], S[g], S[h], S[i], S[j], S[k], S[l]))
                                                                                                    N=N+1
                                                                                                } else {
                                                                                                    for(m in (l+1):(M-K+13)) {
                                                                                                        if(K==13) {
                                                                                                            C[N, ]=sort(c(S[a], S[b], S[c], S[d], S[e], S[f], S[g], S[h], S[i], S[j], S[k], S[l], S[m]))
                                                                                                            N=N+1
                                                                                                        } else {
                                                                                                            for(n in
                                                                                                            (m+1):(M-K+14)) {
                                                                                                                if(K==14) {
                                                                                                                    C[N, ]=
                                                                                                                        sort(c(S[a], S[b], S[c], S[d], S[e], S[f], S[g], S[h], S[i], S[j], S[k], S[l], S[m], S[n]))
                                                                                                                    N=N+1
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for(i in 1:L) {
        if(!all(C[i, ]==sort(C[i, ])))
            warning(paste0('Row ', i, ' is not sorted'))
        if(i<L)
            for(j in (i+1):L)
                if(all(C[i, ]==C[j, ]))
                    warning(paste0('Rows ', i, ' and ', j, ' are the same'))
    }
    return(C)
}

## for(K in 1:14)
##     for(M in K:14)
##         if(!all(comb(M, K)==comb.(M, K)))
##             print(c(M=M, K=K))
