
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

perm=function(M, K, S=1:M) {
    if(!all(sort(S)==unique(sort(S))))
        stop('S must consist of unique integers')
    if(K>M)
        stop('K cannot be greater than M')
    L=factorial(M)/factorial(M-K)
    C=matrix(nrow=L, ncol=K)
    if(K==0) {
        return(C)
    } else if(K==1) {
        C[ , 1]=sort(S[1:M])
        return(C)
    } else {
        a='0'
        for(N in 1:K) {
            b=paste0('i', N)
            if(N==1) {
                cmd=paste0('for(', b, ' in 1:M)\n')
                a=b
            }
            else {
                cmd=paste0(cmd, 'for(', b, ' in (1:M)[-c(', a, ')])\n')
                a=paste0(a, ', ', b)
            }
        }
        cmd=paste0(cmd, '{\nC[N, ]=c(')
        for(N in 1:(K-1)) cmd=paste0(cmd, 'S[', paste0('i', N), '],')
        cmd=paste0(cmd, 'S[', b, '])\n')
        cmd=paste0(cmd, 'N=N+1\n}')
        N=1
        eval(parse(text=cmd))
        for(i in 1:(L-1))
            for(j in (i+1):L)
                if(all(C[i, ]==C[j, ]))
                    warning(paste0('Rows ', i, ' and ', j, ' are the same'))
        return(C)
    }
}

perm.=function(M, K, S=1:M) {
    if(!all(sort(S)==unique(sort(S))))
        stop('S must consist of unique integers')
    if(K>M)
        stop('K cannot be greater than M')
    L=factorial(M)/factorial(M-K)
    C=matrix(nrow=L, ncol=K)
    if(K==0) {
        return(C)
    } else if(K==1) {
        C[ , 1]=S[1:M]
        return(C)
    } else {
        if(K>8)
            stop('K only supported up to 8 at this time')
        N=1
        for(a in 1:M) {
            for(b in (1:M)[-a]) {
                if(K==2) {
                    C[N, ]=c(S[a], S[b])
                    N=N+1
                } else {
                    for(c in (1:M)[-c(a, b)]) {
                        if(K==3) {
                            C[N, ]=c(S[a], S[b], S[c])
                            N=N+1
                        } else {
                            for(d in (1:M)[-c(a, b, c)]) {
                                if(K==4) {
                                    C[N, ]=c(S[a], S[b], S[c], S[d])
                                    N=N+1
                                } else {
                                    for(e in (1:M)[-c(a, b, c, d)]) {
                                        if(K==5) {
                                            C[N, ]=c(S[a], S[b], S[c], S[d], S[e])
                                            N=N+1
                                        } else {
                                            for(f in (1:M)[-c(a, b, c, d, e)]) {
                                                if(K==6) {
                                                    C[N, ]=c(S[a], S[b], S[c], S[d], S[e], S[f])
                                                    N=N+1
                                                } else {
                                                    for(g in (1:M)[-c(a, b, c, d, e, f)]) {
                                                        if(K==7) {
                                                            C[N, ]=c(S[a], S[b], S[c], S[d], S[e], S[f], S[g])
                                                            N=N+1
                                                        } else {
                                                            for(h in (1:M)[-c(a, b, c, d, e, f, g)]) {
                                                                if(K==8) {
                                                                    C[N, ]=c(S[a], S[b], S[c], S[d], S[e], S[f], S[g], S[h])
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
    for(i in 1:(L-1))
            for(j in (i+1):L)
                if(all(C[i, ]==C[j, ]))
                    warning(paste0('Rows ', i, ' and ', j, ' are the same'))
    return(C)
}

## for(K in 1:8)
##     for(M in K:8)
##         if(!all(perm(M, K)==perm.(M, K)))
##             print(c(M=M, K=K))
