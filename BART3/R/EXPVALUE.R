EXPVALUE = function(trees, x.test, S)
{
    H = nrow(x.test)
    M = dim(trees)[1]
    T = dim(trees)[2]
    G = function(n) ## node
    {
        if(trees[i, j, n, 1]==2) return(trees[i, j, n, 4]) ## a leaf
        else { ## a branch
            v=trees[i, j, n, 2]
            c=trees[i, j, n, 3]
            n=2*n
            m=n+1
            if(v %in% S) {
                if(x.test[h, v]<c) return(G(n))
                else return(G(m))
            } else {
                a=trees[i, j, n, 5]
                b=trees[i, j, m, 5]
                return((a*G(n)+b*G(m))/(a+b))
            }
        }
    }
    A = matrix(nrow=M, ncol=H)
    B = matrix(nrow=M, ncol=T)
    for(h in 1:H) { ## settings
        for(i in 1:M) ## samples
            for(j in 1:T) ## trees
                B[i, j]=G(1)
        A[ , h]=apply(B, 1, sum)
    }
    return(A)
}
