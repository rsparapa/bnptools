
    double G(size_t T, size_t node_max,
	     SEXP _trees, SEXP _xtest, SEXP _mask,
	     size_t h, size_t i, size_t j, size_t n)
    {
    Rcpp::NumericVector trees(_trees); // l=0, ..., 4
    // last entry: M*T*node_max*5-1
    Rcpp::NumericMatrix x_test(_xtest);
    Rcpp::IntegerVector mask(_mask);
      size_t k;
      //k=i*T*node_max*5+j*node_max*5+n*5;
      k=((i*T+j)*node_max+n)*5;
      //max k = (((M-1)*T+(T-1))*node_max+node_max-1)*5
      //      = (M-1)*T*node_max*5 + (T-1)*node_max*5 + (node_max-1)*5 + 5 - 5
      //      = M*T*node_max*5 - 5
        //if(trees[i, j, n, 1]==2) return(trees[i, j, n, 4]) ## a leaf
        if(trees[k]==2) return trees[k+3]; //## a leaf
        else { //## a branch
	  size_t v, m;
	  double c;
            v=trees[k+1]-1; // convert from R to C/C++ index
            //v=trees[i, j, n, 2]
            c=trees[k+2];
            //c=trees[i, j, n, 3]
            n=2*(n+1)-1;
            m=n+1;
            if(mask[v]==1) {
                if(x_test(h, v)<c) return G(T, node_max, _trees, _xtest, _mask, h, i, j, n);
                else return G(T, node_max, _trees, _xtest, _mask, h, i, j, m);
            } else {
	      double a, b;
                a=trees[k+4];
                b=trees[k+4];
                //a=trees[i, j, n, 5]
                //b=trees[i, j, m, 5]
                return (a*G(T, node_max, _trees, _xtest, _mask, h, i, j, n)+
			b*G(T, node_max, _trees, _xtest, _mask, h, i, j, m))/(a+b);
            }
        }
    }

RcppExport SEXP cEXPVALUE(SEXP _H, SEXP _M, SEXP _T, SEXP _node_max,
			  SEXP _trees, SEXP _xtest, SEXP _mask)
{
    size_t H = Rcpp::as<int>(_H);
    size_t M = Rcpp::as<int>(_M); // i=0, ..., M-1
    size_t T = Rcpp::as<int>(_T); // j=0, ..., T-1
    size_t node_max = Rcpp::as<int>(_node_max); // n=0, ..., node_max-1
    Rcpp::NumericMatrix A(M, H);

    for(size_t h=0; h<H; ++h) { //## settings
        for(size_t i=0; i<M; ++i) //## samples
            for(size_t j=0; j<T; ++j) //## trees
	      A(i, h)=A(i, h)+G(T, node_max, _trees, _xtest, _mask, h, i, j, 0);
                //B(i, j)=G(1);
        //A[ , h]=apply(B, 1, sum);
    }
    return A;
}

