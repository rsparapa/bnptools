/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2018-2021 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

typedef std::vector<tree> vtree;

#ifdef _OPENMP
void local_getpred(size_t nd, size_t p, size_t m, size_t np, xinfo& xi, 
std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat);
#endif

void getpred(int beg, int end, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat);

RcppExport SEXP cphbart(
   SEXP _itrees,		//treedraws list from fbart
   SEXP _ix,			//x matrix to predict at
   SEXP _itc			//thread count
)
{
   //Rprintf("*****In main of C++ for bart prediction\n");
   //--------------------------------------------------
   //get threadcount
   int tc = Rcpp::as<int>(_itc);
   COUT << "tc (threadcount): " << tc << endl;
   //--------------------------------------------------
   //process trees
   Rcpp::List trees(_itrees);
   Rcpp::CharacterVector ftrees(Rcpp::wrap(trees["f.trees"])),
     strees(Rcpp::wrap(trees["s.trees"])); 
   std::string ftv(ftrees[0]), stv(strees[0]);
   std::stringstream ftss(ftv), stss(stv);

   size_t nd,m,mh,p;
   ftss >> nd >> m >> p;
   stss >> nd >> mh >> p;
   COUT << "number of bart draws: " << nd << endl;
   COUT << "number of trees in f: " << m << endl;
   COUT << "number of trees in s: " << mh << endl;
   COUT << "number of x columns: " << p << endl;
   //--------------------------------------------------
   //process cutpoints (from trees)
   Rcpp::List  ixi(Rcpp::wrap(trees["xicuts"]));
   size_t pp = ixi.size();
   if(p!=pp) COUT << "WARNING: p from trees and p from x don't agree\n";
   xinfo xi;
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      Rcpp::NumericVector cutv(ixi[i]);
      xi[i].resize(cutv.size());
      std::copy(cutv.begin(),cutv.end(),xi[i].begin());
   }
   //--------------------------------------------------
   //process x
   Rcpp::NumericMatrix xpred(_ix);
   size_t np = xpred.ncol();
   COUT << "from x,np,p: " << xpred.nrow() << ", " << xpred.ncol() << endl;
   //--------------------------------------------------
   //read in trees
   std::vector<vtree> fmat(nd), smat(nd);
   for(size_t i=0;i<nd;i++) {
     fmat[i].resize(m);
     for(size_t j=0;j<m;j++) ftss >> fmat[i][j];
     smat[i].resize(mh);
     for(size_t j=0;j<mh;j++) stss >> smat[i][j];
   }
   //--------------------------------------------------
   //get predictions

   Rcpp::NumericMatrix fhat(nd,np), shat(nd,np);
   std::fill(fhat.begin(), fhat.end(), 0.);
   std::fill(shat.begin(), shat.end(), 0.);
   double *px = &xpred(0,0);

   #ifndef _OPENMP
   tc=1;
   #endif
   if(tc==1) {
     COUT << "***using serial code\n"; 
     getpred(0, nd-1, p, m,  np,  xi,  fmat, px,  fhat);
     getpred(0, nd-1, p, mh, np,  xi,  smat, px,  shat);
   }
   #ifdef _OPENMP
   else {
      COUT << "***using parallel code\n";
#pragma omp parallel num_threads(tc)
      local_getpred(nd,p,m, np,xi,fmat,px,fhat);
      local_getpred(nd,p,mh,np,xi,smat,px,shat);
   }
   #endif

   Rcpp::List ret;
   ret["f.test"] = fhat;
   ret["s.test"] = shat;
   return ret;
}

void getpred(int beg, int end, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat)
{
   double *fptemp = new double[np];

   for(int i=beg;i<=end;i++) {
      for(size_t j=0;j<m;j++) {
         fit(tmat[i][j],xi,p,np,px,fptemp);
         for(size_t k=0;k<np;k++) yhat(i,k) += fptemp[k];
      }
   }

   delete [] fptemp;
}
#ifdef _OPENMP
void local_getpred(size_t nd, size_t p, size_t m, size_t np, xinfo& xi, 
std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat)
{
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int h = nd/thread_count; int beg = my_rank*h; int end = beg+h-1;
   
   getpred(beg,end,p,m,np,xi,tmat,px,yhat);
}
#endif
