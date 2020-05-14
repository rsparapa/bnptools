/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2018 Robert McCulloch and Rodney Sparapani
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

RcppExport SEXP chotdeck(
   SEXP _itrain,
   SEXP _itest,
   SEXP _imask,
   SEXP _itrees,		//treedraws list from fbart
   SEXP _itc			//thread count
)
{
   Rprintf("*****In main of C++ for bart prediction\n");
   //random number generation
   arn gen;
   //--------------------------------------------------
   //get threadcount
   int tc = Rcpp::as<int>(_itc);
   cout << "tc (threadcount): " << tc << endl;
   //--------------------------------------------------
   //process trees
   Rcpp::List trees(_itrees);

   Rcpp::CharacterVector itrees(Rcpp::wrap(trees["trees"])); 
   std::string itv(itrees[0]);
   std::stringstream ttss(itv);

   size_t nd,m,p;
   ttss >> nd >> m >> p;
   cout << "number of bart draws: " << nd << endl;
   cout << "number of trees in bart sum: " << m << endl;
   cout << "number of x columns: " << p << endl;
   //--------------------------------------------------
   //process cutpoints (from trees)
   Rcpp::List  ixi(Rcpp::wrap(trees["cutpoints"]));
   size_t pp = ixi.size();
   if(p!=pp) cout << "WARNING: p from trees and p from x don't agree\n";
   xinfo xi;
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      Rcpp::NumericVector cutv(ixi[i]);
      xi[i].resize(cutv.size());
      std::copy(cutv.begin(),cutv.end(),xi[i].begin());
   }
   //--------------------------------------------------
   //process x
   Rcpp::IntegerVector mask(_imask);
   Rcpp::NumericMatrix xtrain(_itrain);
   Rcpp::NumericMatrix xtest(_itest);
   size_t n = xtrain.ncol();
   size_t np = xtest.ncol();
   cout << "from x,np,p: " << xtest.nrow() << ", " << xtest.ncol() << endl;
   //--------------------------------------------------
   //read in trees
   std::vector<vtree> tmat(nd);
   for(size_t i=0;i<nd;i++) tmat[i].resize(m);
   for(size_t i=0;i<nd;i++) {
      for(size_t j=0;j<m;j++) ttss >> tmat[i][j];
   }
   //--------------------------------------------------
   //get predictions

   Rcpp::NumericMatrix yhat(nd,np);
   std::fill(yhat.begin(), yhat.end(), 0.0);
   double *px = &xtest(0,0);

   for(size_t k=0; k<nd; ++k) {   // samples
     for(size_t i=0; i<np; ++i) { // settings
       size_t h;
       h=n*gen.uniform();
       for(size_t j=0; j<p; ++j)  // variables: rows of xpred
	 if(mask[j]==0) xtest(j, i)=xtrain(j, h);
     }
     getpred(k, k, p, m, np, xi, tmat, px, yhat);
   }

   //--------------------------------------------------
   Rcpp::List ret;
   ret["yhat.test"] = yhat;
   return ret;
}
