/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         Hugh A. Chipman and Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * cprnft2.h
 *
 * nftbart is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * nftbart is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author contact information
 * Matthew T. Pratola: mpratola@gmail.com
 * Robert E. McCulloch: robert.e.mculloch@gmail.com
 * Hugh A. Chipman: hughchipman@gmail.com
 * Rodney A. Sparapani: rsparapa@mcw.edu
 *
 */

/*
#include "rn.h"

#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

#include <fstream>
#include <vector>
#include <limits>

#include "ambrt.h"
#include "psbrt.h"
#include "brt.h"
#include "brtfuns.h"
#include "dinfo.h"
#include "mbrt.h"
#include "treefuns.h"
#include "tree.h"
*/

RcppExport SEXP cprnft2(
   SEXP _itrees,		//treedraws list from fbart
   SEXP _ixf,			
   SEXP _ixs,			
   SEXP _ixif,			
   SEXP _ixis,			
   SEXP _itc			//thread count
   //SEXP _idraws
)
{
   //Rprintf("*****In main of C++ for bart prediction\n");
   //--------------------------------------------------
   //get threadcount
   int /*draws = Rcpp::as<int>(_idraws),*/ 
     tc = Rcpp::as<int>(_itc);
   COUT << "tc (threadcount): " << tc << endl;
   //--------------------------------------------------
   //process trees
   Rcpp::List trees(_itrees);
   Rcpp::CharacterVector ftrees(Rcpp::wrap(trees["f.trees"]));
   std::string ftv(ftrees[0]);
   std::stringstream ftss(ftv);

   size_t nd,m,mh,pf,ps;
   ftss >> nd >> m >> pf;
   COUT << "number of bart draws: " << nd << endl;
   COUT << "number of trees in f: " << m << endl;
   COUT << "number of x columns: " << pf << endl;
   //--------------------------------------------------
   //process cutpoints (from trees)
   //Rcpp::List  ixi(Rcpp::wrap(trees["xicuts"]));
   Rcpp::List ixif(_ixif);
   size_t _pf = ixif.size();
   if(pf!=_pf) COUT << "WARNING: pf from trees and xi don't agree\n";
   xinfo xif;
   xif.resize(pf);
   for(size_t i=0;i<pf;i++) {
      Rcpp::NumericVector cutv(ixif[i]);
      xif[i].resize(cutv.size());
      std::copy(cutv.begin(),cutv.end(),xif[i].begin());
   }
   //--------------------------------------------------
   //process x
   Rcpp::NumericMatrix xfpred(_ixf), xspred(_ixs);
   size_t np = xfpred.ncol();
   //COUT << "from x,np,p: " << xpred.nrow() << ", " << xpred.ncol() << endl;
   //--------------------------------------------------
   //read in trees
   std::vector<vtree> fmat(nd);
   for(size_t i=0;i<nd;i++) {
     fmat[i].resize(m);
     for(size_t j=0;j<m;j++) ftss >> fmat[i][j];
   }
   //--------------------------------------------------
   //get predictions

   Rcpp::NumericMatrix fhat(nd,np);
   std::fill(fhat.begin(), fhat.end(), 0.);
   double *xfp = &xfpred(0,0);

   #ifndef _OPENMP
   tc=1;
   #endif
   if(tc==1) {
     COUT << "***using serial code\n"; 
     getpred(0, nd-1, pf, m,  np,  xif,  fmat, xfp,  fhat);
   }
   #ifdef _OPENMP
   else {
      COUT << "***using parallel code\n";
#pragma omp parallel num_threads(tc)
      local_getpred(nd,pf,m, np,xif,fmat,xfp,fhat);
   }
   #endif

   Rcpp::List ret;
   ret["f.test"] = fhat;

//if(draws) 
   {
   Rcpp::CharacterVector strees(Rcpp::wrap(trees["s.trees"])); 
   std::string stv(strees[0]);
   std::stringstream stss(stv);

   stss >> nd >> mh >> ps;
   COUT << "number of bart draws: " << nd << endl;
   COUT << "number of trees in s: " << mh << endl;
   COUT << "number of x columns: " << ps << endl;
   //--------------------------------------------------
   Rcpp::List ixis(_ixis);
   size_t _ps = ixis.size();
   if(ps!=_ps) COUT << "WARNING: ps from trees and xi don't agree\n";
   xinfo xis;
   xis.resize(ps);
   for(size_t i=0;i<ps;i++) {
      Rcpp::NumericVector cutv(ixis[i]);
      xis[i].resize(cutv.size());
      std::copy(cutv.begin(),cutv.end(),xis[i].begin());
   }
   //read in trees
   std::vector<vtree> smat(nd);
   for(size_t i=0;i<nd;i++) {
     smat[i].resize(mh);
     for(size_t j=0;j<mh;j++) stss >> smat[i][j];
   }
   //--------------------------------------------------
   //get predictions

   Rcpp::NumericMatrix shat(nd,np);
   std::fill(shat.begin(), shat.end(), 0.);
   double *xsp = &xfpred(0,0);

   if(tc==1) {
     COUT << "***using serial code\n"; 
     getpred(0, nd-1, ps, mh, np,  xis,  smat, xsp,  shat);
   }
   #ifdef _OPENMP
   else {
      COUT << "***using parallel code\n";
#pragma omp parallel num_threads(tc)
      local_getpred(nd,ps,mh,np,xis,smat,xsp,shat);
   }
   #endif

   ret["s.test"] = shat;
}
   return ret;
}

