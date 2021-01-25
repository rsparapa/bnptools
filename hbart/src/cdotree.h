
#include "hbart.h"

/*
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>

#include "treefuns.h"
#include "tree.h"
#include "rrn.h"

#define COUT Rcpp::Rcout
*/

RcppExport SEXP cdotree(
       SEXP _x,
       SEXP _tmat,
       SEXP _check,
       SEXP _tc
)
{
   Rprintf("In cdotree \n");
   size_t check = Rcpp::as<int>(_check);
   Rprintf("in cdotree, check: %ld\n",check);

   //--------------------------------------------------
   //random number generator
   GetRNGstate();
   rrn gen;

   //--------------------------------------------------
   // process arguments
   // this is supposed to be 2xnc where nc is the number of possible values for each of the 2 x's.
   Rcpp::NumericMatrix x(_x);
   size_t nrx=x.nrow();
   size_t ncx=x.ncol();
   Rprintf("x: nrow: %d, ncol %d\n",nrx,ncx);

   Rcpp::NumericMatrix tmat(_tmat);
   int nr=tmat.nrow();
   Rprintf("tmat: nrow: %d, ncol %d\n",nr,tmat.ncol());

   size_t tc = Rcpp::as<int>(_tc);
   COUT << "In dotree, tc is" << tc << "but currently not used" << "\n";

   //--------------------------------------------------
   // get xinfo, uniform cuts on (0,1) for each of the 2 x's.
   //void makeUnifXinfo(size_t p,size_t nc,xinfo& xi);
   xinfo xi;
   makeUnifXinfo(2,ncx,xi);

   //--------------------------------------------------
   // make tree
   tree t;
   for(int i=0;i<nr;i++) {
      size_t nid = (size_t)tmat(i,0);
      size_t v = (size_t)tmat(i,1);
      size_t cut = (size_t)tmat(i,2);
      double muL = tmat(i,3);
      double muR = tmat(i,4);
      Rprintf("node: %d, v: %d, cut: %d, muL: %lf, muR: %lf\n",nid,v,cut,muL,muR); 
      t.birth(nid,v,cut,muL,muR);
   }

   //print out tree
   t.pr();

   //--------------------------------------------------
   //return data structures using Rcpp
   //draws
   size_t nfit = ncx*ncx;
   Rcpp::NumericMatrix fmat(nfit,3);

   //--------------------------------------------------
   // get fits
   tree::tree_p bp; //pointer to bottom node
   double *xin = new double[2];
   size_t onx=0;
   for(size_t i=0;i!=ncx;i++) {
      for(size_t j=0;j!=ncx;j++) {
         xin[0] = x(0,i);
         xin[1] = x(1,j);
         bp = t.bn(xin,xi);
         //os << x[0] << " " << x[1] << " " << bp->gettheta() << " " << bp->nid() << std::endl;
         fmat(onx,0) = xin[0]; fmat(onx,1) = xin[1]; fmat(onx,2) = bp->gettheta();
         onx += 1;
      }
   }
   delete[] xin;

   //--------------------------------------------------
   // random number generator close
   PutRNGstate();

   //--------------------------------------------------
   //  return info
   Rcpp::List ret;
   ret["fmat"] = fmat;

   return ret;
}
