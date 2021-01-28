//     cmbrt.cpp: Simple 1-tree model.
//     Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
//
//     This file is part of hbart.
//
//     hbart is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 2 of the License, or
//     (at your option) any later version.
//
//     hbart is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//     Author contact information
//     Matthew T. Pratola: mpratola@gmail.com
//     Robert E. McCulloch: robert.e.mculloch@gmail.com
//     Hugh A. Chipman: hughchipman@gmail.com

#include "hbart.h"

/*
#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

#include <fstream>
#include <vector>
#include <limits>

#ifndef NotInR
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#endif

#include "brt.h"
#include "brtfuns.h"
#include "dinfo.h"
#include "mbrt.h"
#include "treefuns.h"
#include "tree.h"
#include "rrn.h"

#ifndef NotInR
#define COUT Rcpp::Rcout 
#else
#define COUT std::cout  << "std::cout "
#endif
*/

RcppExport SEXP cmbrt(
   SEXP _ix,
   SEXP _iy,
   SEXP _ixp,
   SEXP _ind,
   SEXP _iburn,
   SEXP _itau,
   SEXP _ibase,
   SEXP _ipower,
   SEXP _itc,
   SEXP _isigmav,
   SEXP _ichv,
   SEXP _ipbd,
   SEXP _ipb,
   SEXP _istepwpert,
   SEXP _iprobchv,
   SEXP _iminnumbot,
   SEXP _iprintevery
)
{
   Rprintf("*****Into main of mbrt\n");
   //-----------------------------------------------------------
   //random number generation
   GetRNGstate();
   //Rcpp::RNGScope scope;
   //rcpprn gen;
   rrn gen;
   //--------------------------------------------------
   //process args
   //x
   Rcpp::NumericMatrix xm(_ix);
   size_t p = xm.nrow();
   size_t n = xm.ncol();
   double *x = &xm[0];

   //y
   Rcpp::NumericVector yv(_iy);
   double *y = &yv[0];

   //x out-of-sample
   Rcpp::NumericMatrix xpm(_ixp);
   size_t np = xpm.ncol();
   double *xp;
   if(np)  {
      xp = &xpm[0];
   } else {
      xp = nullptr;
   }

   //mu prior
   double tau = Rcpp::as<double>(_itau);

   //nd and burn
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);

   //tree prior
   double alpha = Rcpp::as<double>(_ibase);
   double mybeta = Rcpp::as<double>(_ipower);

   //thread count
   int tc = Rcpp::as<int>(_itc);

   //sigma vector
   Rcpp::NumericVector sigmav(_isigmav);
   double *sig = &sigmav[0];

   //change variable
   Rcpp::NumericMatrix chvm(_ichv);
   COUT << "row, cols chvm: " << chvm.nrow() << ", " << chvm.ncol() << endl;
   std::vector<std::vector<double> > chgv(chvm.nrow());
   for(int i=0;i<chvm.nrow();i++) {
      chgv[i].resize(chvm.ncol());
      for(int j=0;j<chvm.ncol();j++) chgv[i][j]=chvm(i,j);
   }

   //control
   size_t printevery = Rcpp::as<int>(_iprintevery);
   double pbd = Rcpp::as<double>(_ipbd);
   double pb = Rcpp::as<double>(_ipb);
   double stepwpert = Rcpp::as<double>(_istepwpert);
   double probchv = Rcpp::as<double>(_iprobchv);
   size_t minnumbot = Rcpp::as<int>(_iminnumbot);


   //--------------------------------------------------
   //print args
   Rprintf("**********************\n");
   Rprintf("n: %ld\n",n);
   Rprintf("p: %ld\n",p);
   Rprintf("first row: %lf, %lf\n",x[0],x[p-1]);
   Rprintf("second row: %lf, %lf\n",x[p],x[2*p-1]);
   Rprintf("last row: %lf, %lf\n",x[(n-1)*p],x[n*p-1]);
   Rprintf("first and last y: %lf, %lf\n",y[0],y[n-1]);
   if(np) {
      Rprintf("np: %d\n",np);
      Rprintf("first row xp: %lf, %lf\n",xp[0],xp[p-1]);
      Rprintf("second row xp: %lf, %lf\n",xp[p],xp[2*p-1]);
      Rprintf("last row xp : %lf, %lf\n",xp[(np-1)*p],xp[np*p-1]);
   } else {
      Rprintf("no test observations\n");
   }
   Rprintf("tau: %lf\n",tau);
   Rprintf("burn (nskip): %ld\n",burn);
   Rprintf("nd (ndpost): %ld\n",nd);
   Rprintf("tree prior base: %lf\n",alpha);
   Rprintf("tree prior power: %lf\n",mybeta);
   Rprintf("thread count: %ld\n",tc);
   Rprintf("first and last sigmav: %lf, %lf\n",sigmav[0],sigmav[n-1]);
   Rprintf("chgv first row: %lf, %lf\n",chgv[0][0],chgv[0][p-1]);
   Rprintf("chgv last row: %lf, %lf\n",chgv[p-1][0],chgv[p-1][p-1]);
   Rprintf("prob birth/death: %lf\n",pbd);
   Rprintf("prob birth: %lf\n",pb);
   Rprintf("step width pert move: %lf\n",stepwpert);
   Rprintf("prob of a change var move : %lf\n",probchv);
   Rprintf("min num obs in bottom node: %ld\n",minnumbot);
   Rprintf("*****printevery: %d\n",printevery);


   //--------------------------------------------------
   //make xinfo
   xinfo xi;
   size_t nc=1000;
   makexinfo(p,n,&x[0],xi,nc);

   //--------------------------------------------------
   //dinfo
   dinfo di;
   di.n=n;di.p=p,di.x = x; di.y = y;

   //--------------------------------------------------
   //setup mbrt object
   mbrt mbm;

   //cutpoints
   mbm.setxi(&xi);    //set the cutpoints for this model object
   //data objects
   mbm.setdata(&di);  //set the data
   //thread count
   mbm.settc(tc);      //set the number of threads when using OpenMP, etc.
   //tree prior
   mbm.settp(alpha, //the alpha parameter in the tree depth penalty prior
         mybeta     //the beta parameter in the tree depth penalty prior
         );
   //MCMC info
   mbm.setmi(
         pbd,  //probability of birth/death
         pb,  //probability of birth
         minnumbot,    //minimum number of observations in a bottom node
         true, //do perturb/change variable proposal?
         stepwpert,  //initialize stepwidth for perturb proposal.  If no adaptation it is always this.
         probchv,  //probability of doing a change of variable proposal.  perturb prob=1-this.
         &chgv  //initialize the change of variable correlation matrix.
         );
   mbm.setci(tau,sig);

   //--------------------------------------------------
   //return data structures using Rcpp
   size_t nkeeptrain=nd;
   size_t nkeeptest=nd;
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,n);
   double *fp = new double[np];
   dinfo dip;
   dip.x = xp; dip.y=fp;dip.p = p;dip.n=n;
   //means
   Rcpp::NumericVector trmean(n); //train
   for(size_t i=0;i<n;i++) trmean[i]=0.0;

   //--------------------------------------------------
   //run mcmc
   for(size_t i=0;i<burn;i++) mbm.draw(gen);
   for(size_t i=0;i<nd;i++) {
      if((i % printevery) ==0) COUT << "draw " << i << endl;
      mbm.draw(gen);
      for(size_t j=0;j<n;j++) {trmean(j)+=mbm.f(j)/nd; trdraw(i,j)=mbm.f(j);}
      if(np) {
         mbm.predict(&dip);
         for(size_t j=0;j<np;j++) tedraw(i,j) = fp[j];
      }
   }

   Rcpp::List ret;
   ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   ret["yhat.test"]=tedraw;

   //trees
   Rcpp::List treesL;
   //treesL["trees"]=Rcpp::CharacterVector(treess.str());
   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }
   treesL["cutpoints"] = xiret;
   ret["treedraws"] = treesL;
   ret["check"] = "mbrt-check";

   PutRNGstate();

   if(fp) delete [] fp;
   return ret;
}

