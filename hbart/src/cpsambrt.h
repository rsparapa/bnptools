//     cpsambrt.cpp: Implement BART model interface for R.
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

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>

#include "ambrt.h"
#include "psbrt.h"
#include "brt.h"
#include "brtfuns.h"
#include "dinfo.h"
#include "mbrt.h"
#include "treefuns.h"
#include "tree.h"
#include "rrn.h"

using std::cout;
using std::endl;

#ifndef NotInR
#define COUT Rcpp::Rcout 
#else
#define COUT std::cout  << "std::cout "
#endif
*/

RcppExport SEXP cpsambrt(
   SEXP _ix, 		//train x
   SEXP _iy, 		//train y
   SEXP _ixp,		//test x
   SEXP _im, 		//num trees mean
   SEXP _imh,		//num trees var
   SEXP _ind,		//number of draws kept  (adapt, then burn, then draw)
   SEXP _iburn, 	//number of burn draws
   SEXP _inadapt, 	//number of adapt draws 
   SEXP _iadaptevery,	//adapt every adaptevery draws
   SEXP _itau,		//mu_i ~N(0,\tau^2)
   SEXP _ioveralllambda,//sigma_i ~ prod (nu'*lambda')/chisquared^_nu', nu',lambda' from overallnu, overalllambda
   SEXP _ioverallnu,
   SEXP _ibase,		//tree prior for mean
   SEXP _ipower,	//tree prior for mean
   SEXP _ibaseh,	//tree prior for var
   SEXP _ipowerh,	//tree prior for var
   SEXP _itc,		//thread count
   SEXP _isigmav,	//intial values for sigma_i 
   SEXP _ichv,		//change variable (correlation) matrix
   SEXP _ipbd,		//prob of birth/death step, mean trees
   SEXP _ipb,		//pb prob of birth given birth/death
   SEXP _ipbdh,		//prob of birth/death ,var trees
   SEXP _ipbh,		//prob of birth given birth/death, var trees
   SEXP _istepwpert,	//(in (0,1)) shrinkage factor for range of acceptable cutpoints for MH, mean, (adapts)
   SEXP _istepwperth,	//(in (0,1)) shrinkage factor for range of acceptable cutpoints for MH, var, (adapts)
   SEXP _iprobchv,	//prob of trying change variable move, mean, (prob of perturb is 1 minus this)
   SEXP _iprobchvh,	//prob of trying change variable move, var (prob of perturb is 1 minus this)
   SEXP _iminnumbot,	//minimum number of observations in a bottom node, mean
   SEXP _iminnumboth,   //minimum number of observations in a bottom node, var
   SEXP _iprintevery,   //how often you print progress
   SEXP _ixicuts, //variable cutpoints for each predictor
   SEXP _isummarystats  //boolean, do you want summary stats

)
{
   Rprintf("*****Into main of cpsambrt\n");


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
   double *xp = nullptr;
   if(np)  xp = &xpm[0];

   //number of trees
   size_t m = Rcpp::as<int>(_im);
   size_t mh = Rcpp::as<int>(_imh);

   //mu prior (tau, ambrt) and sigma prior (lambda,nu, psbrt)
   double tau = Rcpp::as<double>(_itau);
   double overalllambda = Rcpp::as<double>(_ioveralllambda);
   double overallnu = Rcpp::as<double>(_ioverallnu);

   //nd and burn
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   size_t nadapt = Rcpp::as<int>(_inadapt);
   size_t adaptevery = Rcpp::as<int>(_iadaptevery);

   //tree prior
   double alpha = Rcpp::as<double>(_ibase);
   double mybeta = Rcpp::as<double>(_ipower);
   double alphah = Rcpp::as<double>(_ibaseh);
   double mybetah = Rcpp::as<double>(_ipowerh);

   //thread count
   int tc = Rcpp::as<int>(_itc);

   //sigma vector
   Rcpp::NumericVector sigmav(_isigmav);
   double *sig = &sigmav[0];
   dinfo disig;
   disig.n=n;disig.p=p;disig.x = x; disig.y = sig; disig.tc=tc;

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
   double pbdh = Rcpp::as<double>(_ipbdh);
   double pbh = Rcpp::as<double>(_ipbh);
   double stepwpert = Rcpp::as<double>(_istepwpert);
   double stepwperth = Rcpp::as<double>(_istepwperth);
   double probchv = Rcpp::as<double>(_iprobchv);
   double probchvh = Rcpp::as<double>(_iprobchvh);
   size_t minnumbot = Rcpp::as<int>(_iminnumbot);
   size_t minnumboth = Rcpp::as<int>(_iminnumboth);
   bool dopert=true;
   bool doperth=true;
   if(probchv<0) dopert=false;
   if(probchvh<0) doperth=false;

   //summary statistics 
    unsigned int tmaxd=0;
    unsigned int tmind=0;
    double tavgd=0.0;
   Rcpp::IntegerMatrix fvc(nd, p), svc(nd, p);
   Rcpp::IntegerVector varcount(p);
   bool summarystats = Rcpp::as<bool>(_isummarystats);

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
   Rprintf("number of trees mean: %ld\n",m);
   Rprintf("number of trees stan dev: %ld\n",mh);
   Rprintf("tau: %lf\n",tau);
   Rprintf("overalllambda: %lf\n",overalllambda);
   Rprintf("overallnu: %lf\n",overallnu);
   Rprintf("burn (nskip): %ld\n",burn);
   Rprintf("nd (ndpost): %ld\n",nd);
   Rprintf("nadapt: %ld\n",nadapt);
   Rprintf("adaptevery: %ld\n",adaptevery);
   Rprintf("mean tree prior base: %lf\n",alpha);
   Rprintf("mean tree prior power: %lf\n",mybeta);
   Rprintf("variance tree prior base: %lf\n",alphah);
   Rprintf("variance tree prior power: %lf\n",mybetah);
   Rprintf("thread count: %ld\n",tc);
   Rprintf("first and last sigmav: %lf, %lf\n",sigmav[0],sigmav[n-1]);
   Rprintf("chgv first row: %lf, %lf\n",chgv[0][0],chgv[0][p-1]);
   Rprintf("chgv last row: %lf, %lf\n",chgv[p-1][0],chgv[p-1][p-1]);
   Rprintf("mean trees prob birth/death: %lf\n",pbd);
   Rprintf("mean trees prob birth: %lf\n",pb);
   Rprintf("variance trees prob birth/death: %lf\n",pbdh);
   Rprintf("variance trees prob birth: %lf\n",pbh);
   Rprintf("mean trees initial step width pert move: %lf\n",stepwpert);
   Rprintf("variance trees initial step width pert move: %lf\n",stepwperth);
   Rprintf("mean trees prob of a change var move : %lf\n",probchv);
   Rprintf("variance trees prob of a change var move : %lf\n",probchvh);
   Rprintf("mean trees min num obs in bottom node: %ld\n",minnumbot);
   Rprintf("variance trees min num obs in bottom node: %ld\n",minnumboth);
   Rprintf("*****printevery: %d\n",printevery);

   //--------------------------------------------------
   //make xinfo
   xinfo xi;
   Rcpp::List ixi(_ixicuts);
   if(ixi.length()!= (int)p) {
      Rprintf("Cutpoint definition does not match number of predictor variables!\n");
      return 0;
   }
   xi.resize(p);

   for(size_t i=0;i<p;i++)
      xi[i]=Rcpp::as< std::vector<double> >(ixi[i]);
   COUT << "&&& made xinfo\n";

   //summarize input variables:
   for(size_t i=0;i<p;i++)
   {
      COUT << "Variable " << i << " has numcuts=" << xi[i].size() << " : ";
      COUT << xi[i][0] << " ... " << xi[i][xi[i].size()-1] << endl;
   }


   //--------------------------------------------------
   //dinfo
   dinfo di;
   di.n=n;di.p=p;di.x = x; di.y = y; di.tc=tc;
   //--------------------------------------------------
   // set up ambrt object
   ambrt ambm(m);

   //cutpoints
   ambm.setxi(&xi);    //set the cutpoints for this model object
   //data objects
   ambm.setdata(&di);  //set the data
   //thread count
   ambm.settc(tc);      //set the number of threads when using OpenMP, etc.
   //tree prior
   ambm.settp(alpha, //the alpha parameter in the tree depth penalty prior
         mybeta     //the beta parameter in the tree depth penalty prior
         );
   //MCMC info
   ambm.setmi(
         pbd,  //probability of birth/death
         pb,  //probability of birth
         minnumbot,    //minimum number of observations in a bottom node
         dopert, //do perturb/change variable proposal?
         stepwpert,  //initialize stepwidth for perturb proposal.  If no adaptation it is always this.
         probchv,  //probability of doing a change of variable proposal.  perturb prob=1-this.
         &chgv  //initialize the change of variable correlation matrix.
         );
   ambm.setci(tau,sig);


   //--------------------------------------------------
   //setup psbrt object
   psbrt psbm(mh,overalllambda);

   //make di for psbrt object
   dinfo dips;
   double *r = new double[n];
   for(size_t i=0;i<n;i++) r[i]=sigmav[i];
   dips.x = x; dips.y=r; dips.p=p; dips.n=n; dips.tc=tc;

   double opm=1.0/((double)mh);
   double nu=2.0*pow(overallnu,opm)/(pow(overallnu,opm)-pow(overallnu-2.0,opm));
   double lambda=pow(overalllambda,opm);

   //cutpoints
   psbm.setxi(&xi);    //set the cutpoints for this model object
   //data objects
   psbm.setdata(&dips);  //set the data
   //thread count
   psbm.settc(tc); 
   //tree prior
   psbm.settp(alphah, //the alpha parameter in the tree depth penalty prior
         mybetah     //the beta parameter in the tree depth penalty prior
         );
   psbm.setmi(
         pbdh,  //probability of birth/death
         pbh,  //probability of birth
         minnumboth,    //minimum number of observations in a bottom node
         doperth, //do perturb/change variable proposal?
         stepwperth,  //initialize stepwidth for perturb proposal.  If no adaptation it is always this.
         probchvh,  //probability of doing a change of variable proposal.  perturb prob=1-this.
         &chgv  //initialize the change of variable correlation matrix.
         );
   psbm.setci(nu,lambda);

   //--------------------------------------------------
   //return data structures using Rcpp
   double *fp = new double[np];
   dinfo dip;
   dip.x = xp; dip.y=fp; dip.p = p; dip.n=np; dip.tc=tc;

   //--------------------------------------------------
   //run mcmc
   std::vector<int> onn(nd*m,1);
   std::vector<std::vector<int> > oid(nd*m, std::vector<int>(1));
   std::vector<std::vector<int> > ovar(nd*m, std::vector<int>(1));
   std::vector<std::vector<int> > oc(nd*m, std::vector<int>(1));
   std::vector<std::vector<double> > otheta(nd*m, std::vector<double>(1));
   std::vector<int> snn(nd*mh,1);
   std::vector<std::vector<int> > sid(nd*mh, std::vector<int>(1));
   std::vector<std::vector<int> > svar(nd*mh, std::vector<int>(1));
   std::vector<std::vector<int> > sc(nd*mh, std::vector<int>(1));
   std::vector<std::vector<double> > stheta(nd*mh, std::vector<double>(1));
   brtMethodWrapper fambm(&brt::f,ambm);
   brtMethodWrapper fpsbm(&brt::f,psbm);

   Rprintf("Starting MCMC...\n");
   for(size_t i=0;i<nadapt;i++) { 
      ambm.draw(gen); 
      dips = di;
      dips -= fambm;
      if((i+1)%adaptevery==0) ambm.adapt();
      psbm.draw(gen);
      disig = fpsbm;
      if((i+1)%adaptevery==0) psbm.adapt();      
   }
   COUT << endl;
   for(size_t i=0;i<burn;i++) {
      if((i % printevery) ==0) COUT << "draw  burn " << i << endl;
      ambm.draw(gen);
      dips = di;
      dips -= fambm;      
      psbm.draw(gen);
      disig = fpsbm;
   }
   if(summarystats) {
      ambm.setstats(true);
      psbm.setstats(true);
   }
   for(size_t i=0;i<nd;i++) {
      if((i % printevery) ==0) COUT << "draw  keep " << i << endl;
      ambm.draw(gen);
      dips = di;
      dips -= fambm;
      psbm.draw(gen);
      disig = fpsbm;

      //save tree to vec format
      ambm.savetree(i,m,onn,oid,ovar,oc,otheta);
      psbm.savetree(i,mh,snn,sid,svar,sc,stheta);

      ambm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
      for(size_t k=0;k<p;k++) fvc(i, k)=varcount[k];
      psbm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
      for(size_t k=0;k<p;k++) svc(i, k)=varcount[k];
   }

  std::stringstream mtrees, strees;  
  mtrees=ambm.gettrees(nd,m,onn,oid,ovar,oc,otheta,0.);
  strees=psbm.gettrees(nd,mh,snn,sid,svar,sc,stheta,-1.);

   //Flatten posterior trees to a few (very long) vectors so we can just pass pointers
   //to these vectors back to R (which is much much faster than copying all the data back).
   std::vector<int>* e_ots=new std::vector<int>(nd*m);
   std::vector<int>* e_oid=new std::vector<int>;
   std::vector<int>* e_ovar=new std::vector<int>;
   std::vector<int>* e_oc=new std::vector<int>;
   std::vector<double>* e_otheta=new std::vector<double>;
   std::vector<int>* e_sts=new std::vector<int>(nd*mh);
   std::vector<int>* e_sid=new std::vector<int>;
   std::vector<int>* e_svar=new std::vector<int>;
   std::vector<int>* e_sc=new std::vector<int>;
   std::vector<double>* e_stheta=new std::vector<double>;
   for(size_t i=0;i<nd;i++)
      for(size_t j=0;j<m;j++) {
         e_ots->at(i*m+j)=static_cast<int>(oid[i*m+j].size());
         e_oid->insert(e_oid->end(),oid[i*m+j].begin(),oid[i*m+j].end());
         e_ovar->insert(e_ovar->end(),ovar[i*m+j].begin(),ovar[i*m+j].end());
         e_oc->insert(e_oc->end(),oc[i*m+j].begin(),oc[i*m+j].end());
         e_otheta->insert(e_otheta->end(),otheta[i*m+j].begin(),otheta[i*m+j].end());
      }
   for(size_t i=0;i<nd;i++)
      for(size_t j=0;j<mh;j++) {
         e_sts->at(i*mh+j)=static_cast<int>(sid[i*mh+j].size());
         e_sid->insert(e_sid->end(),sid[i*mh+j].begin(),sid[i*mh+j].end());
         e_svar->insert(e_svar->end(),svar[i*mh+j].begin(),svar[i*mh+j].end());
         e_sc->insert(e_sc->end(),sc[i*mh+j].begin(),sc[i*mh+j].end());
         e_stheta->insert(e_stheta->end(),stheta[i*mh+j].begin(),stheta[i*mh+j].end());
      }

   Rcpp::XPtr< std::vector<int> > extern_ots(e_ots,true);
   Rcpp::XPtr< std::vector<int> > extern_oid(e_oid,true);
   Rcpp::XPtr< std::vector<int> > extern_ovar(e_ovar,true);
   Rcpp::XPtr< std::vector<int> > extern_oc(e_oc,true);
   Rcpp::XPtr< std::vector<double> > extern_otheta(e_otheta,true);
   Rcpp::XPtr< std::vector<int> > extern_sts(e_sts,true);
   Rcpp::XPtr< std::vector<int> > extern_sid(e_sid,true);
   Rcpp::XPtr< std::vector<int> > extern_svar(e_svar,true);
   Rcpp::XPtr< std::vector<int> > extern_sc(e_sc,true);
   Rcpp::XPtr< std::vector<double> > extern_stheta(e_stheta,true);

   Rcpp::List ret = Rcpp::List::create(
				       Rcpp::Named("ots")=extern_ots,
                                       Rcpp::Named("oid")=extern_oid,
                                       Rcpp::Named("ovar")=extern_ovar,
                                       Rcpp::Named("oc")=extern_oc,
                                       Rcpp::Named("otheta")=extern_otheta,
                                       Rcpp::Named("sts")=extern_sts,
                                       Rcpp::Named("sid")=extern_sid,
                                       Rcpp::Named("svar")=extern_svar,
                                       Rcpp::Named("sc")=extern_sc,
                                       Rcpp::Named("stheta")=extern_stheta,
Rcpp::Named("f.trees")=Rcpp::CharacterVector(mtrees.str()),
Rcpp::Named("s.trees")=Rcpp::CharacterVector(strees.str()));

   // summary statistics
   if(summarystats) {
      //COUT << "Calculating summary statistics" << endl;
      //unsigned int varcount[p];
      //std::vector<unsigned int> varcount(p);
/*
      Rcpp::IntegerVector varcount(p);
      for(size_t i=0;i<p;i++) varcount[i]=0;
      unsigned int tmaxd=0;
      unsigned int tmind=0;
      double tavgd=0.0;
*/

      ambm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
      tavgd/=(double)(nd*m);
      ret["mu.tavgd"]=tavgd;
      ret["mu.tmaxd"]=tmaxd;
      ret["mu.tmind"]=tmind;
      //Rcpp::NumericVector vc(p);
      //for(size_t i=0;i<p;i++) vc[i]=varcount[i];
      //ret["mu.varcount"]=vc;
      //ret["mu.varcount"]=Rcpp::clone(varcount);
      ret["mu.varcount"]=fvc;

      //for(size_t i=0;i<p;i++) varcount[i]=0;
      //tmaxd=0; tmind=0; tavgd=0.0;
      psbm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
      tavgd/=(double)(nd*mh);
      ret["sd.tavgd"]=tavgd;
      ret["sd.tmaxd"]=tmaxd;
      ret["sd.tmind"]=tmind;
      //Rcpp::NumericVector sdvc(p);
      //for(size_t i=0;i<p;i++) sdvc[i]=varcount[i];
      //ret["sd.varcount"]=sdvc;
      //ret["sd.varcount"]=varcount;
      ret["sd.varcount"]=svc;
   }

   if(r) delete [] r;
   PutRNGstate();

   return ret;
}


// Draw predictive realizations at the prediciton points, xp.
// Requires the original data locations, x, the number of
// mean and sd trees, m and mh, the number of draws
// saved from the MCMC, nd, the number of cutpoints used
// in the MCMC, numcut, the number of threads to run the
// predicitons in parallel, tc and the R output of the fitted
// MCMC, fit.
RcppExport SEXP cpsambrt_predict(
   SEXP _ix, //training points
   SEXP _ixp, //prediction points
   SEXP _im,  //number of trees in mean model
   SEXP _imh, //number of trees in variance model
   SEXP _ind, //number of draws saved from the posterior
   SEXP _ixicuts, //variable cutpoints for each predictor
   SEXP _itc, //number of parallel compute threads
   SEXP _ifit //saved fitted model object returned from cpsambrt
)
{
   //--------------------------------------------------
   //process args
   //xp prediction points
   Rcpp::NumericMatrix xpm(_ixp);
   size_t np = xpm.ncol();
   size_t p = xpm.nrow();
   double *xp = nullptr;
   if(np>0)  xp = &xpm[0];

   //x training points
   Rcpp::NumericMatrix xm(_ix);
   size_t n = xm.ncol();
   double *x = &xm[0];

   //make xinfo
   xinfo xi;
   Rcpp::List ixi(_ixicuts);
   xi.resize(p);
   for(size_t i=0;i<p;i++)
      xi[i]=Rcpp::as< std::vector<double> >(ixi[i]);

   //number of trees
   size_t m = Rcpp::as<int>(_im);
   size_t mh = Rcpp::as<int>(_imh);

   //nd
   size_t nd = Rcpp::as<int>(_ind);

   //thread count
   int tc = Rcpp::as<int>(_itc);

   // set up ambrt object
   ambrt ambm(m);
   ambm.settc(tc);  //set the number of threads when using OpenMP, etc.
   ambm.setxi(&xi); //set the cutpoints for this model object

   //setup psbrt object
   psbrt psbm(mh);
   psbm.settc(tc);  //set the number of threads when using OpenMP, etc.
   psbm.setxi(&xi); //set the cutpoints for this model object

   //cpsambrt fitted model object (contains posterior MCMC samples)
   Rcpp::List fit(_ifit);
   Rcpp::XPtr< std::vector<int> > e_ots(Rcpp::as<SEXP>(fit["ots"]));
   Rcpp::XPtr< std::vector<int> > e_oid(Rcpp::as<SEXP>(fit["oid"]));
   Rcpp::XPtr< std::vector<int> > e_ovar(Rcpp::as<SEXP>(fit["ovar"]));
   Rcpp::XPtr< std::vector<int> > e_oc(Rcpp::as<SEXP>(fit["oc"]));
   Rcpp::XPtr< std::vector<double> > e_otheta(Rcpp::as<SEXP>(fit["otheta"]));
   Rcpp::XPtr< std::vector<int> > e_sts(Rcpp::as<SEXP>(fit["sts"]));
   Rcpp::XPtr< std::vector<int> > e_sid(Rcpp::as<SEXP>(fit["sid"]));
   Rcpp::XPtr< std::vector<int> > e_svar(Rcpp::as<SEXP>(fit["svar"]));
   Rcpp::XPtr< std::vector<int> > e_sc(Rcpp::as<SEXP>(fit["sc"]));
   Rcpp::XPtr< std::vector<double> > e_stheta(Rcpp::as<SEXP>(fit["stheta"]));

   //objects where we'll store the realizations
   Rcpp::NumericMatrix trdraw(nd,n);
   Rcpp::NumericMatrix tedraw(nd,np);
   Rcpp::NumericMatrix trdrawh(nd,n);
   Rcpp::NumericMatrix tedrawh(nd,np);
   double *f = new double[n], *fp = new double[np];
   dinfo di, dip;
   di.x = x;   di.y=f;   di.p = p;  di.n=n;   di.tc=tc;
   dip.x = xp; dip.y=fp; dip.p = p; dip.n=np; dip.tc=tc;

   // Temporary vectors used for loading one model realization at a time.
   std::vector<int> onn(m,1);
   std::vector<std::vector<int> > oid(m, std::vector<int>(1));
   std::vector<std::vector<int> > ov(m, std::vector<int>(1));
   std::vector<std::vector<int> > oc(m, std::vector<int>(1));
   std::vector<std::vector<double> > otheta(m, std::vector<double>(1));
   std::vector<int> snn(mh,1);
   std::vector<std::vector<int> > sid(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sv(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sc(mh, std::vector<int>(1));
   std::vector<std::vector<double> > stheta(mh, std::vector<double>(1));

   // Draw realizations of the posterior predictive.
   size_t curdx=0;
   size_t cumdx=0;
   // Mean trees first
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<m;j++) {
         onn[j]=e_ots->at(i*m+j);
         oid[j].resize(onn[j]);
         ov[j].resize(onn[j]);
         oc[j].resize(onn[j]);
         otheta[j].resize(onn[j]);
         for(int k=0;k<onn[j];k++) {
            oid[j][k]=e_oid->at(cumdx+curdx+k);
            ov[j][k]=e_ovar->at(cumdx+curdx+k);
            oc[j][k]=e_oc->at(cumdx+curdx+k);
            otheta[j][k]=e_otheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)onn[j];
      }
      cumdx+=curdx;

      ambm.loadtree(0,m,onn,oid,ov,oc,otheta);
      // draw realization
      ambm.predict(&di);
      for(size_t j=0;j<n;j++) trdraw(i,j) = f[j];
      ambm.predict(&dip);
      for(size_t j=0;j<np;j++) tedraw(i,j) = fp[j];
   }
   // Variance trees second
   cumdx=0;
   curdx=0;
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<mh;j++) {
         snn[j]=e_sts->at(i*mh+j);
         sid[j].resize(snn[j]);
         sv[j].resize(snn[j]);
         sc[j].resize(snn[j]);
         stheta[j].resize(snn[j]);
         for(int k=0;k<snn[j];k++) {
            sid[j][k]=e_sid->at(cumdx+curdx+k);
            sv[j][k]=e_svar->at(cumdx+curdx+k);
            sc[j][k]=e_sc->at(cumdx+curdx+k);
            stheta[j][k]=e_stheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)snn[j];
      }
      cumdx+=curdx;

      psbm.loadtree(0,mh,snn,sid,sv,sc,stheta);
      // draw realization
      psbm.predict(&di);
      for(size_t j=0;j<n;j++) trdrawh(i,j) = f[j];
      psbm.predict(&dip);
      for(size_t j=0;j<np;j++) tedrawh(i,j) = fp[j];
   }

   // Save the draws and return to R.
   Rcpp::List ret;
   ret["f.train"]=trdraw;
   ret["s.train"]=trdrawh;
   if(np>0) {
     ret["f.test"]=tedraw;
     ret["s.test"]=tedrawh;
   }

   if(f) delete [] f;
   if(fp) delete [] fp;
   return ret;
}



// Calculate the posterior variable activity realizations.
// Requires the training points, x,
// the number of mean and sd trees, m and mh,
// the number of draws saved from the MCMC, nd, 
// the number of threads to use in performing the calculations, tc,
// and the R output of the fitted MCMC, fit.
RcppExport SEXP cpsambrt_vartivity(
   SEXP _ip, //number of predictor variables
   SEXP _im,  //number of trees in mean model
   SEXP _imh, //number of trees in variance model
   SEXP _ind, //number of draws saved from the posterior
   SEXP _itc, //number of parallel compute threads
   SEXP _ifit //saved fitted model object returned from cpsambrt
)
{
   //--------------------------------------------------
   //process args
   //number of predictor variables
   size_t p = Rcpp::as<int>(_ip);
  
   //number of trees
   size_t m = Rcpp::as<int>(_im);
   size_t mh = Rcpp::as<int>(_imh);

   //nd
   size_t nd = Rcpp::as<int>(_ind);

   //thread count
   int tc = Rcpp::as<int>(_itc);
   COUT << "In vartivity, thread cout is " << tc << "\n";

   //cpsambrt fitted model object (contains posterior MCMC samples)
   Rcpp::List fit(_ifit);
   Rcpp::XPtr< std::vector<int> > e_ots(Rcpp::as<SEXP>(fit["ots"]));
   Rcpp::XPtr< std::vector<int> > e_oid(Rcpp::as<SEXP>(fit["oid"]));
   Rcpp::XPtr< std::vector<int> > e_ovar(Rcpp::as<SEXP>(fit["ovar"]));
   Rcpp::XPtr< std::vector<int> > e_oc(Rcpp::as<SEXP>(fit["oc"]));
   Rcpp::XPtr< std::vector<double> > e_otheta(Rcpp::as<SEXP>(fit["otheta"]));
   Rcpp::XPtr< std::vector<int> > e_sts(Rcpp::as<SEXP>(fit["sts"]));
   Rcpp::XPtr< std::vector<int> > e_sid(Rcpp::as<SEXP>(fit["sid"]));
   Rcpp::XPtr< std::vector<int> > e_svar(Rcpp::as<SEXP>(fit["svar"]));
   Rcpp::XPtr< std::vector<int> > e_sc(Rcpp::as<SEXP>(fit["sc"]));
   Rcpp::XPtr< std::vector<double> > e_stheta(Rcpp::as<SEXP>(fit["stheta"]));

   //objects where we'll store the realizations
   Rcpp::NumericMatrix vdraw(nd,p);
   Rcpp::NumericMatrix vdrawh(nd,p);
   for(size_t i=0;i<nd;i++)
      for(size_t j=0;j<p;j++) {
         vdraw(i,j)=0.0;
         vdrawh(i,j)=0.0;
      }

   // Temporary vectors used for loading one model realization at a time.
   std::vector<int> onn(m,1);
   std::vector<std::vector<int> > oid(m, std::vector<int>(1));
   std::vector<std::vector<int> > ov(m, std::vector<int>(1));
   std::vector<std::vector<int> > oc(m, std::vector<int>(1));
   std::vector<std::vector<double> > otheta(m, std::vector<double>(1));
   std::vector<int> snn(mh,1);
   std::vector<std::vector<int> > sid(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sv(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sc(mh, std::vector<int>(1));
   std::vector<std::vector<double> > stheta(mh, std::vector<double>(1));

   // Draw realizations of the posterior predictive.
   size_t curdx=0;
   size_t cumdx=0;
   size_t cid=0;
   bool haschild=false;
   // Mean trees first
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<m;j++) {
         onn[j]=e_ots->at(i*m+j);
         oid[j].resize(onn[j]);
         ov[j].resize(onn[j]);
         oc[j].resize(onn[j]);
         otheta[j].resize(onn[j]);
         for(int k=0;k<onn[j];k++) {
            oid[j][k]=e_oid->at(cumdx+curdx+k);
            ov[j][k]=e_ovar->at(cumdx+curdx+k);
            oc[j][k]=e_oc->at(cumdx+curdx+k);
            otheta[j][k]=e_otheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)onn[j];
      }
      cumdx+=curdx;

      for(size_t j=0;j<m;j++)
         for(int k=0;k<onn[j];k++) {
            cid=2*oid[j][k];
            for(int l=0;l<onn[j];l++)
               if(oid[j][l]== (int)cid)
                  haschild=true;

            if(haschild) vdraw(i,ov[j][k])++;
            haschild=false;
         }
   }

   // Variance trees second
   cumdx=0;
   curdx=0;
   cid=0;
   haschild=false;
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<mh;j++) {
         snn[j]=e_sts->at(i*mh+j);
         sid[j].resize(snn[j]);
         sv[j].resize(snn[j]);
         sc[j].resize(snn[j]);
         stheta[j].resize(snn[j]);
         for(int k=0;k<snn[j];k++) {
            sid[j][k]=e_sid->at(cumdx+curdx+k);
            sv[j][k]=e_svar->at(cumdx+curdx+k);
            sc[j][k]=e_sc->at(cumdx+curdx+k);
            stheta[j][k]=e_stheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)snn[j];
      }
      cumdx+=curdx;

      for(size_t j=0;j<mh;j++)
         for(int k=0;k<snn[j];k++) {
            cid=2*sid[j][k];
            for(int l=0;l<snn[j];l++)
               if(sid[j][l]== (int)cid)
                  haschild=true;

            if(haschild) vdrawh(i,sv[j][k])++;
            haschild=false;
         }
   }

   // Save the draws and return to R.
   Rcpp::List ret;
   ret["vdraws"]=vdraw;
   ret["vdrawsh"]=vdrawh;

   return ret;
}


// Save BART posterior model realizations.
// Requires the original data locations, x, the number of
// mean and sd trees, m and mh, the number of draws
// saved from the MCMC, nd, the number of cutpoints used
// in the MCMC, numcut, the number of threads to run the
// predicitons in parallel, tc and the R output of the fitted
// MCMC, fit.
RcppExport SEXP cpsambrt_save(
   SEXP _ifilename, //filename
   SEXP _ix, //training points
   SEXP _im,  //number of trees in mean model
   SEXP _imh, //number of trees in variance model
   SEXP _ind, //number of draws saved from the posterior
   SEXP _ixicuts, //variable cutpoints for each predictor
   SEXP _ifit //saved fitted model object returned from cpsambrt
)
{
   //--------------------------------------------------
   //process args
   // filename
   std::string fname = Rcpp::as<std::string>(_ifilename);

   //x training points
   Rcpp::NumericMatrix xm(_ix);
   size_t n = xm.ncol();
   size_t p = xm.nrow();
   //double *x = &xm[0];

   //make xinfo
   xinfo xi;
   Rcpp::List ixi(_ixicuts);
   xi.resize(p);
   for(size_t i=0;i<p;i++)
      xi[i]=Rcpp::as< std::vector<double> >(ixi[i]);

   //number of trees
   size_t m = Rcpp::as<int>(_im);
   size_t mh = Rcpp::as<int>(_imh);

   //nd
   size_t nd = Rcpp::as<int>(_ind);

   // set up ambrt object
   ambrt ambm(m);
   ambm.settc(1);  //set the number of threads when using OpenMP, etc.
   ambm.setxi(&xi); //set the cutpoints for this model object

   //setup psbrt object
   psbrt psbm(mh);
   psbm.settc(1);  //set the number of threads when using OpenMP, etc.
   psbm.setxi(&xi); //set the cutpoints for this model object

   //cpsambrt fitted model object (contains posterior MCMC samples)
   Rcpp::List fit(_ifit);
   Rcpp::XPtr< std::vector<int> > e_ots(Rcpp::as<SEXP>(fit["ots"]));
   Rcpp::XPtr< std::vector<int> > e_oid(Rcpp::as<SEXP>(fit["oid"]));
   Rcpp::XPtr< std::vector<int> > e_ovar(Rcpp::as<SEXP>(fit["ovar"]));
   Rcpp::XPtr< std::vector<int> > e_oc(Rcpp::as<SEXP>(fit["oc"]));
   Rcpp::XPtr< std::vector<double> > e_otheta(Rcpp::as<SEXP>(fit["otheta"]));
   Rcpp::XPtr< std::vector<int> > e_sts(Rcpp::as<SEXP>(fit["sts"]));
   Rcpp::XPtr< std::vector<int> > e_sid(Rcpp::as<SEXP>(fit["sid"]));
   Rcpp::XPtr< std::vector<int> > e_svar(Rcpp::as<SEXP>(fit["svar"]));
   Rcpp::XPtr< std::vector<int> > e_sc(Rcpp::as<SEXP>(fit["sc"]));
   Rcpp::XPtr< std::vector<double> > e_stheta(Rcpp::as<SEXP>(fit["stheta"]));

   // Temporary vectors used for loading one model realization at a time.
   std::vector<int> onn(m,1);
   std::vector<std::vector<int> > oid(m, std::vector<int>(1));
   std::vector<std::vector<int> > ov(m, std::vector<int>(1));
   std::vector<std::vector<int> > oc(m, std::vector<int>(1));
   std::vector<std::vector<double> > otheta(m, std::vector<double>(1));
   std::vector<int> snn(mh,1);
   std::vector<std::vector<int> > sid(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sv(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sc(mh, std::vector<int>(1));
   std::vector<std::vector<double> > stheta(mh, std::vector<double>(1));

   // Save fitted model to the output file.  The format is as follows:
   // n<endl>
   // p<endl>
   // m<endl>
   // mh<endl>
   // numcut<endl>
   // xinfo<endl>
   // draw1 meantree1<endl>
   // ...
   // draw1 meantreem<endl>
   // ..
   // drawnd vartreemh<endl>
   std::ofstream treeoutf(fname.c_str()); //file to write tree to.
   treeoutf << n << endl;
   treeoutf << p << endl;
   treeoutf << m << endl;
   treeoutf << mh << endl;
   for(size_t i=0;i<p;i++)
      treeoutf << xi[i].size() << endl;  //numcuts for each predictor
   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<xi[i].size();j++)
         treeoutf << xi[i][j] << " ";
      treeoutf << endl;
   }


   // Save realizations of the posterior predictive.
   size_t curdx=0;
   size_t cumdx=0;
   COUT << "Mean trees...";
   // Mean trees first
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<m;j++) {
         onn[j]=e_ots->at(i*m+j);
         oid[j].resize(onn[j]);
         ov[j].resize(onn[j]);
         oc[j].resize(onn[j]);
         otheta[j].resize(onn[j]);
         for(int k=0;k<onn[j];k++) {
            oid[j][k]=e_oid->at(cumdx+curdx+k);
            ov[j][k]=e_ovar->at(cumdx+curdx+k);
            oc[j][k]=e_oc->at(cumdx+curdx+k);
            otheta[j][k]=e_otheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)onn[j];
      }
      cumdx+=curdx;

      ambm.loadtree(0,m,onn,oid,ov,oc,otheta);
      // save these m trees
      for(size_t j=0;j<m;j++)
         treeoutf << *ambm.gettree(j);
   }
   // Variance trees second
   COUT << " Variance trees... ";
   cumdx=0;
   curdx=0;
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<mh;j++) {
         snn[j]=e_sts->at(i*mh+j);
         sid[j].resize(snn[j]);
         sv[j].resize(snn[j]);
         sc[j].resize(snn[j]);
         stheta[j].resize(snn[j]);
         for(int k=0;k<snn[j];k++) {
            sid[j][k]=e_sid->at(cumdx+curdx+k);
            sv[j][k]=e_svar->at(cumdx+curdx+k);
            sc[j][k]=e_sc->at(cumdx+curdx+k);
            stheta[j][k]=e_stheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)snn[j];
      }
      cumdx+=curdx;

      psbm.loadtree(0,mh,snn,sid,sv,sc,stheta);
      // draw these mh trees
      for(size_t j=0;j<mh;j++)
         treeoutf << *psbm.gettree(j);
   }

   treeoutf.close();
   COUT << "done." << endl;

   Rcpp::List ret;
   ret["saved"]=fname;

   return ret;
}

// Save BART posterior model realizations.
// Requires the original data locations, x, the number of
// mean and sd trees, m and mh, the number of draws
// saved from the MCMC, nd, the number of cutpoints used
// in the MCMC, numcut, the number of threads to run the
// predicitons in parallel, tc and the R output of the fitted
// MCMC, fit.
RcppExport SEXP cpsambrt_Rexport(
   SEXP _ix, //training points
   SEXP _im,  //number of trees in mean model
   SEXP _imh, //number of trees in variance model
   SEXP _ind, //number of draws saved from the posterior
   SEXP _ixicuts, //variable cutpoints for each predictor
   SEXP _ifit //saved fitted model object returned from cpsambrt
)
{
   //--------------------------------------------------
   //process args

   //x training points
   Rcpp::NumericMatrix xm(_ix);
   //size_t n = xm.ncol();
   size_t p = xm.nrow();
   //double *x = &xm[0];

   //make xinfo
   xinfo xi;
   Rcpp::List ixi(_ixicuts);
   xi.resize(p);
   for(size_t i=0;i<p;i++)
      xi[i]=Rcpp::as< std::vector<double> >(ixi[i]);

   //number of trees
   size_t m = Rcpp::as<int>(_im);
   size_t mh = Rcpp::as<int>(_imh);

   //nd
   size_t nd = Rcpp::as<int>(_ind);

   // set up ambrt object
   ambrt ambm(m);
   ambm.settc(1);  //set the number of threads when using OpenMP, etc.
   ambm.setxi(&xi); //set the cutpoints for this model object

   //setup psbrt object
   psbrt psbm(mh);
   psbm.settc(1);  //set the number of threads when using OpenMP, etc.
   psbm.setxi(&xi); //set the cutpoints for this model object

   //cpsambrt fitted model object (contains posterior MCMC samples)
   Rcpp::List fit(_ifit);
   Rcpp::XPtr< std::vector<int> > e_ots(Rcpp::as<SEXP>(fit["ots"]));
   Rcpp::XPtr< std::vector<int> > e_oid(Rcpp::as<SEXP>(fit["oid"]));
   Rcpp::XPtr< std::vector<int> > e_ovar(Rcpp::as<SEXP>(fit["ovar"]));
   Rcpp::XPtr< std::vector<int> > e_oc(Rcpp::as<SEXP>(fit["oc"]));
   Rcpp::XPtr< std::vector<double> > e_otheta(Rcpp::as<SEXP>(fit["otheta"]));
   Rcpp::XPtr< std::vector<int> > e_sts(Rcpp::as<SEXP>(fit["sts"]));
   Rcpp::XPtr< std::vector<int> > e_sid(Rcpp::as<SEXP>(fit["sid"]));
   Rcpp::XPtr< std::vector<int> > e_svar(Rcpp::as<SEXP>(fit["svar"]));
   Rcpp::XPtr< std::vector<int> > e_sc(Rcpp::as<SEXP>(fit["sc"]));
   Rcpp::XPtr< std::vector<double> > e_stheta(Rcpp::as<SEXP>(fit["stheta"]));

   // Temporary vectors used for loading one model realization at a time.
   std::vector<int> onn(m,1);
   std::vector<std::vector<int> > oid(m, std::vector<int>(1));
   std::vector<std::vector<int> > ov(m, std::vector<int>(1));
   std::vector<std::vector<int> > oc(m, std::vector<int>(1));
   std::vector<std::vector<double> > otheta(m, std::vector<double>(1));
   std::vector<int> snn(mh,1);
   std::vector<std::vector<int> > sid(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sv(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sc(mh, std::vector<int>(1));
   std::vector<std::vector<double> > stheta(mh, std::vector<double>(1));

   // Save fitted model to the output stringstream.  The format is as follows:
   // draw1 meantree1<endl>
   // ...
   // draw1 meantreem<endl>
   // ..
   // drawnd vartreemh<endl>

   std::stringbuf outbuf;
   std::ostream soutbuf(&outbuf);

   // Important: need to set the precision of the stringstream so that 
   // it will record doubles accurately.  Otherwise it will round, introducing
   // a noticeable amount of error when we subsequently import saved trees.
   soutbuf << std::setprecision (std::numeric_limits<double>::digits10 + 2);

   // Save realizations of the posterior predictive.
   size_t curdx=0;
   size_t cumdx=0;
   // Mean trees first
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<m;j++) {
         onn[j]=e_ots->at(i*m+j);
         oid[j].resize(onn[j]);
         ov[j].resize(onn[j]);
         oc[j].resize(onn[j]);
         otheta[j].resize(onn[j]);
         for(int k=0;k<onn[j];k++) {
            oid[j][k]=e_oid->at(cumdx+curdx+k);
            ov[j][k]=e_ovar->at(cumdx+curdx+k);
            oc[j][k]=e_oc->at(cumdx+curdx+k);
            otheta[j][k]=e_otheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)onn[j];
      }
      cumdx+=curdx;

      ambm.loadtree(0,m,onn,oid,ov,oc,otheta);
      // save these m trees
      for(size_t j=0;j<m;j++)
         soutbuf << *ambm.gettree(j);
   }
   // Variance trees second
   cumdx=0;
   curdx=0;
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<mh;j++) {
         snn[j]=e_sts->at(i*mh+j);
         sid[j].resize(snn[j]);
         sv[j].resize(snn[j]);
         sc[j].resize(snn[j]);
         stheta[j].resize(snn[j]);
         for(int k=0;k<snn[j];k++) {
            sid[j][k]=e_sid->at(cumdx+curdx+k);
            sv[j][k]=e_svar->at(cumdx+curdx+k);
            sc[j][k]=e_sc->at(cumdx+curdx+k);
            stheta[j][k]=e_stheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)snn[j];
      }
      cumdx+=curdx;

      psbm.loadtree(0,mh,snn,sid,sv,sc,stheta);
      // draw these mh trees
      for(size_t j=0;j<mh;j++)
         soutbuf << *psbm.gettree(j);
   }

   Rcpp::List ret;
   ret["poststr"]=outbuf.str();

   return ret;
}


// Load BART posterior model realizations that were previously
// saved using cpsambrt_Rexport.
// Input is the string saved from cpsambrt_Rexport.
RcppExport SEXP cpsambrt_Rimport(
   SEXP _ix, //training points
   SEXP _im,  //number of trees in mean model
   SEXP _imh, //number of trees in variance model
   SEXP _ind, //number of draws saved from the posterior
   SEXP _ixicuts, //variable cutpoints for each predictor
   SEXP _ipoststr //saved posterior tree draws
)
{
   //--------------------------------------------------
   //process args
   std::stringbuf inbuf(Rcpp::as<std::string>(_ipoststr));
   std::istream sinbuf(&inbuf);

   //x training points
   Rcpp::NumericMatrix xm(_ix);
   //size_t n = xm.ncol();
   size_t p = xm.nrow();

   //make xinfo
   xinfo xi;
   Rcpp::List ixi(_ixicuts);
   xi.resize(p);
   for(size_t i=0;i<p;i++)
      xi[i]=Rcpp::as< std::vector<double> >(ixi[i]);

   //number of trees
   size_t m = Rcpp::as<int>(_im);
   size_t mh = Rcpp::as<int>(_imh);

   //nd
   size_t nd = Rcpp::as<int>(_ind);

   // set up ambrt object
   ambrt ambm(m);
   ambm.settc(1);  //set the number of threads when using OpenMP, etc.
   ambm.setxi(&xi); //set the cutpoints for this model object

   //setup psbrt object
   psbrt psbm(mh);
   psbm.settc(1);  //set the number of threads when using OpenMP, etc.
   psbm.setxi(&xi); //set the cutpoints for this model object

   // Temporary vectors used for loading posterior model realizations.
   std::vector<int> onn(nd*m,1);
   std::vector<std::vector<int> > oid(nd*m, std::vector<int>(1));
   std::vector<std::vector<int> > ov(nd*m, std::vector<int>(1));
   std::vector<std::vector<int> > oc(nd*m, std::vector<int>(1));
   std::vector<std::vector<double> > otheta(nd*m, std::vector<double>(1));
   std::vector<int> snn(nd*mh,1);
   std::vector<std::vector<int> > sid(nd*mh, std::vector<int>(1));
   std::vector<std::vector<int> > sv(nd*mh, std::vector<int>(1));
   std::vector<std::vector<int> > sc(nd*mh, std::vector<int>(1));
   std::vector<std::vector<double> > stheta(nd*mh, std::vector<double>(1));

   // Load realizations of the posterior.
   // Mean trees first
   for(size_t i=0;i<nd;i++) {
      for(size_t j=0;j<m;j++) {
         sinbuf >> *ambm.gettree(j);

      }
      ambm.savetree(i,m,onn,oid,ov,oc,otheta);
   }
   // Variance trees second
   for(size_t i=0;i<nd;i++) {
      for(size_t j=0;j<mh;j++) {
         sinbuf >> *psbm.gettree(j);
      }
      psbm.savetree(i,mh,snn,sid,sv,sc,stheta);
   }

   //Flatten posterior trees to a few (very long) vectors so we can just pass pointers
   //to these vectors back to R (which is much much faster than copying all the data back).
   std::vector<int>* e_ots=new std::vector<int>(nd*m);
   std::vector<int>* e_oid=new std::vector<int>;
   std::vector<int>* e_ovar=new std::vector<int>;
   std::vector<int>* e_oc=new std::vector<int>;
   std::vector<double>* e_otheta=new std::vector<double>;
   std::vector<int>* e_sts=new std::vector<int>(nd*mh);
   std::vector<int>* e_sid=new std::vector<int>;
   std::vector<int>* e_svar=new std::vector<int>;
   std::vector<int>* e_sc=new std::vector<int>;
   std::vector<double>* e_stheta=new std::vector<double>;
   for(size_t i=0;i<nd;i++)
      for(size_t j=0;j<m;j++) {
         e_ots->at(i*m+j)=static_cast<int>(oid[i*m+j].size());
         e_oid->insert(e_oid->end(),oid[i*m+j].begin(),oid[i*m+j].end());
         e_ovar->insert(e_ovar->end(),ov[i*m+j].begin(),ov[i*m+j].end());
         e_oc->insert(e_oc->end(),oc[i*m+j].begin(),oc[i*m+j].end());
         e_otheta->insert(e_otheta->end(),otheta[i*m+j].begin(),otheta[i*m+j].end());
      }
   for(size_t i=0;i<nd;i++)
      for(size_t j=0;j<mh;j++) {
         e_sts->at(i*mh+j)=static_cast<int>(sid[i*mh+j].size());
         e_sid->insert(e_sid->end(),sid[i*mh+j].begin(),sid[i*mh+j].end());
         e_svar->insert(e_svar->end(),sv[i*mh+j].begin(),sv[i*mh+j].end());
         e_sc->insert(e_sc->end(),sc[i*mh+j].begin(),sc[i*mh+j].end());
         e_stheta->insert(e_stheta->end(),stheta[i*mh+j].begin(),stheta[i*mh+j].end());
      }

   Rcpp::XPtr< std::vector<int> > extern_ots(e_ots,true);
   Rcpp::XPtr< std::vector<int> > extern_oid(e_oid,true);
   Rcpp::XPtr< std::vector<int> > extern_ovar(e_ovar,true);
   Rcpp::XPtr< std::vector<int> > extern_oc(e_oc,true);
   Rcpp::XPtr< std::vector<double> > extern_otheta(e_otheta,true);
   Rcpp::XPtr< std::vector<int> > extern_sts(e_sts,true);
   Rcpp::XPtr< std::vector<int> > extern_sid(e_sid,true);
   Rcpp::XPtr< std::vector<int> > extern_svar(e_svar,true);
   Rcpp::XPtr< std::vector<int> > extern_sc(e_sc,true);
   Rcpp::XPtr< std::vector<double> > extern_stheta(e_stheta,true);

   Rcpp::List ret = Rcpp::List::create(Rcpp::Named("ots")=extern_ots,
                                       Rcpp::Named("oid")=extern_oid,
                                       Rcpp::Named("ovar")=extern_ovar,
                                       Rcpp::Named("oc")=extern_oc,
                                       Rcpp::Named("otheta")=extern_otheta,
                                       Rcpp::Named("sts")=extern_sts,
                                       Rcpp::Named("sid")=extern_sid,
                                       Rcpp::Named("svar")=extern_svar,
                                       Rcpp::Named("sc")=extern_sc,
                                       Rcpp::Named("stheta")=extern_stheta);


   return ret;
}
