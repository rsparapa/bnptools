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

#include "rrn.h"
#include "Dp.h"
#include "DpMuTau.h"
#include "DpMuSigma.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"

double lreg(int n, double *x, double *y, double sigma, double betabar, double Abeta, rn& gen);

RcppExport SEXP cbiv(
   SEXP _zx,
   SEXP _x,
   SEXP _T,
   SEXP _Y,
   SEXP _burn,
   SEXP _nd,
   SEXP _burnf,
   SEXP _burnh,
   SEXP _m,
   SEXP _nc,
   SEXP _power,
   SEXP _base,
   SEXP _tauf,
   SEXP _tauh,
   SEXP _betabar,
   SEXP _Abeta,
   SEXP _v,
   SEXP _nu,
   SEXP _a,
   SEXP _ag,
   SEXP _priag,
   SEXP _centermeans,
   SEXP _fs,
   SEXP _hs,
   SEXP _betas,
   SEXP _mTs,
   SEXP _mYs,
   SEXP _sTs,
   SEXP _gammas,
   SEXP _sYs,
   SEXP _printevery
)
{
   Rprintf("*****************************************************************\n");
   Rprintf("*****Into main of cbiv\n");

   //-----------------------------------------------------------
   //random number generation
   GetRNGstate();
   rrn gen;

   //--------------------------------------------------
   //process arguments

   //zx
   Rcpp::NumericMatrix zxm(_zx);
   size_t nzx = zxm.ncol();
   size_t pzx = zxm.nrow();
   double *zx = &zxm[0];

   //x
   Rcpp::NumericMatrix xm(_x);
   size_t nx = xm.ncol();
   size_t px = xm.nrow();
   double *x = &xm[0];

   // T,Y ------------
   Rcpp::NumericVector Tv(_T);
   double *T = &Tv[0];
   size_t nT = Tv.size();

   Rcpp::NumericVector Yv(_Y);
   double *Y = &Yv[0];
   size_t nY = Yv.size();

   // burn, nd
   size_t burn = Rcpp::as<size_t>(_burn);
   size_t nd = Rcpp::as<size_t>(_nd);
   size_t burnf = Rcpp::as<size_t>(_burnf);
   size_t burnh = Rcpp::as<size_t>(_burnh);

   //bart prior
   size_t m = Rcpp::as<size_t>(_m);
   size_t nc = Rcpp::as<size_t>(_nc);
   double mybeta = Rcpp::as<double>(_power);
   double alpha = Rcpp::as<double>(_base);
   double tauf = Rcpp::as<double>(_tauf);
   double tauh = Rcpp::as<double>(_tauh);

   //beta prior
   double betabar = Rcpp::as<double>(_betabar);
   double Abeta = Rcpp::as<double>(_Abeta);

   //base prior -------------
   double v = Rcpp::as<double>(_v);
   double nu = Rcpp::as<double>(_nu);
   double a = Rcpp::as<double>(_a);

   //alpha prior--------
   Rcpp::NumericVector agv(_ag);
   Rcpp::NumericVector priagv(_priag);
   size_t nag=agv.size();
   Dp::dv ag(nag);
   Dp::dv priag(nag);
   for(size_t i=0;i<nag;i++) { // should use STL
      ag[i]=agv[i];
      priag[i]=priagv[i];
   }

   //should the means be centered after each DpMuSigma draw
   bool centermeans = Rcpp::as<bool>(_centermeans);

   //starting values
   //fs,hs
   Rcpp::NumericVector fsv(_fs);
   double *fs = &fsv[0];
   Rcpp::NumericVector hsv(_hs);
   double *hs = &hsv[0];

   double betas = Rcpp::as<double>(_betas);
   double mTs = Rcpp::as<double>(_mTs);
   double mYs = Rcpp::as<double>(_mYs);
   double sTs = Rcpp::as<double>(_sTs);
   double gammas = Rcpp::as<double>(_gammas);
   double sYs = Rcpp::as<double>(_sYs);

   //other
   size_t printevery = Rcpp::as<size_t>(_printevery);

   size_t n = nT;

   //--------------------------------------------------
   // print args

   Rprintf("***burn, nd, burnf, burnh: %ld, %ld, %ld\n",burn,nd,burnf,burnh);

   Rprintf("*** Prior:\n");
   Rprintf("m (num trees), nc (num cut points), %ld, %ld\n",m,nc);
   Rprintf("****power, base, tauf, tauh: %lf, %lf, %lf, %lf\n",mybeta,alpha,tauf,tauh);
   Rprintf("betabar: %lf\n",betabar);
   Rprintf("Abeta: %lf\n",Abeta);
   Rprintf("v: %lf\n",v);
   Rprintf("nu: %lf\n",nu);
   Rprintf("a: %lf\n",a);
   Rprintf("alpha prior: grid size: %ld\n",ag.size());
   Rprintf("\tag, first, last: %lf, %lf\n",ag[0],ag[nag-1]);
   Rprintf("\tpriag, first, last: %lf, %lf\n",priag[0],priag[nag-1]);

   Rprintf("*** Data:\n");
   Rprintf("nzx,pzx: %ld, %ld\n",nzx,pzx);
   Rprintf("nx,px: %ld, %ld\n",nx,px);
   Rprintf("nT: %ld\n",nT);
   Rprintf("first and last T: %lf, %lf\n",T[0],T[nT-1]);
   Rprintf("nY: %ld\n",nY);
   Rprintf("first and last Y: %lf, %lf\n",Y[0],Y[nY-1]);

   Rprintf("*** starting values, n is: %ld\n",n);
   Rprintf("\t first and last fs: %lf, %lf\n",fs[0],fs[n-1]);
   Rprintf("\t first and last hs: %lf, %lf\n",hs[0],hs[n-1]);
   Rprintf(" starting for beta, mT, mY, sT, gamma, sY: %lf, %lf, %lf, %lf, %lf, %lf\n",betas,mTs,mYs,sTs,gammas,sYs);

   Rprintf("***other\n");
   Rprintf("printevery: %ld\n",printevery);



   //-------------------------------------------------
   //-------------------------------------------------
   // bart f setup
   //--------------------------------------------------
   // first make zx2 which is zx with the columns repeated
   cout << "@@@@@@@@@@@@@@@@@write out zx2 to check\n";
   cout << "n, pzc: " << n << ", " << pzx << std::endl;
   double *zx2 = new double[n*pzx*2];
   size_t start=0;
   for(size_t i=0;i<n;i++) {
      start = i*2*pzx;
      for(size_t j=0;j<pzx;j++) {
         zx2[start+j] = zx[i*pzx+j];
      }
      start += pzx;
      for(size_t j=0;j<pzx;j++) {
         zx2[start+j] = zx[i*pzx+j];
      }
   }
   //write out 6 rows of zx2 to check
   for(size_t i=0;i< 6;i++) {
      for(size_t j=0;j<pzx;j++) {
         cout << zx2[i*pzx+j] << " ";
      }
      cout << std::endl;
   }
   heterbart bmf(m);
   bmf.setprior(alpha,mybeta,tauf);
   double *ytempf = new double[n*2];  //y for h bart
   double *svecf = new double[n*2];   // sigma_i for h bart
   bmf.setdata(pzx,n*2,zx2,ytempf,nc);

   Rcpp::NumericMatrix dfburn(burnf,n); //h draws on train
   Dp::dv fhatb(n,0.0);
   if(burnf) {
   for(size_t i=0; i<n; i++) {
      ytempf[2*i] = T[i] - mTs;
      svecf[2*i] = sTs;
      ytempf[2*i+1] = T[i] - mTs;
      svecf[2*i+1] = sTs;
   }
   for(size_t i=0;i<burnf;i++) {
      if(i%printevery==0) Rprintf("burnf: done %d (out of %d)\n",i,burnf);
      bmf.draw(svecf,gen);
      for(size_t j=0;j<n;j++) dfburn(i,j) = bmf.f(2*j);
      for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j] + bmf.f(2*j);
   }
   for(size_t j=0;j<n;j++) fhatb[j] = fhatb[j]/burnf;
   }
   //bmf data from starting
   double R;
   for(size_t i=0; i<n; i++) {
      ytempf[2*i] = T[i] - mTs;
      svecf[2*i] = sTs;
      R = (betas*sTs+gammas)*(T[i]-mTs) - sTs*(Y[i] - mYs - hs[i] - betas*mTs);
      ytempf[2*i+1] = R/gammas;
      svecf[2*i+1] = (sTs*sYs)/gammas;
   }
   // bmf ouput storage
   Rcpp::NumericMatrix df(nd,n); //f draws on train

   //-------------------------------------------------
   //-------------------------------------------------
   // bart h setup
   //--------------------------------------------------
   heterbart bmh(m);
   bmh.setprior(alpha,mybeta,tauh);
   double *ytemp = new double[n];  //y for h bart
   double *svec = new double[n];   // sigma_i for h bart
   bmh.setdata(px,n,x,ytemp,nc);


   //h burn-in
   Rcpp::NumericMatrix dhburn(burnh,n); //h draws on train
   double ZZ1=0.0;
   for(size_t i=0;i<burnh;i++) {
      if(i%printevery==0) Rprintf("burnh: done %d (out of %d)\n",i,burnh);

      // h conditional -----------------------------------
      // update ----- 
      for(size_t j=0;j<n;j++) {
         ZZ1 = (T[j] - mTs - bmf.f(2*j))/sTs;
         ytemp[j] = Y[j] - mYs - betas*T[j] - gammas * ZZ1;
         svec[j] = sYs; 
      }
      //draw -----
      bmh.draw(svec,gen);

      for(size_t j=0;j<n;j++) {
         dhburn(i,j) = bmh.f(j);
      }
   }

   //bmh data from starting
   for(size_t i=0;i<n;i++) {
      ZZ1 = (T[i]-mTs-fs[i])/sTs;
      ytemp[i] = Y[i] - mYs - betas*T[i] - gammas * ZZ1;
      svec[i] = sYs; 
   }


   Rcpp::NumericMatrix dh(nd,n); //h draws on train


   //--------------------------------------------------
   //--------------------------------------------------
   // beta setup
   double *yb = new double[n];
   double *xb = new double[n];
   double Z1=0.0;
   for(size_t i=0;i<n;i++) {
      Z1 = (T[i]-mTs-fs[i])/sTs;
      yb[i] = (Y[i] - mYs - hs[i] - gammas * Z1)/sYs;
      xb[i] = T[i]/sYs;
   }
   double betad = betas;
   Rcpp::NumericVector dbeta(nd); //storage for output beta draws

   //--------------------------------------------------
   //--------------------------------------------------
   // DpMuSigma setup
   cout << "\n*** making a DpMuSigma object:\n";
   
   Dp::dv itheta(5,0);
   itheta[0] = mTs; itheta[1] = mYs; itheta[2] = sTs; itheta[3] = gammas; itheta[4] = sYs; 
   cout << "itheta:\n";
   printdv(itheta);

   //data for errors
   double *yS = new double[2*n];
   //intialize using starting values
   for(size_t i=0;i<n;i++) {
      yS[2*i] = T[i] - fs[i];
      yS[2*i+1] = Y[i] - betas*T[i] - hs[i];
   }
   cout << "check yS: " << yS[0] << ", " << yS[2*n-1] << std::endl;
   
   DpMuSigma dpmS(n,itheta);
   dpmS.setData(yS);
   dpmS.setPrior(v,nu,a);
   dpmS.setAlphaPrior(ag,priag);
   dpmS.setalpha(1.0);

   // to screen
   cout << "dpmS.toscreen()\n";
   dpmS.toscreen();

   //output storage for DpMuSigma
   Rcpp::IntegerVector dnpart(nd);
   Rcpp::NumericVector dalpha(nd);
   Rcpp::NumericMatrix dmu1(nd,n);
   Rcpp::NumericMatrix dsigma1(nd,n);
   Rcpp::NumericMatrix dmu2(nd,n);
   Rcpp::NumericMatrix dsigma2(nd,n);

   bool dobeta = true;
   bool dof = true;
   bool doh = true;
   bool dodp = true;

   bool doup = true;

   //--------------------------------------------------
   //MCMC
   time_t tp;
   int time1 = time(&tp);
   Dp::dvv tmat; 
   tmat = dpmS.thetaMatrix();

   //check tmat
   cout << "tmat.size: " << tmat.size() << std::endl;
   cout << "first 5 rows of initial tmat:\n";
   for(size_t i=0;i<5;i++) {
      for(size_t j=0 ; j< 5;j++) cout << tmat[i][j] << " ";
      cout << std::endl;
   }

   double mT=mTs,mY=mYs,sT=sTs,gamma=gammas,sY=sYs; //temporary values for mean and Lchol of Sigma
   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) Rprintf("done %d (out of %d)\n",i,nd+burn);

      // beta conditional -----------------------------------
      if(dobeta) {
      // update -----
      if(doup) {
      for(size_t j=0;j<n;j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         Z1 = (T[j] - mT - bmf.f(2*j))/sT;
         yb[j] = (Y[j] - mY - bmh.f(j) - gamma * Z1)/sY;
         xb[j] = T[j]/sY;
      }
      } else {
      for(size_t j=0;j<n;j++) {
         Z1 = (T[j] - mTs - fs[j])/sTs;
         yb[j] = (Y[j] - mYs - hs[j] - gammas * Z1)/sYs;
         xb[j] = T[j]/sYs;
      }
      }
      //draw -----
      betad = lreg(n,xb,yb,1.0,betabar,Abeta,gen);
      }

      // h conditional -----------------------------------
      if(doh) {
      // update ----- 
      if(doup) {
      for(size_t j=0;j<n;j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         ZZ1 = (T[j] - mT - bmf.f(2*j))/sT;
         ytemp[j] = Y[j] - mY - betad*T[j] - gamma * ZZ1;
         svec[j] = sY; 
      }
      } else {
      for(size_t j=0;j<n;j++) {
         ZZ1 = (T[j] - mTs - fs[j])/sTs;
         ytemp[j] = Y[j] - mYs - betas*T[j] - gammas * ZZ1;
         svec[j] = sYs; 
      }
      }
      //draw -----
      bmh.draw(svec,gen);
      }

      // f conditional --------------------------------
      if(dof) {
      // update -----
      if(doup) {
      for(size_t j=0; j<n; j++) {
         mT = tmat[j][0]; mY = tmat[j][1]; sT = tmat[j][2]; gamma = tmat[j][3]; sY = tmat[j][4];
         ytempf[2*j] = T[j] - mT;
         svecf[2*j] = sT;
         R = (betad*sT+gamma)*(T[j]-mT) - sT*(Y[j] - mY - bmh.f(j) - betad*mT);
         ytempf[2*j+1] = R/gamma;
         svecf[2*j+1] = (sT*sY)/gamma;
      }
      } else {
      for(size_t j=0; j<n; j++) {
         ytempf[2*j] = T[j] - mTs;
         svecf[2*j] = sTs;
         R = (betas*sTs+gammas)*(T[j]-mTs) - sTs*(Y[j] - mYs - hs[j] - betas*mTs);
         ytempf[2*j+1] = R/gammas;
         svecf[2*j+1] = (sTs*sYs)/gammas;
      }
      }
      // draws -----
      bmf.draw(svecf,gen);
      }

      // Sigma conditional -----------------------------
      //cout << "\n&&&&&&about to do Sigma\n";
      if(dodp) {
      //update -----
      for(size_t j=0;j<n;j++) {
         yS[2*j] = T[j] - bmf.f(2*j);
         yS[2*j+1] = Y[j] - betad*T[j] - bmh.f(j);
      }
      //draw -----
      dpmS.draw(gen);
      if(centermeans) dpmS.center();
      tmat = dpmS.thetaMatrix();
      }

      if(i >= burn) {
         size_t j=i-burn;

         if(dodp)
         dnpart[j] = dpmS.npart();
         dalpha[j] = dpmS.getalpha();
         for(size_t k=0;k<n;k++) {
            dmu1(j,k) = tmat[k][0];
            dsigma1(j,k) = tmat[k][2];
            dmu2(j,k) = tmat[k][1];
            dsigma2(j,k) = sqrt(tmat[k][3]*tmat[k][3] + tmat[k][4]*tmat[k][4]);
         }

         if(dobeta) {dbeta[j] = betad;}

         if(doh) {
         for(size_t k=0;k<n;k++) {
            dh(j,k) = bmh.f(k);
         }
         }

         if(dof) {
         for(size_t k=0;k<n;k++) {
            df(j,k) = bmf.f(2*k);
         }
         }
      }
   }
   int time2 = time(&tp);
   Rprintf("time: %d\n",time2-time1);


   //--------------------------------------------------
   // clear RNG at end
   PutRNGstate();

   Rcpp::List ret;
   ret["check"] = "biv";
   ret["dnpart"]=dnpart;
   ret["dalpha"]=dalpha;
   ret["dmu1"]=dmu1;
   ret["dsigma1"]=dsigma1;
   ret["dmu2"]=dmu2;
   ret["dsigma2"]=dsigma2;
   ret["dbeta"] = dbeta;
   ret["dh"] = dh;
   ret["df"] = df;
   ret["dfburn"] = dfburn;
   ret["dhburn"] = dhburn;

   //-------------------------------------------------
   // free
   if(yS) delete [] yS;
   if(yb) delete [] yb;
   if(xb) delete [] xb;
   if(ytemp) delete [] ytemp;
   if(svec) delete [] svec;
   if(zx2) delete [] zx2;
   if(ytempf) delete [] ytempf;
   if(svecf) delete [] svecf;

   return ret;
}
double lreg(int n, double *x, double *y, double sigma, double betabar, double Abeta, rn& gen)
{
   double sxx=0.0,sxy=0.0;
   for(int i=0;i<n;i++) {
      sxx += x[i]*x[i];
      sxy += x[i]*y[i];
   }
   double sig2 = sigma*sigma;
   double v = 1.0/((sxx/sig2) + Abeta);
   double m = v*((sxy/sig2) + Abeta*betabar);

   return (m + sqrt(v)*gen.normal());
}
