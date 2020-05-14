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

#include <BART3.h>
#include "qbart.h"


#ifndef NoRcpp

#define TRDRAW1(a, b) trdraw1(a, b)
#define TEDRAW1(a, b) tedraw1(a, b)
#define QDRAW(a, b) qdraw(a, b)
#define STDRAW(a, b) stdraw(a, b)
#define PTDRAW(a, b) ptdraw(a, b)
#define TRDRAW2(a, b) trdraw2(a, b)
#define TEDRAW2(a, b) tedraw2(a, b)

RcppExport SEXP cqbart(
   //SEXP _type,          //1:wbart, 2:pbart, 3:lbart
   SEXP _in,            //number of observations in training data
   SEXP _ip1,            //dimension of x1
   SEXP _ip2,            //dimension of x2
   SEXP _inp,           //number of observations in test data
   SEXP _ix1,            //x1, train for cure status,  pxn (transposed so rows are contiguous in memory)
   SEXP _ix2,           //x2, train for y
   SEXP _iy,            //y, train,  nx1
   SEXP _imaxy,
   SEXP _idelta,        //censoring indicator
   SEXP _iq0,             //initial cure status
   SEXP _ixp1,           //x, test1, pxnp (transposed so rows are contiguous in memory)
   SEXP _ixp2,           //x, test2 for y
   SEXP _im,            //number of trees
   SEXP _inc1,           //number of cut points for x1
   SEXP _inc2,           //number of cut points for x2
   SEXP _ind,           //number of kept draws (except for thinnning ..)
   SEXP _iburn,         //number of burn-in draws skipped
   SEXP _ithin,         //thinning
   SEXP _ipower,
   SEXP _ibase,
   SEXP _binaryOffset,  //center of uncured rate
   SEXP _Offset,        //center of log(time)
   SEXP _itau,
   SEXP _inu,
   SEXP _ilambda,
   SEXP _isigest,
   SEXP _iw,
   SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP _itheta,
   SEXP _iomega,
   //SEXP _igrp,
   SEXP _ia,            //param a for sparsity prior
   SEXP _ib,            //param b for sparsity prior
   SEXP _irho1,          //param rho1 for sparsity prior (default to p1)
   SEXP _irho2,          //param rho2 for sparsity prior (default to p2)
   SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP _inprintevery,
   SEXP _X1info,
   SEXP _X2info
)
{
   //process args
   //int type = Rcpp::as<int>(_type);
   size_t n = Rcpp::as<int>(_in);
   size_t p1 = Rcpp::as<int>(_ip1);
   size_t p2 = Rcpp::as<int>(_ip2);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv1(_ix1);
   double *ix1 = &xv1[0];
   Rcpp::NumericVector  xv2(_ix2);
   double *ix2 = &xv2[0];
   Rcpp::NumericVector  yv(_iy); 
   double *iy = &yv[0];
   double maxy = Rcpp::as<double>(_imaxy);
   Rcpp::IntegerVector  deltav(_idelta); 
   int *delta = &deltav[0];
   Rcpp::IntegerVector  qv(_iq0); 
   int *q0 = &qv[0];
   Rcpp::NumericVector  xpv1(_ixp1);
   double *ixp1 = &xpv1[0];
   Rcpp::NumericVector  xpv2(_ixp2);
   double *ixp2 = &xpv2[0];
   size_t m = Rcpp::as<int>(_im);
   Rcpp::IntegerVector _nc1(_inc1);
   int *numcut1 = &_nc1[0];
   Rcpp::IntegerVector _nc2(_inc2);
   int *numcut2 = &_nc2[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   size_t thin = Rcpp::as<int>(_ithin);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   double binaryOffset = Rcpp::as<double>(_binaryOffset);
   double Offset = Rcpp::as<double>(_Offset);
   double tau = Rcpp::as<double>(_itau);
   double nu = Rcpp::as<double>(_inu);
   double lambda = Rcpp::as<double>(_ilambda);
   double sigma=Rcpp::as<double>(_isigest);
   Rcpp::NumericVector  wv(_iw); 
   double *iw = &wv[0];
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho1 = Rcpp::as<double>(_irho1);
   double rho2 = Rcpp::as<double>(_irho2);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   //Rcpp::IntegerVector _grp(_igrp);
   //int *grp = &_grp[0];
   size_t nkeeptrain = nd/thin;     //Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = nd/thin;      //Rcpp::as<int>(_inkeeptest);
   size_t nkeeptreedraws = nd/thin; //Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
   Rcpp::NumericMatrix varprb1(nkeeptreedraws,p1);
   Rcpp::IntegerMatrix varcnt1(nkeeptreedraws,p1);
   Rcpp::NumericMatrix varprb2(nkeeptreedraws,p2);
   Rcpp::IntegerMatrix varcnt2(nkeeptreedraws,p2);
   Rcpp::NumericMatrix X1info(_X1info);
   Rcpp::NumericMatrix X2info(_X2info);
   Rcpp::NumericVector sdraw(nd/thin);
   Rcpp::NumericMatrix trdraw1(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw1(nkeeptest,np);
   Rcpp::NumericMatrix qdraw(nkeeptrain,n);
   Rcpp::NumericMatrix stdraw(nkeeptrain,n);
   Rcpp::NumericMatrix ptdraw(nkeeptrain,n);
   Rcpp::NumericMatrix trdraw2(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw2(nkeeptest,np);

   //random number generation
   arn gen;

   //probit BART for cure status
   bart bm1(m);

   if(X1info.size()>0) {
     xinfo _x1i;
     _x1i.resize(p1);
     for(size_t i=0;i<p1;i++) {
       _x1i[i].resize(numcut1[i]);
       for(size_t j=0;j<numcut1[i];j++) _x1i[i][j]=X1info(i, j);
     }
     bm1.setxinfo(_x1i);
   }

   //log-normal BART for y
   qbart bm2(m);

   if(X2info.size()>0) {
     xinfo _x2i;
     _x2i.resize(p2);
     for(size_t i=0;i<p2;i++) {
       _x2i[i].resize(numcut2[i]);
       for(size_t j=0;j<numcut2[i];j++) _x2i[i][j]=X2info(i, j);
     }
     bm2.setxinfo(_x2i);
   }
#else

#define TRDRAW1(a, b) trdraw1[a][b]
#define TEDRAW1(a, b) tedraw1[a][b]
#define QDRAW(a, b) qdraw[a][b]
#define STDRAW(a, b) stdraw[a][b]
#define PTDRAW(a, b) ptdraw[a][b]
#define TRDRAW2(a, b) trdraw2[a][b]
#define TEDRAW2(a, b) tedraw2[a][b]

void cqbart(
   //int type,            //1:wbart, 2:pbart, 3:lbart
   size_t n,            //number of observations in training data
   size_t p1,		//dimension of x1
   size_t p2,           //dimension of x2
   size_t np,		//number of observations in test data
   double* ix1,		//x1, train for cure status,  pxn (transposed so rows are contiguous in memory)
   double* ix2,         //x2, train for y
   double* iy,		//y, train,  nx1
   double maxy,
   int* delta,          //censoring indicator
   int* q0,              //initial cure status
   double* ixp1,		//x, test1, pxnp (transposed so rows are contiguous in memory)
   double* ixp2,           //x, test2 for y
   size_t m,		//number of trees
   int *numcut1,		//number of cut points for x1
   int* numcut2,           //number of cut points for x2
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   size_t thin,		//thinning
   double mybeta,
   double alpha,
   double binaryOffset,
   double Offset,
   double tau,
   double nu,
   double lambda,
   double sigma,
   double* iw,
   bool dart,           //dart prior: true(1)=yes, false(0)=no   
   double theta,
   double omega, 
   //int* grp,
   double a,		//param a for sparsity prior   
   double b,		//param b for sparsity prior          
   double rho1,		//param rho1 for sparsity prior (default to p1)  
   double rho2,		//param rho2 for sparsity prior (default to p2)
   bool aug,		//categorical strategy: true(1)=data augment false(0)=degenerate trees
//   size_t nkeeptrain, //   size_t nkeeptest, //   size_t nkeeptreedraws,
   size_t printevery,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
   double* sdraw,
   double* _trdraw1,
   double* _tedraw1,
   double* _qdraw,
   double* _stdraw,
   double* _ptdraw,
   double* _trdraw2,
   double* _tedraw2
)
{
   //return data structures (using C++)
   size_t nkeeptrain=nd/thin, nkeeptest=nd/thin, nkeeptreedraws=nd/thin;
   std::vector<double*> trdraw1(nkeeptrain);
   std::vector<double*> tedraw1(nkeeptest);
   std::vector<double*> qdraw(nkeeptrain);
   std::vector<double*> stdraw(nkeeptrain);
   std::vector<double*> ptdraw(nkeeptrain);
   std::vector<double*> trdraw2(nkeeptrain);
   std::vector<double*> tedraw2(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) {trdraw1[i]=&_trdraw1[i*n]; trdraw2[i]=&_trdraw2[i*n];qdraw[i]=&_qdraw[i*n];stdraw[i]=&_stdraw[i*n];ptdraw[i]=&_ptdraw[i*n];}
   for(size_t i=0; i<nkeeptest; ++i) {tedraw1[i]=&_tedraw1[i*np]; tedraw2[i]=&_tedraw2[i*np];}

   //matrix to return dart posteriors (counts and probs)
   std::vector< std::vector<size_t> > varcnt1;
   std::vector< std::vector<double> > varprb1;
   std::vector< std::vector<size_t> > varcnt2;
   std::vector< std::vector<double> > varprb2;

   //random number generation
   arn gen(n1, n2);

   bart bm1(m);
   qbart bm2(m);
#endif

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p1 << endl;

   printf("*****Calling qbart: \n");

   size_t skiptr=thin, skipte=thin, skiptreedraws=thin;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p1,p2,np: %zu, %zu, %zu, %zu\n",n,p1,p2,np);
   printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   printf("x11,x1[n*p1]: %lf, %lf\n",ix1[0],ix1[n*p1-1]);
   printf("x21,x2[n*p2]: %lf, %lf\n",ix2[0],ix2[n*p2-1]);
   if(np) printf("xp11,xp1[np*p1]: %lf, %lf\n",ixp1[0],ixp1[np*p1-1]);
   if(np) printf("xp21,xp2[np*p2]: %lf, %lf\n",ixp2[0],ixp2[np*p2-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points for BART1: %d ... %d\n", numcut1[0], numcut1[p1-1]);
   printf("*****Number of Cut Points for BART2: %d ... %d\n", numcut2[0], numcut2[p2-1]);
   printf("*****burn,nd,thin: %zu,%zu,%zu\n",burn,nd,thin);
// printf("Prior:\nbeta,alpha,tau,nu,lambda,offset: %lf,%lf,%lf,%lf,%lf,%lf\n",
//                    mybeta,alpha,tau,nu,lambda,Offset);
   cout << "*****Prior:beta,alpha,tau,nu,lambda,offset: " 
	<< mybeta << ',' << alpha << ',' << tau << ',' 
        << nu << ',' << lambda << ',' << Offset << endl;
//if(type==1) {
//   printf("*****sigma: %lf\n",sigma);
//   printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
//}
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho1 << ',' << rho2 << ',' << aug << endl;
   //printf("*****nkeeptrain,nkeeptest: %zu, %zu\n",nkeeptrain,nkeeptest);
   printf("*****printevery: %zu\n",printevery);

   //--------------------------------------------------
   //create temporaries
   double *z1 = new double[n];
   double *z2 = new double[n];
   int *q = new int[n];
   double *pz1 = new double[n];
   double *pt = new double[n];
   double *st = new double[n];
   //double *svec = new double[n]; 
   
   //double *sign;
   //if(type!=1) sign = new double[n]; 
   
   for(size_t k=0; k<n; k++) {
     if(q0[k]==0) z1[k] = -rtnorm(0., binaryOffset, 1., gen);
     else z1[k] = rtnorm(0., -binaryOffset, 1., gen);
     z2[k] = iy[k]; 
     q[k] = q0[k];
   }
   //--------------------------------------------------
   //set up BART1 model
   bm1.setprior(alpha,mybeta,tau);
   bm1.setdata(p1,n,ix1,z1,numcut1);
   bm1.setdart(a,b,rho1,aug,dart);

   
   //set up BART2 model
   bm2.setprior(alpha,mybeta,tau);
   bm2.setdata(p2,n,ix2,z2,q,numcut2);
   bm2.setdart(a,b,rho2,aug,dart);

   // dart iterations
   std::vector<double> ivarprb1 (p1,0.);
   std::vector<size_t> ivarcnt1 (p1,0);
   std::vector<double> ivarprb2 (p2,0.);
   std::vector<size_t> ivarcnt2 (p2,0);
   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest1=0; 
   double* fhattest2=0;
   double* fhattr2= new double[n];
   if(np) { fhattest1 = new double[np]; fhattest2 = new double[np]; }

   //--------------------------------------------------
   //mcmc
   printf("\nMCMC\n");
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   bool keeptest,/*keeptestme*/keeptreedraw;

   time_t tp;
   int time1 = time(&tp), total=nd+burn;
   xinfo& xi1 = bm1.getxinfo();
   xinfo& xi2 = bm2.getxinfo();

   for(size_t i=0;i<total;i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,total);
      //if(i%printevery==0) printf("%22zu/%zu\r",i,total);
      if(i==(burn/2)&&dart) {bm1.startdart(); bm2.startdart();}
      
      //draw bart2
      bm2.draw(sigma, gen);

      //predict
      bm2.predict(p2,n,ix2,fhattr2);

      //draw sigma
      double rss=0.;
      double df=nu;
      for(size_t k=0;k<n;k++) {
	if (q[k]==1) {
	  rss += pow((iy[k]-fhattr2[k])/(iw[k]), 2.); 
	  df += 1;
	}
      }
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
      //sdraw[i]=sigma;

	
      for (size_t k=0; k<n; k++){
	if (delta[k] == 0){
	  #ifndef NoRcpp
	  pz1[k] = R::pnorm(bm1.f(k)+binaryOffset,0,1,1,0);
	  st[k] = R::pnorm(iy[k], fhattr2[k], sigma, 0, 0);
	  #else
	  pz1[k] = ::pnorm(bm1.f(k)+binaryOffset,0,1,1,0);
	  st[k] = ::pnorm(iy[k], fhattr2[k], sigma, 0, 0);
	  #endif
	  pt[k] = pz1[k]*st[k]/(1-pz1[k]+pz1[k]*st[k]);
	  //draw q
	  #ifndef NoRcpp
	  q[k] = R::rbinom(1, pt[k]);
	  #else
	  q[k] = ::rbinom(1, pt[k]);
	  #endif
	}
	else {q[k] = 1; st[k] = 1; pt[k] = 1;}
	//draw z1
	if (q[k]==0) z1[k] = -rtnorm(-bm1.f(k), binaryOffset, 1., gen);
	else  {
	  z1[k] = rtnorm(bm1.f(k), -binaryOffset, 1., gen);
	  if (delta[k] == 0) {
	    //draw z2
	    z2[k] = rtnorm(fhattr2[k], iy[k], sigma, gen);
	    z2[k] = std::min(maxy, z2[k]);  //maximum cap
	  }
	}
      }
      //draw bart1
      bm1.draw(1., gen);

//printf("binaryOffset: %lf\n", binaryOffset);
      
      if(i>=burn) {
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            for(size_t k=0;k<n;k++) {
	      TRDRAW1(trcnt,k)=binaryOffset+bm1.f(k);
	      qdraw(trcnt,k)=q[k];
	      stdraw(trcnt,k)=Offset+bm2.f(k);
	      ptdraw(trcnt,k)=pt[k];
	      TRDRAW2(trcnt,k)=Offset+fhattr2[k];
	    }
	    sdraw[trcnt]=sigma;
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         if(keeptest) {
	   bm1.predict(p1,np,ixp1,fhattest1);
	   bm2.predict(p2,np,ixp2,fhattest2);
	   for(size_t k=0;k<np;k++) {
	     TEDRAW1(tecnt,k)=binaryOffset+fhattest1[k];
	     TEDRAW2(tecnt,k)=Offset+fhattest2[k];
	   }
	   tecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
	   for(size_t j=0;j<m;j++) {
	     treess << bm1.gettree(j);
             #ifndef NoRcpp
	     ivarcnt1=bm1.getnv();
	     ivarprb1=bm1.getpv();
	     ivarcnt2=bm2.getnv();
	     ivarprb2=bm2.getpv();
	     size_t k=(i-burn)/skiptreedraws;
	     for(size_t j=0;j<p1;j++){
	      varcnt1(k,j)=ivarcnt1[j];
	      varprb1(k,j)=ivarprb1[j];
	     }
	     for(size_t j=0;j<p2;j++){
	      varcnt2(k,j)=ivarcnt2[j];
	      varprb2(k,j)=ivarprb2[j];
	     }
             #else
	     varcnt1.push_back(bm1.getnv());
	     varprb1.push_back(bm1.getpv());
	     varcnt2.push_back(bm2.getnv());
	     varprb2.push_back(bm2.getpv());
             #endif
	   }
         }
      }
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);

   if(fhattest1) delete[] fhattest1; if(fhattest2) delete[] fhattest2;
   delete[] z1; delete[] z2; delete[] pz1; delete[] pt; delete[] st; delete[] fhattr2;
   //delete[] svec;
   

#ifndef NoRcpp
   //return list
   Rcpp::List ret;
   ret["sigma"]=sdraw;
   ret["y1hat.train"]=trdraw1;
   ret["y1hat.test"]=tedraw1;
   ret["varcount1"]=varcnt1;
   ret["varprob1"]=varprb1;
   ret["y2hat.train"]=trdraw2;
   ret["y2hat.test"]=tedraw2;
   ret["varcount2"]=varcnt2;
   ret["varprob2"]=varprb2;
   ret["qdraws"]=qdraw;
   ret["bm2fdraws"]=stdraw;
   ret["ptdraws"]=ptdraw;

   Rcpp::List xi1ret(xi1.size());
   for(size_t i=0;i<xi1.size();i++) {
      Rcpp::NumericVector vtemp(xi1[i].size());
      std::copy(xi1[i].begin(),xi1[i].end(),vtemp.begin());
      xi1ret[i] = Rcpp::NumericVector(vtemp);
   }
   Rcpp::List xi2ret(xi2.size());
   for(size_t i=0;i<xi2.size();i++) {
      Rcpp::NumericVector vtemp(xi2[i].size());
      std::copy(xi2[i].begin(),xi2[i].end(),vtemp.begin());
      xi2ret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   treesL["cutpoints1"] = xi1ret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   ret["treedraws"] = treesL;
   
   return ret;
#else

#endif

}
