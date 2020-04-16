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

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cqbart(
   //SEXP _type,          //1:wbart, 2:pbart, 3:lbart
   SEXP _in,            //number of observations in training data
   SEXP _ip1,            //dimension of x1
   SEXP _ip2,            //dimension of x2
   SEXP _inp,           //number of observations in test data
   SEXP _ix1,            //x1, train for cure status,  pxn (transposed so rows are contiguous in memory)
   SEXP _ix2,           //x2, train for y
   SEXP _iy,            //y, train,  nx1
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
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);

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

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cqbart(
   //int type,            //1:wbart, 2:pbart, 3:lbart
   size_t n,            //number of observations in training data
   size_t p1,		//dimension of x1
   size_t p2,           //dimension of x2
   size_t np,		//number of observations in test data
   double* ix1,		//x1, train for cure status,  pxn (transposed so rows are contiguous in memory)
   double* ix2,         //x2, train for y
   double* iy,		//y, train,  nx1
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
   double* _trdraw,
   double* _tedraw
)
{
   //return data structures (using C++)
   size_t nkeeptrain=nd/thin, nkeeptest=nd/thin, nkeeptreedraws=nd/thin;
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

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
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   printf("*****Calling qbart: \n");

   size_t skiptr=thin, skipte=thin, skiptreedraws=thin;
/*
   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;
*/

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p1,p2,np: %zu, %zu, %zu, %zu\n",n,p1,p2,np);
   printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   printf("x11,x1[n*p1]: %lf, %lf\n",ix1[0],ix1[n*p1-1]);
   printf("x21,x2[n*p2]: %lf, %lf\n",ix2[0],ix2[n*p2-1]);
   if(np1) printf("xp11,xp1[np*p1]: %lf, %lf\n",ixp1[0],ixp1[np*p1-1]);
   if(np2) printf("xp21,xp2[np*p2]: %lf, %lf\n",ixp2[0],ixp2[np*p2-1]);
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
   double df=n+nu;
   double *z1 = new double[n];
   double *z2 = new double[n];
   int *q = new int[n];
   double *pz1 = new double[n];
   double *pt = new double[n];
   //double *svec = new double[n]; 
   
   //double *sign;
   //if(type!=1) sign = new double[n]; 

   for(size_t k=0; k<n; k++) {
     //if(delta[k]==0) z1[k] = -rtnorm(0., binaryOffset, 1., gen);
     z1[k] = rtnorm(0., -binaryOffset, 1., gen);
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
   std::vector<double> ivarprb1 (p,0.);
   std::vector<size_t> ivarcnt1 (p,0);
   std::vector<double> ivarprb2 (p,0.);
   std::vector<size_t> ivarcnt2 (p,0);
   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; 
   if(np) { fhattest = new double[np]; }

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

      //draw sigma
      double rss=0.;
      for(size_t k=0;k<n;k++) rss += pow((iy[k]-bm2.f(k))/(iw[k]), 2.); 
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
      sdraw[i]=sigma;

      //
      for (size_t k=0; k<n; k++){
	if (q[k]==0) z1[k] = -rtnorm(-bm1.f(k), binaryOffset, 1., gen);
	else {
	  z1[k] = rtnorm(bm1.f(k), -binaryOffset, 1., gen);
	}
	#ifndef NoRcpp
	pz1[k] = R::pnorm(z1[k]+binaryOffset,0,1,1,0);
	#else
	pz1[k] = ::pnorm(z1[k]+binaryOffset,0,1,1,0);
	#endif
	if (delta[k]==0) pt[k] = pt[k]*
      }

      //for(size_t k=0; k<n; k++) {
	//if(type==1) {
	  //svec[k]=iw[k]*sigma;
	  if(delta[k]==0) z[k]= rtnorm(bm.f(k), iy[k], sigma, gen);
	//}
	//else {
	  //z[k]= sign[k]*rtnorm(sign[k]*bm.f(k), -sign[k]*Offset, sigma, gen);
	  //if(type==3) 
	    //svec[k]=sqrt(draw_lambda_i(pow(svec[k], 2.), sign[k]*bm.f(k), 1000, 1, gen));
	  //}
      //}
      if(i>=burn) {
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=Offset+bm.f(k);
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=Offset+fhattest[k];
            tecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);

	      #ifndef NoRcpp
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	    }
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif
	    }
         }
      }

      bm1.draw(1.,gen);
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);

   if(fhattest) delete[] fhattest;
   delete[] z;
   //delete[] svec;
   if(type!=1) delete[] sign;

#ifndef NoRcpp
   //return list
   Rcpp::List ret;
   if(type==1) ret["sigma"]=sdraw;
   ret["yhat.train"]=trdraw;
   ret["yhat.test"]=tedraw;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   ret["treedraws"] = treesL;
   
   return ret;
#else

#endif

}
