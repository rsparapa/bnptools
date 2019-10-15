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

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cabart(
   SEXP _type,          //1:wbart, 2:pbart, 3:lbart
   SEXP _in,            //number of observations in training data
   SEXP _ip,            //dimension of x
   SEXP _inp,           //number of observations in test data
   SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,            //y, train,  nx1
   SEXP _idelta,        //censoring indicator
   SEXP _ixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _im,            //number of trees
   SEXP _inc,           //number of cut points
   SEXP _ind,           //number of kept draws (except for thinnning ..)
   SEXP _iburn,         //number of burn-in draws skipped
   SEXP _ithin,         //thinning
   SEXP _ipower,
   SEXP _ibase,
   SEXP _Offset,
   SEXP _itau,
   SEXP _inu,
   SEXP _ilambda,
   SEXP _isigest,
   SEXP _iw,
   SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,            //param a for sparsity prior
   SEXP _ib,            //param b for sparsity prior
   SEXP _irho,          //param rho for sparsity prior (default to p)
   SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP _inprintevery,
   SEXP _Xinfo
)
{
   //process args
   int type = Rcpp::as<int>(_type);
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::NumericVector  yv(_iy); 
   double *iy = &yv[0];
   Rcpp::IntegerVector  deltav(_idelta); 
   int *delta = &deltav[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   size_t thin = Rcpp::as<int>(_ithin);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
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
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
   int *grp = &_grp[0];
   size_t nkeeptrain = nd/thin;     //Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = nd/thin;      //Rcpp::as<int>(_inkeeptest);
   size_t nkeeptreedraws = nd/thin; //Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);

   //random number generation
   arn gen;

   heterbart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cabart(
   int type,            //1:wbart, 2:pbart, 3:lbart
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* iy,		//y, train,  nx1
   int* delta,          //censoring indicator
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int *numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   size_t thin,		//thinning
   double mybeta,
   double alpha,
   double Offset,
   double tau,
   double nu,
   double lambda,
   double sigma,
   double* iw,
   bool dart,           //dart prior: true(1)=yes, false(0)=no   
   double theta,
   double omega, 
   int* grp,
   double a,		//param a for sparsity prior                                          
   double b,		//param b for sparsity prior                                          
   double rho,		//param rho for sparsity prior (default to p)                         
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
   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2);

   heterbart bm(m);
#endif

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   printf("*****Calling abart: type=%d\n", type);

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
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn,nd,thin: %zu,%zu,%zu\n",burn,nd,thin);
// printf("Prior:\nbeta,alpha,tau,nu,lambda,offset: %lf,%lf,%lf,%lf,%lf,%lf\n",
//                    mybeta,alpha,tau,nu,lambda,Offset);
   cout << "*****Prior:beta,alpha,tau,nu,lambda,offset: " 
	<< mybeta << ',' << alpha << ',' << tau << ',' 
        << nu << ',' << lambda << ',' << Offset << endl;
if(type==1) {
   printf("*****sigma: %lf\n",sigma);
   printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
}
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   //printf("*****nkeeptrain,nkeeptest: %zu, %zu\n",nkeeptrain,nkeeptest);
   printf("*****printevery: %zu\n",printevery);

   //--------------------------------------------------
   //create temporaries
   double df=n+nu;
   double *z = new double[n]; 
   double *svec = new double[n]; 
   double *sign;
   if(type!=1) sign = new double[n]; 

   for(size_t i=0; i<n; i++) {
     if(type==1) {
       svec[i] = iw[i]*sigma; 
       z[i]=iy[i]; 
     }
     else {
       svec[i] = 1.;
       if(iy[i]==0) sign[i] = -1.;
       else sign[i] = 1.;
       z[i] = sign[i];
     }
   }
   //--------------------------------------------------
   //set up BART model
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,z,numcut);
   bm.setdart(a,b,rho,aug,dart);

   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);
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
   xinfo& xi = bm.getxinfo();

   for(size_t i=0;i<total;i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      //if(i%printevery==0) printf("%22zu/%zu\r",i,total);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(svec,gen);

      if(type==1) {
      //draw sigma
	double rss=0.;
	for(size_t k=0;k<n;k++) rss += pow((iy[k]-bm.f(k))/(iw[k]), 2.); 
	sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
	sdraw[i]=sigma;
      }

      for(size_t k=0; k<n; k++) {
	if(type==1) {
	  svec[k]=iw[k]*sigma;
	  if(delta[k]==0) z[k]= rtnorm(bm.f(k), iy[k], svec[k], gen);
	}
	else {
	  z[k]= sign[k]*rtnorm(sign[k]*bm.f(k), -sign[k]*Offset, svec[k], gen);
	  if(type==3) 
	    svec[k]=sqrt(draw_lambda_i(pow(svec[k], 2.), sign[k]*bm.f(k), 1000, 1, gen));
	  }
      }
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
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);

   if(fhattest) delete[] fhattest;
   delete[] z;
   delete[] svec;
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

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cgbart(
   SEXP _type,          //1:wbart, 2:pbart, 3:lbart
   SEXP _in,            //number of observations in training data
   SEXP _ip,            //dimension of x
   SEXP _inp,           //number of observations in test data
   SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,            //y, train,  nx1
   SEXP _ixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _im,            //number of trees
   SEXP _inc,           //number of cut points
   SEXP _ind,           //number of kept draws (except for thinnning ..)
   SEXP _iburn,         //number of burn-in draws skipped
   SEXP _ithin,         //thinning
   SEXP _ipower,
   SEXP _ibase,
   SEXP _Offset,
   SEXP _itau,
   SEXP _inu,
   SEXP _ilambda,
   SEXP _isigest,
   SEXP _iw,
   SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,            //param a for sparsity prior
   SEXP _ib,            //param b for sparsity prior
   SEXP _irho,          //param rho for sparsity prior (default to p)
   SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP _inprintevery,
   SEXP _Xinfo
)
{
   //process args
   int type = Rcpp::as<int>(_type);
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::NumericVector  yv(_iy); 
   double *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   size_t thin = Rcpp::as<int>(_ithin);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
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
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
   int *grp = &_grp[0];
   size_t nkeeptrain = nd/thin;     //Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = nd/thin;      //Rcpp::as<int>(_inkeeptest);
   size_t nkeeptreedraws = nd/thin; //Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);

   //random number generation
   arn gen;

   heterbart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cgbart(
   int type,            //1:wbart, 2:pbart, 3:lbart
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int *numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   size_t thin,		//thinning
   double mybeta,
   double alpha,
   double Offset,
   double tau,
   double nu,
   double lambda,
   double sigma,
   double* iw,
   bool dart,           //dart prior: true(1)=yes, false(0)=no   
   double theta,
   double omega, 
   int* grp,
   double a,		//param a for sparsity prior                                          
   double b,		//param b for sparsity prior                                          
   double rho,		//param rho for sparsity prior (default to p)                         
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
   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2);

   heterbart bm(m);
#endif

   /* multiple imputation hot deck implementation
      this will cause trouble for multiple threads
      so it must be done prior to calling C++ */
/*
   bool hotdeck=false;

   std::vector<int> _missing(n*p);
   std::vector<int*> missing(n);
   std::vector<double*> x(n);
   Rcpp::NumericMatrix X(n, p); 

   for(size_t i=0; i<n; ++i) {
     missing[i]=&_missing[i*p];
     x[i]=&ix[i*p];
     for(size_t j=0; j<p; ++j) 
       if(x[i][j]!=x[i][j]) {
	 hotdeck=true;
	 missing[i][j]=1;
       }
       else missing[i][j]=0;
   }     

   if(hotdeck) {
     for(size_t i=0; i<n; ++i) {
       for(size_t j=0; j<p; ++j) { 
	 if(missing[i][j]==1) {
	   while(x[i][j]!=x[i][j]) {
	     size_t k=n*gen.uniform();
	     x[i][j]=x[k][j];
	   }
	 }
	 X(i, j)=x[i][j]; 
       }
     }    
   } 
*/

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   printf("*****Calling gbart: type=%d\n", type);

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
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
//   if(hotdeck) 
//printf("warning: missing elements in x multiply imputed with hot decking\n");
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn,nd,thin: %zu,%zu,%zu\n",burn,nd,thin);
// printf("Prior:\nbeta,alpha,tau,nu,lambda,offset: %lf,%lf,%lf,%lf,%lf,%lf\n",
//                    mybeta,alpha,tau,nu,lambda,Offset);
   cout << "*****Prior:beta,alpha,tau,nu,lambda,offset: " 
	<< mybeta << ',' << alpha << ',' << tau << ',' 
        << nu << ',' << lambda << ',' << Offset << endl;
if(type==1) {
   printf("*****sigma: %lf\n",sigma);
   printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
}
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   //printf("*****nkeeptrain,nkeeptest: %zu, %zu\n",nkeeptrain,nkeeptest);
   printf("*****printevery: %zu\n",printevery);

   //--------------------------------------------------
   //create temporaries
   double df=n+nu;
   double *z = new double[n]; 
   double *svec = new double[n]; 
   double *sign;
   if(type!=1) sign = new double[n]; 

   for(size_t i=0; i<n; i++) {
     if(type==1) {
       svec[i] = iw[i]*sigma; 
       z[i]=iy[i]; 
     }
     else {
       svec[i] = 1.;
       if(iy[i]==0) sign[i] = -1.;
       else sign[i] = 1.;
       z[i] = sign[i];
     }
   }
   //--------------------------------------------------
   //set up BART model
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,z,numcut);
   bm.setdart(a,b,rho,aug,dart);

   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);
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
   bool keeptest,/*keeptestme*/keeptreedraw,type1sigest=(type==1 && lambda!=0.);

   time_t tp;
   int time1 = time(&tp), total=nd+burn;
   xinfo& xi = bm.getxinfo();

   for(size_t i=0;i<total;i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      //if(i%printevery==0) printf("%22zu/%zu\r",i,total);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(svec,gen);

      if(type1sigest) {
      //draw sigma
	double rss=0.;
	for(size_t k=0;k<n;k++) rss += pow((iy[k]-bm.f(k))/(iw[k]), 2.); 
	sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
	sdraw[i]=sigma;
      }

      for(size_t k=0; k<n; k++) {
	if(type==1) svec[k]=iw[k]*sigma;
	else {
	  z[k]= sign[k]*rtnorm(sign[k]*bm.f(k), -sign[k]*Offset, svec[k], gen);
	  if(type==3) 
	    svec[k]=sqrt(draw_lambda_i(pow(svec[k], 2.), sign[k]*bm.f(k), 1000, 1, gen));
	  }
      }

/*
      if(hotdeck) {
	//draw x
	for(size_t h=0; h<n; ++h) {
	  for(size_t j=0; j<p; ++j) {
	    if(missing[h][j]==1) {
	      size_t k=n*gen.uniform();
	      x[h][j]=x[k][j];
	    }
	  }
	}    
      } 
*/

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
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);

   if(fhattest) delete[] fhattest;
   delete[] z;
   delete[] svec;
   if(type!=1) delete[] sign;

#ifndef NoRcpp
   //return list
   Rcpp::List ret;
//   ret["X"]=X; 
   if(type1sigest) ret["sigma"]=sdraw;
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

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)
#define UDRAW(a, b) udraw(a, b)

RcppExport SEXP cgbmm(
		      SEXP _type,    //1:wbart, 2:pbart, 3:lbart
		      SEXP _in,      //number of observations in training data
		      SEXP _ip,      //dimension of x
		      SEXP _inp,     //number of observations in test data
		      SEXP _ix,      //x, train,  pxn 
                      // (transposed so rows are contiguous in memory)
		      SEXP _iy,      //y, train,  nx1
		      SEXP _ixp,     //x, test, pxnp 
                      //(transposed so rows are contiguous in memory)
		      SEXP _iu_train,
		      SEXP _in_train,
		      SEXP _iu,
		      SEXP _iL,
		      SEXP _iB,
		      SEXP _im,      //number of trees
		      SEXP _inc,     //number of cut points
		      SEXP _ind,     //number of kept draws (except thinnning)
		      SEXP _iburn,   //number of burn-in draws skipped
		      SEXP _ithin,   //thinning
		      SEXP _ipower,
		      SEXP _ibase,
		      SEXP _Offset,
		      SEXP _itau,
		      SEXP _inu,
		      SEXP _ilambda,
		      SEXP _isigest,
		      //   SEXP _iw,
		      SEXP _idart,   //dart prior: true(1)=yes, false(0)=no
		      SEXP _itheta,
		      SEXP _iomega,
		      SEXP _igrp,
		      SEXP _ia,      //param a for sparsity prior
		      SEXP _ib,      //param b for sparsity prior
		      SEXP _irho,    //param rho for sparsity prior (default=p)
		      SEXP _iaug,    //category strategy:  true(1)=data augment 
		      //false(0)=degenerate trees
		      SEXP _inprintevery,
		      SEXP _Xinfo
		      )
{
  //process args
  int type = Rcpp::as<int>(_type);
  size_t n = Rcpp::as<int>(_in);
  size_t p = Rcpp::as<int>(_ip);
  size_t np = Rcpp::as<int>(_inp);
  Rcpp::NumericVector  xv(_ix);
  double *ix = &xv[0];
  Rcpp::NumericVector  yv(_iy); 
  double *iy = &yv[0];
  Rcpp::NumericVector  xpv(_ixp);
  double *ixp = &xpv[0];
  Rcpp::IntegerVector _u_train(_iu_train);
  int *u_train = &_u_train[0];
  Rcpp::IntegerVector _n_train(_in_train);
  int *n_train = &_n_train[0];
  Rcpp::NumericVector _u(_iu);
  double *u = &_u[0];
  size_t L = Rcpp::as<int>(_iL);
  double B = Rcpp::as<double>(_iB);
  size_t m = Rcpp::as<int>(_im);
  Rcpp::IntegerVector _nc(_inc);
  int *numcut = &_nc[0];
  size_t nd = Rcpp::as<int>(_ind);
  size_t burn = Rcpp::as<int>(_iburn);
  size_t thin = Rcpp::as<int>(_ithin);
  double mybeta = Rcpp::as<double>(_ipower);
  double alpha = Rcpp::as<double>(_ibase);
  double Offset = Rcpp::as<double>(_Offset);
  double tau = Rcpp::as<double>(_itau);
  double nu = Rcpp::as<double>(_inu);
  double lambda = Rcpp::as<double>(_ilambda);
  double sigma=Rcpp::as<double>(_isigest);
  //   Rcpp::NumericVector  wv(_iw); 
  //   double *iw = &wv[0];
  bool dart;
  if(Rcpp::as<int>(_idart)==1) dart=true;
  else dart=false;
  double a = Rcpp::as<double>(_ia);
  double b = Rcpp::as<double>(_ib);
  double rho = Rcpp::as<double>(_irho);
  bool aug;
  if(Rcpp::as<int>(_iaug)==1) aug=true;
  else aug=false;
  double theta = Rcpp::as<double>(_itheta);
  double omega = Rcpp::as<double>(_iomega);
  Rcpp::IntegerVector _grp(_igrp);
  int *grp = &_grp[0];
  size_t nkeeptrain = nd/thin;     //Rcpp::as<int>(_inkeeptrain);
  size_t nkeeptest = nd/thin;      //Rcpp::as<int>(_inkeeptest);
  size_t nkeeptreedraws = nd/thin; //Rcpp::as<int>(_inkeeptreedraws);
  size_t printevery = Rcpp::as<int>(_inprintevery);
  Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
  Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
  Rcpp::NumericMatrix Xinfo(_Xinfo);
  Rcpp::NumericVector sdraw(nd+burn);
  Rcpp::NumericVector sdudraw(nd+burn);
  Rcpp::NumericMatrix trdraw(nkeeptrain,n);
  Rcpp::NumericMatrix tedraw(nkeeptest,np);
  Rcpp::NumericMatrix udraw(nkeeptrain,L);

  //random number generation
  arn gen;

  //heterbart bm(m);
  bart bm(m);

  if(Xinfo.size()>0) {
    xinfo _xi;
    _xi.resize(p);
    for(size_t i=0;i<p;i++) {
      _xi[i].resize(numcut[i]);
      for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
    }
    bm.setxinfo(_xi);
  }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]
#define UDRAW(a, b) udraw[a][b]

  void cgbmm(
	     int type,    //1:wbart, 2:pbart, 3:lbart
	     size_t n,    //number of observations in training data
	     size_t p,	  //dimension of x
	     size_t np,	  //number of observations in test data
	     double* ix,  //x, train,  pxn 
             //(transposed so rows are contiguous in memory)
	     double* iy,  //y, train,  nx1
	     double* ixp, //x, test, pxnp 
             //(transposed so rows are contiguous in memory)
	     int *u_train,
	     int *n_train,
	     double *u,
	     size_t L,
	     double B,
	     size_t m,	  //number of trees
	     int *numcut, //number of cut points
	     size_t nd,	  //number of kept draws (except for thinnning ..)
	     size_t burn, //number of burn-in draws skipped
	     size_t thin, //thinning
	     double mybeta,
	     double alpha,
	     double Offset,
	     double tau,
	     double nu,
	     double lambda,
	     double sigma,
	     //   double* iw,
	     bool dart,   //dart prior: true(1)=yes, false(0)=no   
	     double theta,
	     double omega, 
	     int* grp,
	     double a,	  //param a for sparsity prior                         
	     double b,	  //param b for sparsity prior                        
	     double rho,  //param rho for sparsity prior (default to p)   
	     bool aug,    //categorical strategy: true(1)=data augment  
	     //false(0)=degenerate trees
	     //   size_t nkeeptrain, //   size_t nkeeptest, //   size_t nkeeptreedraws,
	     size_t printevery,
	     unsigned int n1, // additional parameters needed to call from C++
	     unsigned int n2,
	     double* sdraw,
	     double* sdudraw,
	     double* _trdraw,
	     double* _tedraw,
	     double* _udraw
	     )
  {
  //return data structures (using C++)
  size_t nkeeptrain=nd/thin, nkeeptest=nd/thin, nkeeptreedraws=nd/thin;
  std::vector<double*> trdraw(nkeeptrain);
  std::vector<double*> tedraw(nkeeptest);
  std::vector<double*> udraw(nkeeptrain);

  for(size_t i=0; i<nkeeptrain; ++i) {
    trdraw[i]=&_trdraw[i*n];
    udraw[i]=&_udraw[i*L];
  }
  for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

  //matrix to return dart posteriors (counts and probs)
  std::vector< std::vector<size_t> > varcnt;
  std::vector< std::vector<double> > varprb;

  //random number generation
  arn gen(n1, n2);

  //heterbart bm(m);
  bart bm(m);
#endif

  /* multiple imputation hot deck implementation
     this will cause trouble for multiple threads
     so it must be done prior to calling C++ */
  /*
    bool hotdeck=false;

    std::vector<int> _missing(n*p);
    std::vector<int*> missing(n);
    std::vector<double*> x(n);
    Rcpp::NumericMatrix X(n, p); 

    for(size_t i=0; i<n; ++i) {
    missing[i]=&_missing[i*p];
    x[i]=&ix[i*p];
    for(size_t j=0; j<p; ++j) 
    if(x[i][j]!=x[i][j]) {
    hotdeck=true;
    missing[i][j]=1;
    }
    else missing[i][j]=0;
    }     

    if(hotdeck) {
    for(size_t i=0; i<n; ++i) {
    for(size_t j=0; j<p; ++j) { 
    if(missing[i][j]==1) {
    while(x[i][j]!=x[i][j]) {
    size_t k=n*gen.uniform();
    x[i][j]=x[k][j];
    }
    }
    X(i, j)=x[i][j]; 
    }
    }    
    } 
  */

  std::stringstream treess;  //string stream to write trees to
  treess.precision(10);
  treess << nkeeptreedraws << " " << m << " " << p << endl;

  printf("*****Calling gbmm: type=%d\n", type);

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
  printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
  printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
  printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
  //   if(hotdeck) 
  //printf("warning: missing elements in x multiply imputed with hot decking\n");
  if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
  printf("*****Number of Trees: %zu\n",m);
  printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
  printf("*****burn,nd,thin: %zu,%zu,%zu\n",burn,nd,thin);
  // printf("Prior:\nbeta,alpha,tau,nu,lambda,offset: %lf,%lf,%lf,%lf,%lf,%lf\n",
  //                    mybeta,alpha,tau,nu,lambda,Offset);
  cout << "*****Prior:beta,alpha,tau,nu,lambda,offset: " 
       << mybeta << ',' << alpha << ',' << tau << ',' 
       << nu << ',' << lambda << ',' << Offset << endl;
  if(type==1) {
    printf("*****sigma: %lf\n",sigma);
    //   printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
  }
  cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
       << dart << ',' << theta << ',' << omega << ',' << a << ',' 
       << b << ',' << rho << ',' << aug << endl;
  //printf("*****nkeeptrain,nkeeptest: %zu, %zu\n",nkeeptrain,nkeeptest);
  printf("*****printevery: %zu\n",printevery);

  //create temporaries
  double df=n+nu;
  double *z = new double[n]; 
  //double *svec = new double[n]; 
  double *sign;
  if(type!=1) sign = new double[n]; 
  for(size_t i=0; i<n; i++) {
    if(type==1) {
      //svec[i] = iw[i]*sigma; 
      //svec[i] = sigma; 
      z[i]=iy[i]-Offset; 
    }
    else {
      //svec[i] = 1.;
      if(iy[i]==0) sign[i] = -1.;
      else sign[i] = 1.;
      z[i] = sign[i];
    }
  }

  //double *u = new double[L]; 
  double invB2=pow(B, -2.), invB2_16=16.*invB2, sd_u=B*0.5, 
    tau_u=pow(sd_u, -2.);

  if(u[0]!=u[0]) 
    for(size_t j=0; j<L; j++) {
      u[j]=sd_u*gen.normal();
    }

  //set up BART model
  bm.setprior(alpha,mybeta,tau);
  bm.setdata(p,n,ix,z,numcut);
  bm.setdart(a,b,rho,aug,dart);

  // dart iterations
  std::vector<double> ivarprb (p,0.);
  std::vector<size_t> ivarcnt (p,0);
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
  bool keeptest,keeptreedraw,type1sigest=(type==1 && lambda!=0.);

  time_t tp;
  int time1 = time(&tp), total=nd+burn;
  xinfo& xi = bm.getxinfo();

  for(size_t i=0;i<total;i++) {
    if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
    //if(i%printevery==0) printf("%22zu/%zu\r",i,total);
    if(i==(burn/2)&&dart) bm.startdart();
    //draw bart
    //bm.draw(svec,gen);
    for(size_t k=0;k<n;k++) z[k]=iy[k]-Offset-u[u_train[k]]; 
    bm.draw(sigma,gen);

    if(type1sigest) {
      //draw sigma
      double rss=0.;
      for(size_t k=0;k<n;k++) 
	rss += pow((iy[k]-Offset-bm.f(k)-u[u_train[k]]), 2.); 
      //rss += pow((iy[k]-bm.f(k)-u[u_train[k]])/(iw[k]), 2.); 
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
      sdraw[i]=sigma;
    }

    if(type==2)
      for(size_t k=0; k<n; k++) 
	z[k]= sign[k]*rtnorm(sign[k]*bm.f(k), 
			     -sign[k]*(Offset+u[u_train[k]]), sigma, gen);
    /*
      for(size_t k=0; k<n; k++) {
      if(type==1) svec[k]=sigma;
      //svec[k]=iw[k]*sigma;
      else {
      z[k]= sign[k]*rtnorm(sign[k]*bm.f(k), -sign[k]*(Offset+u[u_train[k]]), sigma, gen);
      if(type==3) 
      svec[k]=sqrt(draw_lambda_i(pow(svec[k], 2.), sign[k]*bm.f(k), 1000, 1, gen));
      }
      }
    */

    // draw tau_u
    double sum_u2;
    sum_u2=0.;
    for(size_t k=0; k<L; k++) sum_u2 += pow(u[k], 2.);
    tau_u=rtgamma(0.5*(L-1.), 0.5*sum_u2, invB2, gen); 
    //if(i<burn) tau_u=std::min(tau_u, invB2_16);
      
    // draw u
    size_t h, n_i;
    double mu_u_i, sd_u_i, prec=pow(sigma, -2.);

    h=0;
    for(size_t k=0; k<L; k++) {
      n_i=n_train[k];
      mu_u_i=0.;
      sd_u_i=pow(tau_u+n_i*prec, -0.5);
      for(size_t j=0; j<n_i; j++) {
	mu_u_i += (iy[h]-Offset-bm.f(h));
	h++;
      }
      mu_u_i *= prec*pow(sd_u_i, 2.);
      u[k]=gen.normal()*sd_u_i+mu_u_i;
    }

    sdudraw[i]=pow(tau_u, -0.5);
    /*
      if(hotdeck) {
      //draw x
      for(size_t h=0; h<n; ++h) {
      for(size_t j=0; j<p; ++j) {
      if(missing[h][j]==1) {
      size_t k=n*gen.uniform();
      x[h][j]=x[k][j];
      }
      }
      }    
      } 
    */

    if(i>=burn) {
      if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
	for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=Offset+bm.f(k);
	for(size_t k=0;k<L;k++) UDRAW(trcnt,k)=u[k];
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
  }
  int time2 = time(&tp);
  printf("time: %ds\n",time2-time1);
  printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);

  if(fhattest) delete[] fhattest;
  delete[] z;
  //delete[] svec;
  if(type!=1) delete[] sign;
  // delete[] u;

#ifndef NoRcpp
  //return list
  Rcpp::List ret;
  //   ret["X"]=X; 
  if(type1sigest) ret["sigma"]=sdraw;
  ret["yhat.train"]=trdraw;
  ret["yhat.test"]=tedraw;
  ret["varcount"]=varcnt;
  ret["varprob"]=varprb;
  ret["u.train"]=udraw;
  ret["sd.u"]=sdudraw;

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

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP clbart(
   SEXP _in,            //number of observations in training data
   SEXP _ip,            //dimension of x
   SEXP _inp,           //number of observations in test data
   SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,            //y, train,  nx1
   SEXP _ixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _im,            //number of trees
   SEXP _inc,           //number of cut points
   SEXP _ind,           //number of kept draws (except for thinnning ..)
   SEXP _iburn,         //number of burn-in draws skipped
   SEXP _ipower,
   SEXP _ibase,
   SEXP _binaryOffset,
   SEXP _itau,
   SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP _ia,            //param a for sparsity prior
   SEXP _ib,            //param b for sparsity prior
   SEXP _irho,          //param rho for sparsity prior (default to p)
   SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
//   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery,
//   SEXP _treesaslists,
   SEXP _Xinfo
)
{

   //--------------------------------------------------
   //process args
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::IntegerVector  yv(_iy); // binary
   int *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   //size_t nc = Rcpp::as<int>(_inc);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   // lbart does not currently employ the binaryOffset trick
   double binaryOffset = Rcpp::as<double>(_binaryOffset);
   double tau = Rcpp::as<double>(_itau);
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
//   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
//   int treesaslists = Rcpp::as<int>(_treesaslists);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
//   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

   //return data structures (using Rcpp)
/*
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector temean(np);
*/
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
//   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);

   //random number generation
   arn gen;

   heterbart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       //Rcpp::IntegerVector cutpts(Xinfo[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void clbart(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   int* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int *numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   double binaryOffset,
   double tau,
   bool dart,           //dart prior: true(1)=yes, false(0)=no                                
   double a,		//param a for sparsity prior                                          
   double b,		//param b for sparsity prior                                          
   double rho,		//param rho for sparsity prior (default to p)                         
   bool aug,		//categorical strategy: true(1)=data augment false(0)=degenerate trees
   size_t nkeeptrain,
   size_t nkeeptest,
//   size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
//   int treesaslists,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
/*
   double* trmean,
   double* temean,
*/
   double* _trdraw,
   double* _tedraw
)
{

   //return data structures (using C++)
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

   //matrix to return dart posteriors (counts and probs)
   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2);

   heterbart bm(m);
#endif

/*
   for(size_t i=0; i<n; i++) trmean[i]=0.0;
   for(size_t i=0; i<np; i++) temean[i]=0.0;
*/

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   printf("*****Into main of lbart\n");

   size_t skiptr,skipte,/*skipteme,*/skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
//   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
//   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %d, %d\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
//   printf("*****Prior:\nbeta,alpha,tau,nu,lambda: %lf,%lf,%lf,%lf,%lf\n",
   printf("*****Prior:\nbeta,alpha,tau: %lf,%lf,%lf\n",
                   mybeta,alpha,tau);
 //  printf("*****sigma: %lf\n",sigma);
   cout << "*****Dirichlet:sparse,a,b,rho,augment: " << dart << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   printf("*****nkeeptrain,nkeeptest: %zu, %zu\n",nkeeptrain,nkeeptest);
   printf("*****printevery: %zu\n",printevery);

   //--------------------------------------------------
   //create logit latents
   //z = f(x) + eps, eps ~N(0,lambda), f ~ BART
   double *z = new double[n]; //latent z's
   double *lambda = new double [n]; //latent lambda's
//   double *yf = new double[n]; //??
   double *svec = new double[n]; //vector of standard dev for bart = sqrt(lambda)
   for(unsigned int i=0; i<n; i++) {
      if(iy[i]>0) z[i] = 1.0;
      else z[i]=-1.0;
      //iy[i]=z[i]; //iy is already +/- 1
      lambda[i] = 1.0;
      svec[i]=1.0; //square root of 1 is 1.
   }
   //newRNGstates();
   //--------------------------------------------------
   //set up BART model
   //heterbart bm(m);
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,z,numcut);
   bm.setdart(a,b,rho,aug,dart);


   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);
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
//   size_t temecnt=0; //count test draws into posterior mean
//   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,/*keeptestme*/keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = bm.getxinfo();

   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %zu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(svec,gen);

      for(size_t k=0; k<n; k++) {
	z[k]= iy[k]*rtnorm(iy[k]*bm.f(k), -iy[k]*binaryOffset, svec[k], gen);
	lambda[k]=draw_lambda_i(lambda[k], iy[k]*bm.f(k), 1000, 1, gen);
//lambda[k]=draw_lambda_i(lambda[k], iy[k]*bm.f(k), 1000, 1, states[0]);
	svec[k] = sqrt(lambda[k]);
      }
/*
      for(unsigned int j=0; j<n; j++) yf[j] = iy[j]*bm.f(j);
      draw_z(n, yf, lambda, z);
      //for(unsigned int j=0; j<n; j++) z[j] *= iy[j];
      draw_lambda(n, yf, 1000, 1, lambda);
      for(unsigned int j=0; j<n; j++) {
	z[j] *= iy[j];
	svec[j] = sqrt(lambda[j]); 
      }
*/
      if(i>=burn) {
         //for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            //index = trcnt*n;;
            //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
 //        keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
//         if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
            //index=tecnt*np;
            //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
            tecnt+=1;
         }
         // if(keeptestme) {
         //    for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
         //    temecnt+=1;
         // }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
	   // #ifndef NoRcpp
	   // Rcpp::List lists(m*treesaslists);
	   // #endif

            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);

	      #ifndef NoRcpp
	      //varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
	      //if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
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
//	    #ifndef NoRcpp
	    //if(treesaslists) list_of_lists(treedrawscnt)=lists;
//	    #endif
//            treedrawscnt +=1;
         }
      }
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
//   for(size_t k=0;k<n;k++) trmean[k]/=nd;
//   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
   printf("check counts\n");
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);
   //printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   //--------------------------------------------------

   if(fhattest) delete[] fhattest;
   delete[] z;
//   delete[] yf;
   delete[] lambda;
   delete[] svec;
   //deleteRNGstates();

#ifndef NoRcpp
   //--------------------------------------------------
   //return
   Rcpp::List ret;
//   ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
//   ret["yhat.test.mean"]=temean;
   ret["yhat.test"]=tedraw;
//   ret["varcount"]=varcount;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;

   //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
   //}

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
   //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
   //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
//   if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   return ret;
#else

#endif

}

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cpbart(
   SEXP _in,            //number of observations in training data
   SEXP _ip,		//dimension of x
   SEXP _inp,		//number of observations in test data
   SEXP _ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,		//y, train,  nx1
   SEXP _ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _im,		//number of trees
   SEXP _inc,		//number of cut points
   SEXP _ind,		//number of kept draws (except for thinnning ..)
   SEXP _iburn,		//number of burn-in draws skipped
   SEXP _ipower,
   SEXP _ibase,
   SEXP _binaryOffset,
   SEXP _itau,
//   SEXP _iM, // number of shards for Modified LISA MCMC
   SEXP _idart,
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,
   SEXP _ib,
   SEXP _irho,
   SEXP _iaug,
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
//   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery,
//   SEXP _treesaslists,
   SEXP _Xinfo
)
{

   //--------------------------------------------------
   //process args
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::IntegerVector  yv(_iy); // binary
   int *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   //size_t nc = Rcpp::as<int>(_inc);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   double binaryOffset = Rcpp::as<double>(_binaryOffset);
   double tau = Rcpp::as<double>(_itau);
//   double rootM = sqrt(Rcpp::as<double>(_iM));
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
   int *grp = &_grp[0];
   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
//   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
//   int treesaslists = Rcpp::as<int>(_treesaslists);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
   //Rcpp::List Xinfo(_Xinfo);
//   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

   //return data structures (using Rcpp)
/*
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector temean(np);
*/
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
//   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);

   //random number generation
   arn gen;

   bart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       //Rcpp::IntegerVector cutpts(Xinfo[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cpbart(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   int* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int *numcut, //size_t nc,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   double binaryOffset,
   double tau,
   bool dart,
   double theta,
   double omega,
   int *grp,
   double a,
   double b,
   double rho,
   bool aug,
   size_t nkeeptrain,
   size_t nkeeptest,
   //size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
//   int treesaslists,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
/*
   double* trmean,
   double* temean,
*/
   double* _trdraw,
   double* _tedraw
)
{

   //return data structures (using C++)
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2);

   bart bm(m);
#endif

/*
   for(size_t i=0; i<n; ++i) trmean[i]=0.0;
   for(size_t i=0; i<np; ++i) temean[i]=0.0;
*/

   double* iz = new double[n];

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;
   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   printf("*****Into main of pbart\n");

   size_t skiptr,skipte,skiptreedraws;
   //size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
/*
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
*/
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %d, %d\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
   printf("*****Prior:mybeta,alpha,tau: %lf,%lf,%lf\n",
                   mybeta,alpha,tau);
   printf("*****binaryOffset: %lf\n",binaryOffset);
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   printf("*****nkeeptrain,nkeeptest,nkeeptreedraws: %zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptreedraws);
/*
   printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
*/
   printf("*****printevery: %zu\n",printevery);
   //printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);
   printf("*****skiptr,skipte,skiptreedraws: %zu,%zu,%zu\n",skiptr,skipte,skiptreedraws);

   //--------------------------------------------------
   //bart bm(m);
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,iz,numcut);
   bm.setdart(a,b,rho,aug,dart,theta,omega);
   //--------------------------------------------------
//init
  for(size_t k=0; k<n; k++) {
    if(iy[k]==0) iz[k]= -rtnorm(0., binaryOffset, 1., gen);
    else iz[k]=rtnorm(0., -binaryOffset, 1., gen);
/*
    if(iy[k]==0) iz[k]= -r_lefttruncnorm(0., binaryOffset, 1., gen);
    else iz[k]=r_lefttruncnorm(0., -binaryOffset, 1., gen);
*/
  }

   //--------------------------------------------------

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
//   size_t temecnt=0; //count test draws into posterior mean
//   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,/*keeptestme,*/keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = bm.getxinfo();

   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(1., gen);
      //bm.draw(rootM, gen);

  for(size_t k=0; k<n; k++) {
    if(iy[k]==0) iz[k]= -rtnorm(-bm.f(k), binaryOffset, 1., gen);
    else iz[k]=rtnorm(bm.f(k), -binaryOffset, 1., gen);
/*
    if(iy[k]==0) iz[k]= -r_lefttruncnorm(-bm.f(k), binaryOffset, 1., gen);
    else iz[k]=r_lefttruncnorm(bm.f(k), -binaryOffset, 1., gen);
*/
  }

      if(i>=burn) {

         //for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            //index = trcnt*n;;
            //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         //keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
         //if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
            //index=tecnt*np;
            //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
            tecnt+=1;
         }
/*
         if(keeptestme) {
            for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
            temecnt+=1;
         }
*/
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
//	   #ifndef NoRcpp
//	   Rcpp::List lists(m*treesaslists);
//	   #endif

            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);
/*
	      #ifndef NoRcpp
	      varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
	      if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
	      #endif
*/
	    }
	    #ifndef NoRcpp
//	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      //varcnt(i-burn,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	      //varprb(i-burn,j)=ivarprb[j];
	    }
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif

            //treedrawscnt +=1;
         }
      }
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
/*
   for(size_t k=0;k<n;k++) trmean[k]/=nd;
   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
*/
   printf("check counts\n");
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);
   //printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   //--------------------------------------------------
   //PutRNGstate();

   if(fhattest) delete[] fhattest;
   delete[] iz;

#ifndef NoRcpp
   //--------------------------------------------------
   //return
   Rcpp::List ret;
   //ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   //ret["yhat.test.mean"]=temean;
   ret["yhat.test"]=tedraw;
   //ret["varcount"]=varcount;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;

   //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
   //}

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
   //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
   //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   //if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   return ret;
#else

#endif

}

typedef std::vector<tree> vtree;

#ifdef _OPENMP
void local_getpred(size_t nd, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat);
//void local_getpred(size_t nd, size_t p, size_t m, size_t np, xinfo xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat);
#endif

void getpred(int beg, int end, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat);

RcppExport SEXP cpwbart(
   SEXP _itrees,		//treedraws list from fbart
   SEXP _ix,			//x matrix to predict at
   SEXP _itc			//thread count
)
{
   Rprintf("*****In main of C++ for bart prediction\n");

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
   Rcpp::NumericMatrix xpred(_ix);
   size_t np = xpred.ncol();
   cout << "from x,np,p: " << xpred.nrow() << ", " << xpred.ncol() << endl;
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
   double *px = &xpred(0,0);

   #ifndef _OPENMP
   cout << "***using serial code\n";
   getpred(0, nd-1, p, m, np,  xi,  tmat, px,  yhat);
   #else
   if(tc==1) {cout << "***using serial code\n"; getpred(0, nd-1, p, m, np,  xi,  tmat, px,  yhat);}
   else {
      cout << "***using parallel code\n";
#pragma omp parallel num_threads(tc)
      local_getpred(nd,p,m,np,xi,tmat,px,yhat);
   }
   #endif

   //--------------------------------------------------
   Rcpp::List ret;
   ret["yhat.test"] = yhat;
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
//void local_getpred(size_t nd, size_t p, size_t m, size_t np, xinfo xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat)
void local_getpred(size_t nd, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat)
{


   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int h = nd/thread_count; int beg = my_rank*h; int end = beg+h-1;
   
   getpred(beg,end,p,m,np,xi,tmat,px,yhat);
}
#endif

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cwbart(
   SEXP _in,            //number of observations in training data
   SEXP _ip,		//dimension of x
   SEXP _inp,		//number of observations in test data
   SEXP _ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,		//y, train,  nx1
   SEXP _ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _im,		//number of trees
   SEXP _inc,		//number of cut points
   SEXP _ind,		//number of kept draws (except for thinnning ..)
   SEXP _iburn,		//number of burn-in draws skipped
   SEXP _ipower,
   SEXP _ibase,
   SEXP _itau,
   SEXP _inu,
   SEXP _ilambda,
   SEXP _isigest,
   SEXP _iw,
   SEXP _idart,
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,
   SEXP _ib,
   SEXP _irho,
   SEXP _iaug,
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery,
//   SEXP _treesaslists,
   SEXP _Xinfo
)
{

   //--------------------------------------------------
   //process args
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::NumericVector  yv(_iy); 
   double *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   size_t m = Rcpp::as<int>(_im);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   //size_t nc = Rcpp::as<int>(_inc);
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
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
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
   int *grp = &_grp[0];
   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
//   int treesaslists = Rcpp::as<int>(_treesaslists);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
//   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

   //return data structures (using Rcpp)
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector temean(np);
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
//   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);

   //random number generation
   arn gen;

   heterbart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       //Rcpp::IntegerVector cutpts(Xinfo[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cwbart(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int* numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   double tau,
   double nu,
   double lambda,
   double sigma,
   double* iw,
   bool dart,
   double theta,
   double omega,
   int *grp,
   double a,
   double b,
   double rho,
   bool aug,
   size_t nkeeptrain,
   size_t nkeeptest,
   size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
//   int treesaslists,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
   double* trmean,
   double* temean,
   double* sdraw,
   double* _trdraw,
   double* _tedraw
)
{

   //return data structures (using C++)
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2); 

   heterbart bm(m);
#endif

   for(size_t i=0;i<n;i++) trmean[i]=0.0;
   for(size_t i=0;i<np;i++) temean[i]=0.0;

   printf("*****Into main of wbart\n");
   //-----------------------------------------------------------

   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
   printf("*****Prior:beta,alpha,tau,nu,lambda: %lf,%lf,%lf,%lf,%lf\n",
                   mybeta,alpha,tau,nu,lambda);
   printf("*****sigma: %lf\n",sigma);
   printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
   printf("*****printevery: %zu\n",printevery);
   printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);

   //--------------------------------------------------
   //heterbart bm(m);
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,iy,numcut);
   bm.setdart(a,b,rho,aug,dart,theta,omega);

   //--------------------------------------------------
   //sigma
   //gen.set_df(n+nu);
   double *svec = new double[n];
   for(size_t i=0;i<n;i++) svec[i]=iw[i]*sigma;

   //--------------------------------------------------

   std::stringstream treess;  //string stream to write trees to  
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;
   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; //posterior mean for prediction
   if(np) { fhattest = new double[np]; }
   double restemp=0.0,rss=0.0;


   //--------------------------------------------------
   //mcmc
   printf("\nMCMC\n");
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   size_t temecnt=0; //count test draws into posterior mean
   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,keeptestme,keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = bm.getxinfo();

   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(svec,gen);
      //draw sigma
      rss=0.0;
      for(size_t k=0;k<n;k++) {restemp=(iy[k]-bm.f(k))/(iw[k]); rss += restemp*restemp;}
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(n+nu));
      for(size_t k=0;k<n;k++) svec[k]=iw[k]*sigma;
      sdraw[i]=sigma;
      if(i>=burn) {
         for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            //index = trcnt*n;;
            //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k);
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
         if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
         if(keeptest) {
            //index=tecnt*np;
            //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
            tecnt+=1;
         }
         if(keeptestme) {
            for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
            temecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
//	   #ifndef NoRcpp
//	   Rcpp::List lists(m*treesaslists);
//	   #endif

            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);
/*      
	      #ifndef NoRcpp
	      varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
	      if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
	      #endif
*/
	    }
            #ifndef NoRcpp
//	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      //varcnt(i-burn,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	      //varprb(i-burn,j)=ivarprb[j];
	    }
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif

            treedrawscnt +=1;
         }
      }
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   for(size_t k=0;k<n;k++) trmean[k]/=nd;
   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
   printf("check counts\n");
   printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   //--------------------------------------------------
   //PutRNGstate();

   if(fhattest) delete[] fhattest;
   if(svec) delete [] svec;

   //--------------------------------------------------
   //return
#ifndef NoRcpp
   Rcpp::List ret;
   ret["sigma"]=sdraw;
   ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   ret["yhat.test.mean"]=temean;
   ret["yhat.test"]=tedraw;
   //ret["varcount"]=varcount;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;

   //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
   //}

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
   //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
   //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
//   if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   return ret;
#else

#endif

}

#ifndef NoRcpp

RcppExport SEXP mc_cores_openmp() {

#else

int mc_cores_openmp() {

#endif

#ifdef _OPENMP

int mc_cores_openmp=omp_get_num_threads();

#else

int mc_cores_openmp=0;

#endif

#ifndef NoRcpp

return Rcpp::wrap(mc_cores_openmp);

#else

return mc_cores_openmp;

#endif

}

