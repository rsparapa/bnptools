/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
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

#include <ctime>
#include "common.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"
//#include "latent.h"
//#include "rand_draws.h"
#include "rtnorm.h"
#include "lambda.h"
#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP cltrtbart(
    SEXP _in,            //number of observations in training data
    SEXP _ip1,		//dimension of x1 
    SEXP _inp,		//number of observations in test data
    SEXP _ix1,		//x, train1,  pxn (transposed so rows are contiguous in memory)
    SEXP _itrt,	// Treatment train
    SEXP _iy,		 //y, train,  nx1
    SEXP _ixp1,		//x, test1,  
    SEXP _itrtp,		//Treatment,test 
    SEXP _imub0,   // Prior mean for beta
    SEXP _is2b0,   // Prior variance for beta
    SEXP _im,		 //number of trees
    SEXP _inc1,		//number of cut points for x1 
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
    SEXP _inkeeptreedraws,
    SEXP _inprintevery,
    SEXP _X1info 
)
{
  
  //--------------------------------------------------
  //process args
  size_t n = Rcpp::as<int>(_in); 
  size_t p1 = Rcpp::as<int>(_ip1); 
  size_t np = Rcpp::as<int>(_inp);
  
  Rcpp::NumericVector  xv1(_ix1);
  double *ix1 = &xv1[0];
  Rcpp::NumericVector  trtv(_itrt);
  double *itrt = &trtv[0]; 
  Rcpp::NumericVector  yv(_iy); 
  double *iy = &yv[0]; 
  Rcpp::NumericVector  xpv1(_ixp1);
  double *ixp1 = &xpv1[0]; 
  Rcpp::NumericVector  trtpv(_itrtp);
  double *itrtp = &trtpv[0]; 
  size_t m = Rcpp::as<int>(_im);
  
  Rcpp::IntegerVector _nc1(_inc1);
  int *numcut1 = &_nc1[0]; 
  
  double mu0 = Rcpp::as<double>(_imub0);
  double sigmab20 = Rcpp::as<double>(_is2b0);
  
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
  
  
  Rcpp::NumericMatrix varprb1(nkeeptreedraws,p1);
  Rcpp::IntegerMatrix varcnt1(nkeeptreedraws,p1); 
  Rcpp::NumericMatrix X1info(_X1info); 
  
  
  //   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);
  
  //return data structures (using Rcpp)
  /*
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector temean(np);
   */
  Rcpp::NumericVector betadraw(nd+burn);
  Rcpp::NumericMatrix trdraw(nkeeptrain,n);
  Rcpp::NumericMatrix tedraw(nkeeptest,np);
  //   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
  
  //random number generation
  arn gen;
  
  // BART 
  heterbart bm1(m);
  
  if(X1info.size()>0) {
    xinfo _x1i;
    _x1i.resize(p1);
    for(size_t i=0;i<p1;i++) {
      _x1i[i].resize(numcut1[i]);
      //Rcpp::IntegerVector cutpts(Xinfo[i]);
      for(size_t j=0;j<numcut1[i];j++) _x1i[i][j]=X1info(i, j);
    }
    bm1.setxinfo(_x1i);
  }
  
  
#else
  
#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]
  
  void cltrtbart(
      size_t n,            //number of observations in training data
      size_t p1,		//dimension of x1 
      size_t np,		//number of observations in test data
      double* ix1,		//x, train1,  p1xn (transposed so rows are contiguous in memory)
      double* itrt,		// Treatment,train  
      int* iy,		//y, train,  nx1
      double* ixp1,		//x, test1, p1xnp (transposed so rows are contiguous in memory)
      double* itrtp,	 // Treatment,train  
      double* mu0,
      double* sigmab20,
      size_t m,		//number of trees
      int *numcut1,		//number of cut points for x1 
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
    std::vector< std::vector<size_t> > varcnt1;
    std::vector< std::vector<double> > varprb1; 
    //random number generation
    arn gen(n1, n2);
    
    heterbart bm1(m); 
#endif
    
    /*
     for(size_t i=0; i<n; i++) trmean[i]=0.0;
     for(size_t i=0; i<np; i++) temean[i]=0.0;
     */
    
    
    std::stringstream treess1;  //string stream to write trees to  
    treess1.precision(10);
    treess1<< nkeeptreedraws << " " << m << " " << p1 << endl;
    
    printf("*****Into main of ltrtbart\n");
    
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
    Rprintf("*****Data:\n");
    printf("data:n,p1,np: %zu,%zu,%zu\n",n,p1,np);
    printf("y1,yn: %d, %d\n",iy[0],iy[n-1]); 
    Rprintf("x11,x[n*p1]: %lf, %lf\n",ix1[0],ix1[n*p1-1]);
    printf("trt1,trtn: %lf, %lf\n",itrt[0],itrt[n-1]);
    if(np) Rprintf("xp1,xp1[np*p1]: %lf, %lf\n",ixp1[0],ixp1[np*p1-1]);
    if(np) Rprintf("trtp,trtp[n]: %lf, %lf\n",itrtp[0],itrtp[np-1]);
    printf("*****Number of Trees: %zu\n",m);
    printf("*****Number of Cut Points for BART  : %d ... %d\n", numcut1[0], numcut1[p1-1]); 
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
   // double *yf = new double[n]; //??
    double *svec = new double[n]; //vector of standard dev for bart = sqrt(lambda)
    
    for(unsigned int i=0; i<n; i++) {
      if(iy[i]>0) z[i] = 1.0;
      else z[i]=-1.0;
      //iy[i]=z[i]; //iy is already +/- 1
      lambda[i] = 1.0;
      svec[i]=1.0; //square root of 1 is 1.
    }
    
    
    // Initialize BART1
    //Fit a parametric model to trt and draw beta
    // Get posterior mean and standard deviation for beta based on Holmes and Held
    
    double muhatb=0.0,sighatb=0.0;
    double  AWZ=0.0,AWA=0.0;
    double beta=0.0;
    
    //Create a diagonal matrix
    
    Rcpp::NumericMatrix W(n,n); 
    for(size_t i=0;i<n;i++) {
      for(size_t j=0;j<n;j++){
        W(i,j)=0;
      if(i==j){W(i,j)= 1/lambda[i];}
      }
      }
    
 
    for(size_t k=0;k<n;k++){
      // AWZ += (itrt[k]*W(k,k)*z[k]);
      // AWA += (itrt[k]*W(k,k)*itrt[k]);
      // 
      AWZ += (itrt[k]*(1/lambda[k])*z[k]);
      AWA += (itrt[k]*(1/lambda[k])*itrt[k]);
    }
    
    double *z1 = new double[n]; //latent z's  
    double V = 0.0;  
    
    V = 1/((1/sigmab20) + AWA);
    sighatb = sqrt(V);
    muhatb = V*((1/sigmab20)*mu0 + AWZ);
    
    beta =  sighatb*gen.normal() + muhatb ;
    
    for(unsigned int j=0; j<n; j++) {z1[j]=z[j]-beta*itrt[j];}
    
    //  newRNGstates();
    //--------------------------------------------------
    //set up BART1 model
    //heterbart bm(m);
    bm1.setprior(alpha,mybeta,tau);
    bm1.setdata(p1,n,ix1,z1,numcut1);
   // bm1.setdart(a,b,rho,aug,dart);
    bm1.draw(svec,gen);
    
    
    // dart iterations
    std::vector<double> ivarprb1 (p1,0.);
    std::vector<size_t> ivarcnt1 (p1,0); 
    //--------------------------------------------------
    //temporary storage
    //out of sample fit 
    double* fhattest1=0; //posterior mean for prediction  
    if(np){ 
      fhattest1 = new double[np]; 
    }
    
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
    xinfo& x1i  = bm1.getxinfo(); 
    
    for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %zu)\n",i,nd+burn);
      //if(i==(burn/2)&&dart) bm1.startdart(); 
      //draw bart
      //Get latent z
      
      for(size_t k=0; k<n; k++) {
        z[k]= iy[k]*rtnorm(iy[k]*(beta*itrt[k] + bm1.f(k)), -iy[k]*binaryOffset, svec[k], gen);
        lambda[k]=draw_lambda_i(lambda[k], iy[k]*(beta*itrt[k] + bm1.f(k)), 1000, 1, gen);
        //lambda[k]=draw_lambda_i(lambda[k], iy[k]*bm.f(k), 1000, 1, states[0]);
        svec[k] = sqrt(lambda[k]);
      }
      
      //Frida
      // Draw samples for parametric model
      
      
     // double muhatb=0.0,sighatb=0.0;
      double  AWZ=0.0,AWA=0.0;
    // double beta=0.0;
      
      //Create a diagonal matrix
      
      for(size_t i=0;i<n;i++) {
        for(size_t j=0;j<n;j++){
          W(i,j)=0;
        if(i==j){W(i,j)= 1/lambda[i];}
      }
      }
      
      for(size_t k=0;k<n;k++){
        // AWZ += (itrt[k]*W(k,k)*(z[k]-bm1.f(k)));
        // AWA += (itrt[k]*W(k,k)*itrt[k]);
        // 
        AWZ += (itrt[k]*(1/lambda[k])*(z[k]-bm1.f(k)));
        AWA += (itrt[k]*(1/lambda[k])*itrt[k]);
        
        
      }
      
     // double V = 0.0;  
      
      V = pow(pow(sigmab20,-1) + AWA,-1);
      sighatb = sqrt(V);
      muhatb = V*(pow(sigmab20,-1)*mu0 + AWZ);
      
      beta =  sighatb*gen.normal() + muhatb ;
       
      betadraw[i] = beta;
      //       for(unsigned int j=0; j<n; j++) yf[j] = iy[j]*(bm1.f(j)+bm2.f(j));
      //       draw_z(n, yf, lambda, z);
      //       //for(unsigned int j=0; j<n; j++) z[j] *= iy[j];
      //       draw_lambda(n, yf, 1000, 1, lambda);
      //       for(unsigned int j=0; j<n; j++) {
      // 	       z[j] *= iy[j];
      //         	svec[j] = sqrt(lambda[j]); // this line was missing
      //       }
      
      //draw bart 1 
      for(unsigned int k=0;k<n;k++){z1[k]= z[k]-beta*itrt[k];}
      //bm1.setdata(p1,n,ix1,z1,numcut1);
      bm1.draw(svec,gen);
      
      
      if(i>=burn) {
        //for(size_t k=0;k<n;k++) trmean[k]+=(bm1.f(k)+bm2.f(k));
        if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
          //index = trcnt*n;;
          // for(size_t k=0;k<n;k++) trdraw(trcnt,k)=(beta*itrt[k]+bm1.f(k));
          for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=(beta*itrt[k]+bm1.f(k));
          trcnt+=1;
        }
        keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
        //        keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
        //         if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
        
        
        if(keeptest) {
          bm1.predict(p1,np,ixp1,fhattest1); 
          //index=tecnt*np;
          // for(size_t k=0;k<np;k++) tedraw(tecnt,k)=(fhattest1[k] + fhattest2[k]);
          for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=(beta*itrtp[k]+fhattest1[k]);
          tecnt+=1;
        }
        // if(keeptestme) {
        //    for(size_t k=0;k<np;k++) temean[k]+=(fhattest1[k] + fhattest2[k]);
        //    temecnt+=1;
        // }
        keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
        if(keeptreedraw) {
          // #ifndef NoRcpp
          // Rcpp::List lists(m*treesaslists);
          // #endif
          
          for(size_t j=0;j<m;j++) {
            treess1 << bm1.gettree(j); 
#ifndef NoRcpp
            //varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
            //if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
            ivarcnt1=bm1.getnv();
            ivarprb1=bm1.getpv(); 
            size_t k=(i-burn)/skiptreedraws;
            for(size_t j=0;j<p1;j++){
              varcnt1(k,j)=ivarcnt1[j];
              varprb1(k,j)=ivarprb1[j]; 
            }
            
            
#else
            varcnt1.push_back(bm1.getnv());
            varprb1.push_back(bm1.getpv()); 
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
    //PutRNGstate();
    
    if(fhattest1) delete[] fhattest1; 
    delete[] z;
    delete[] z1; 
   // delete[] yf;
    delete[] lambda;
    delete[] svec;
    // deleteRNGstates();
    
#ifndef NoRcpp
    //--------------------------------------------------
    //return
    Rcpp::List ret;
    //   ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    ret["beta"]=betadraw;
    //   ret["yhat.test.mean"]=temean;
    ret["yhat.test"]=tedraw;
    //   ret["varcount"]=varcount;
    ret["varcount1"]=varcnt1;
    ret["varprob1"]=varprb1; 
    //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
    //}
    
    Rcpp::List x1iret(x1i.size());
    for(size_t i=0;i<x1i.size();i++) {
      Rcpp::NumericVector vtemp(x1i[i].size());
      std::copy(x1i[i].begin(),x1i[i].end(),vtemp.begin());
      x1iret[i] = Rcpp::NumericVector(vtemp);
    }
    
    Rcpp::List treesL1; 
    //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
    //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
    //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
    treesL1["cutpoints"] = x1iret; 
    treesL1["trees"]=Rcpp::CharacterVector(treess1.str()); 
    //   if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws1"] = treesL1; 
    
    return ret;
#else
    
#endif
    
  }
  
