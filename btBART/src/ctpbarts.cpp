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
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "rtnorm.h"
#include "common.h"
//#include "rtruncnorm.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

RcppExport SEXP ctpbarts(
    SEXP _in,            //number of observations in training data
    SEXP _ip1,		//dimension of x
    SEXP _ip2,		//dimension of x
    SEXP _inp,		//number of observations in test data
    SEXP _ix1,		//x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP _ix2,		//x, train,  pxn (transposed so rows are contiguous in memory)
    SEXP _iy,		//y, train,  nx1
    SEXP _ixp1,		//x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP _ixp2,		//x, test, pxnp (transposed so rows are contiguous in memory)
    SEXP _im,		//number of trees
    SEXP _inc1,		//number of cut points
    SEXP _inc2,		//number of cut points
    SEXP _ind,		//number of kept draws (except for thinnning ..)
    SEXP _iburn,		//number of burn-in draws skipped
    SEXP _ipower,
    SEXP _ibase,
    SEXP _binaryOffset,
    SEXP _itau,
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
    SEXP _inkeeptreedraws,
    SEXP _inprintevery, 
    SEXP _X1info,
    SEXP _X2info
)
{
  
  //--------------------------------------------------
  //process args
  size_t n = Rcpp::as<int>(_in);
  size_t p1 = Rcpp::as<int>(_ip1);
  size_t p2 = Rcpp::as<int>(_ip2);
  size_t np = Rcpp::as<int>(_inp);
  
  Rcpp::NumericVector  xv1(_ix1);
  double *ix1 = &xv1[0];
  Rcpp::NumericVector  xv2(_ix2);
  double *ix2 = &xv2[0];
  Rcpp::IntegerVector  yv(_iy); // binary
  int *iy = &yv[0];
  
  Rcpp::NumericVector  xpv1(_ixp1);
  double *ixp1 = &xpv1[0];
  Rcpp::NumericVector  xpv2(_ixp2);
  double *ixp2 = &xpv2[0];
  
  size_t m = Rcpp::as<int>(_im);
  //size_t nc = Rcpp::as<int>(_inc);
  Rcpp::IntegerVector _nc1(_inc1);
  int *numcut1 = &_nc1[0];
  Rcpp::IntegerVector _nc2(_inc2);
  int *numcut2 = &_nc2[0];
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
  Rcpp::NumericMatrix X1info(_X1info);
  Rcpp::NumericMatrix X2info(_X2info);
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
  
  Rcpp::NumericMatrix varprb1(nkeeptreedraws,p1);
  Rcpp::IntegerMatrix varcnt1(nkeeptreedraws,p1);
  Rcpp::NumericMatrix varprb2(nkeeptreedraws,p2);
  Rcpp::IntegerMatrix varcnt2(nkeeptreedraws,p2);
  
  //random number generation
  arn gen;
  
  // BART 1
  bart bm1(m);
  
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
  
  // BART 2
  bart bm2(m);
  
  if(X2info.size()>0) {
    xinfo _x2i;
    _x2i.resize(p2);
    for(size_t i=0;i<p2;i++) {
      _x2i[i].resize(numcut2[i]);
      //Rcpp::IntegerVector cutpts(Xinfo[i]);
      for(size_t j=0;j<numcut2[i];j++) _x2i[i][j]=X2info(i, j);
    }
    bm2.setxinfo(_x2i);
  }
#else
  
#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]
  
  void ctpbarts(
      size_t n,            //number of observations in training data
      size_t p1,		//dimension of x
      size_t p2,		//dimension of x
      size_t np,		//number of observations in test data
      double* ix1,		//x, train,  pxn (transposed so rows are contiguous in memory)
      double* ix2,		//x, train,  pxn (transposed so rows are contiguous in memory)
      int* iy,		//y, train,  nx1
      double* ixp1,		//x, test, pxnp (transposed so rows are contiguous in memory)
      double* ixp2,		//x, test, pxnp (transposed so rows are contiguous in memory)
      size_t m,		//number of trees
      int *numcut1, //size_t nc,		//number of cut points
      int *numcut2, //size_t nc,		//number of cut points
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
    
    std::vector< std::vector<size_t> > varcnt1;
    std::vector< std::vector<double> > varprb1;
    std::vector< std::vector<size_t> > varcnt2;
    std::vector< std::vector<double> > varprb2;
    
    //random number generation
    arn gen(n1, n2);
    
    bart bm1(m);
    bart bm2(m);
#endif
    
    /*
     for(size_t i=0; i<n; ++i) trmean[i]=0.0;
     for(size_t i=0; i<np; ++i) temean[i]=0.0;
     */
    
    
    
    double* iz = new double[n];
    double* iz1 = new double[n];
    double* iz2 = new double[n];
    
    std::stringstream treess1;  //string stream to write trees to
    treess1.precision(10);
    treess1 << nkeeptreedraws << " " << m << " " << p1 << endl;
    
    std::stringstream treess2;  //string stream to write trees to
    treess2.precision(10);
    treess2 << nkeeptreedraws << " " << m << " " << p2 << endl;
    // dart iterations
    std::vector<double> ivarprb1 (p1,0.);
    std::vector<size_t> ivarcnt1 (p1,0);
    std::vector<double> ivarprb2 (p2,0.);
    std::vector<size_t> ivarcnt2 (p2,0);
    
    printf("*****Into main of tpbarts\n");
    
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
    printf("data:n,p1,p2,np: %zu, %zu, %zu\n",n,p1,p2,np);
    printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
    printf("x11,x1[n*p1]: %lf, %lf\n",ix1[0],ix1[n*p1-1]);
    printf("x21,x2[n*p2]: %lf, %lf\n",ix2[0],ix2[n*p2-1]);
    if(np) Rprintf("xp1,xp1[np*p1]: %lf, %lf\n",ixp1[0],ixp1[np*p1-1]);
    if(np) Rprintf("xp2,xp2[np*p2]: %lf, %lf\n",ixp2[0],ixp2[np*p2-1]);
    printf("*****Number of Cut Points for BART 1 : %d ... %d\n", numcut1[0], numcut1[p1-1]);
    printf("*****Number of Cut Points for BART 2: %d ... %d\n", numcut2[0], numcut2[p2-1]);
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
    //init
    for(size_t k=0; k<n; k++){
      if(iy[k]==0) iz[k]= -rtnorm(0., binaryOffset, 1., gen);
      else iz[k]=rtnorm(0., -binaryOffset, 1., gen);
         }
    
    //BART 1
    
    for(size_t k=0; k<n; k++){iz1[k] = iz[k];}
    bm1.setprior(alpha,mybeta,tau);
    bm1.setdata(p1,n,ix1,iz1,numcut1);
   // bm1.setdart(a,b,rho,aug,dart,theta,omega);
    bm1.draw(1.,gen);
    
    // BART 2  
    for(size_t k=0; k<n; k++){iz2[k] = iz[k]-bm1.f(k);}
    
    bm2.setprior(alpha,mybeta,tau);
    bm2.setdata(p2,n,ix2,iz2,numcut2);
  //  bm2.setdart(a,b,rho,aug,dart,theta,omega);
    bm2.draw(1.,gen);
    
    
    //--------------------------------------------------
    
    //--------------------------------------------------
    //temporary storage
    //out of sample fit
    double* fhattest1=0; 
    double* fhattest2=0; 
    if(np) { fhattest1 = new double[np]; 
      fhattest2 = new double[np];}
    
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
    xinfo& x1i = bm1.getxinfo();
    xinfo& x2i = bm2.getxinfo();
    
    for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
       // if(i==(burn/2)&&dart){ bm1.startdart();}
       // if(i==(burn/2)&&dart){ bm2.startdart();}
        
      for(size_t k=0; k<n; k++) {
        if(iy[k]==0) iz[k]= -rtnorm(-(bm1.f(k)+bm2.f(k)), binaryOffset, 1., gen);
        else iz[k]=rtnorm((bm1.f(k)+bm2.f(k)), -binaryOffset, 1., gen);
      }
        //draw bart
        
        for(size_t k=0;k<n;k++){iz1[k]=iz[k]-bm2.f(k);}
      // bm1.setdata(p1,n,ix1,iz1,numcut1);
        bm1.draw(1.,gen);
        //draw bart 2  
        
        for(size_t k=0;k<n;k++){iz2[k]=iz[k]-bm1.f(k);}
       // bm2.setdata(p2,n,ix2,iz2,numcut2);
        bm2.draw(1.,gen);
        
        //bm.draw(rootM, gen);
        
        
        /*
         if(iy[k]==0) iz[k]= -r_lefttruncnorm(-bm.f(k), binaryOffset, 1., gen);
         else iz[k]=r_lefttruncnorm(bm.f(k), -binaryOffset, 1., gen);
         */
      
      
      if(i>=burn) {
        
        //for(size_t k=0;k<n;k++) trmean[k]+=bm.f(k);
        if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
          //index = trcnt*n;;
          //for(size_t k=0;k<n;k++) trdraw[index+k]=bm.f(k);
          for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=(bm1.f(k)+bm2.f(k));
          trcnt+=1;
        }
        keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
        //keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;
        //if(keeptest || keeptestme) bm.predict(p,np,ixp,fhattest);
        if(keeptest) {
          bm1.predict(p1,np,ixp1,fhattest1);
          bm2.predict(p2,np,ixp2,fhattest2);
          //index=tecnt*np;
          //for(size_t k=0;k<np;k++) tedraw[index+k]=fhattest[k];
          for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=(fhattest1[k]+fhattest2[k]);
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
            treess1 << bm1.gettree(j);
            treess2 << bm2.gettree(j);
            /*
#ifndef NoRcpp
             varcount.row(treedrawscnt)=varcount.row(treedrawscnt)+bm.gettree(j).tree2count(p);
             if(treesaslists) lists(j)=bm.gettree(j).tree2list(xi, 0., 1.);
#endif
             */
          }
#ifndef NoRcpp
          //	    if(treesaslists) list_of_lists(treedrawscnt)=lists;
          ivarcnt1=bm1.getnv();
          ivarprb1=bm1.getpv();
          ivarcnt2=bm2.getnv();
          ivarprb2=bm2.getpv();
          size_t k=(i-burn)/skiptreedraws;
          for(size_t j=0;j<p1;j++){
            varcnt1(k,j)=ivarcnt1[j];
            varprb1(k,j)=ivarprb1[j]; 
            //varcnt(i-burn,j)=ivarcnt[j]; 
            //varprb(i-burn,j)=ivarprb[j];
          }
          
          for(size_t j=0;j<p2;j++){ 
            varcnt2(k,j)=ivarcnt2[j];
            //varcnt(i-burn,j)=ivarcnt[j];
            varprb2(k,j)=ivarprb2[j];
            //varprb(i-burn,j)=ivarprb[j];
          }
          
#else
          varcnt1.push_back(bm1.getnv());
          varprb1.push_back(bm1.getpv());
          varcnt2.push_back(bm2.getnv());
          varprb2.push_back(bm2.getpv());
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
    
    if(fhattest1) delete[] fhattest1;
    if(fhattest2) delete[] fhattest2;
    delete[] iz;
    delete[] iz1;
    delete[] iz2;
#ifndef NoRcpp
    //--------------------------------------------------
    //return
    Rcpp::List ret;
    //ret["yhat.train.mean"]=trmean;
    ret["yhat.train"]=trdraw;
    //ret["yhat.test.mean"]=temean;
    ret["yhat.test"]=tedraw;
    //ret["varcount"]=varcount;
    ret["varcount1"]=varcnt1;
    ret["varprob1"]=varprb1;
    ret["varcount2"]=varcnt2;
    ret["varprob2"]=varprb2;
    //for(size_t i=0;i<m;i++) {
    //  bm.gettree(i).pr();
    //}
    
    Rcpp::List x1iret(x1i.size());
    for(size_t i=0;i<x1i.size();i++) {
      Rcpp::NumericVector vtemp1(x1i[i].size());
      std::copy(x1i[i].begin(),x1i[i].end(),vtemp1.begin());
      x1iret[i] = Rcpp::NumericVector(vtemp1);
    }
    
    Rcpp::List x2iret(x2i.size());
    for(size_t i=0;i<x2i.size();i++) {
      Rcpp::NumericVector vtemp2(x2i[i].size());
      std::copy(x2i[i].begin(),x2i[i].end(),vtemp2.begin());
      x2iret[i] = Rcpp::NumericVector(vtemp2);
    }
    
    Rcpp::List treesL1;
    Rcpp::List treesL2;
    //treesL["nkeeptreedraws"] = Rcpp::wrap<int>(nkeeptreedraws); //in trees
    //treesL["ntree"] = Rcpp::wrap<int>(m); //in trees
    //treesL["numx"] = Rcpp::wrap<int>(p); //in cutpoints
    treesL1["cutpoints"] = x1iret;
    treesL1["cutpoints"] = x2iret;
    treesL1["trees"]=Rcpp::CharacterVector(treess1.str());
    treesL2["trees"]=Rcpp::CharacterVector(treess2.str());
    //if(treesaslists) treesL["lists"]=list_of_lists;
    ret["treedraws1"] = treesL1;
    ret["treedraws2"] = treesL2;
    return ret;
#else
    
#endif
    
  }
  
