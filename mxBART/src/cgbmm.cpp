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

/*
#include <ctime>
#include "common.h"
#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "mixedbart.h"
//#include "heterbart.h"
#include "rtnorm.h"
#include "rtgamma.h"
#include "lambda.h"
*/
#include <BART3.h>
#include "mixedbart.h"
#include "RcppEigen.h"

typedef std::vector<double> v1d;
typedef std::vector<v1d> v2d;
typedef std::vector<v2d> v3d;

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)

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
		      SEXP _iid_train, //matrix of cluster / group ids
		      SEXP _iz_train,
		      SEXP _iz_cols,
		      SEXP _in_train,  //matrix of cluster / group sizes (no of replicates)
		      SEXP _iid_train_no,
		      SEXP _iL,
		      SEXP _ipriorNo, //prior of sigma^2_u
		      SEXP _ipriorDf,    //df value for prior parameter
		      SEXP _ipriorScale,
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
		      SEXP _itheta,  //value of sparsity parameter
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
  Rcpp::IntegerVector _z_cols(_iz_cols);
  int *z_cols = &_z_cols[0];
  size_t sum_z_cols=0;
  size_t id_train_no = Rcpp::as<int>(_iid_train_no);
  for(size_t i=0;i<id_train_no;i++) sum_z_cols+=z_cols[i];
  Eigen::Map<Eigen::MatrixXd> z_trainall = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(_iz_train);
  Eigen::Map<Eigen::MatrixXd> n_trainall = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(_in_train);
  Eigen::Map<Eigen::MatrixXd> id_trainall = Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(_iid_train);
  Rcpp::IntegerVector _L(_iL);
  int *L = &_L[0];
  Rcpp::IntegerVector _priorNo(_ipriorNo);
  int *priorNo = &_priorNo[0];
  Rcpp::NumericVector _priorDf(_ipriorDf);
  double *priorDf = &_priorDf[0];
  Rcpp::NumericVector _priorScalev(_ipriorScale);
  double *priorScalep = &_priorScalev[0];
  v2d priorScale;
  priorScale.resize(id_train_no);
  size_t tmp_it=0;
  for(size_t lev=0;lev<id_train_no;lev++){
    priorScale[lev].resize(z_cols[lev]);
    for(size_t re=0;re<z_cols[lev];re++){
      priorScale[lev][re]=priorScalep[tmp_it];
      tmp_it++;
    }
  }
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
  bool useDart;
  if(Rcpp::as<int>(_idart)==1) useDart=true;
  else useDart=false;
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
  std::vector<std::vector<Eigen::MatrixXd> > redraws;
  std::vector<std::vector<Eigen::MatrixXd> > varcovdraws;
  redraws.resize(id_train_no);
  varcovdraws.resize(id_train_no);
  for(size_t lev=0;lev<id_train_no;lev++) {
    redraws[lev].resize(nkeeptrain);
    varcovdraws[lev].resize(nkeeptrain);
    for(size_t i=0;i<nkeeptrain;i++) {
      redraws[lev][i].resize(L[lev],z_cols[lev]);
      varcovdraws[lev][i].resize(z_cols[lev],z_cols[lev]);
    }
  }
  Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
  Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
  Rcpp::NumericMatrix Xinfo(_Xinfo);
  Rcpp::NumericVector sdraw(nkeeptrain);
  Rcpp::NumericVector thetadraws(nkeeptrain);
  Rcpp::NumericMatrix trdraw(nkeeptrain,n);
  Rcpp::NumericMatrix tedraw(nkeeptest,np);

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
	     int *id_trainv,
	     double *z_train,
	     int *z_cols,
	     int id_train_no,
	     int *n_train,
	     int *L,
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
	     bool useDart,   //dart prior: true(1)=yes, false(0)=no   
	     double theta,
	     double omega, 
	     int* priorNo,
	     double* priorDf,
	     double* priorScale,
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
	     double* _trdraw,
	     double* _tedraw,
	     double* _udraw
	     )
  {
  //return data structures (using C++)
  size_t nkeeptrain=nd/thin, nkeeptest=nd/thin, nkeeptreedraws=nd/thin;
  std::vector<double*> trdraw(nkeeptrain);
  std::vector<double*> tedraw(nkeeptest);
  

  for(size_t i=0; i<nkeeptrain; ++i) {
    trdraw[i]=&_trdraw[i*n];
  }
  for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

  //matrix to return dart posteriors (counts and probs)
  v2d varcnt;
  v2d varprb;

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

  // Random effect stuff, initial values are draw from prior
  size_t cur_col=0;
  size_t cur_clust=0;
  std::vector<Eigen::MatrixXd> current_randomEffects(id_train_no);
  std::vector<Eigen::MatrixXd> current_varianceComponents(id_train_no);
  std::vector<Eigen::MatrixXd> z_train(id_train_no);
  std::vector<Eigen::VectorXd> id_train(id_train_no);
  std::vector<Eigen::VectorXd> n_train(id_train_no);
  std::vector<std::vector<Eigen::MatrixXd> > RE_zzt(id_train_no);
  std::vector<std::vector<Eigen::VectorXd> > RE_z(id_train_no);
  Eigen::VectorXd RE_sumztre(n);
  RE_sumztre.setZero();
  // First loop over levels
  cout << "*****Mixed BART Settings\n";
  cout << "*****Number of levels: " << id_train_no+1 << '\n';
  for(size_t lev=0;lev<id_train_no;lev++){
    cout << "**Level: " << lev+2 << '\n';
    cout << "Units: " << L[lev] << '\n';
    cout << "Random Effects (per unit): " << z_cols[lev] << '\n';
    cout << "prior,df " << priorNo[lev] << " " << priorDf[lev] << '\n';
    cout << "scale: ";
    for(size_t k=0;k<z_cols[lev];k++) cout << priorScale[lev][k] << ' ';
    cout << '\n';
    current_randomEffects[lev].resize(L[lev],z_cols[lev]);   
    // Set initial random effect values to 0
    current_randomEffects[lev].setZero(); 
    // Set initial variance component values to 1
    current_varianceComponents[lev].resize(z_cols[lev],z_cols[lev]);
    current_varianceComponents[lev].setZero();
    for(size_t j=0;j<z_cols[lev];j++) current_varianceComponents[lev](j,j)=1;
    z_train[lev]=z_trainall.block(0,cur_col,n,z_cols[lev]);
    id_train[lev]=id_trainall.col(lev);
    n_train[lev]=n_trainall.col(lev);
    RE_zzt[lev].resize(n);
    RE_z[lev].resize(n);
    // Now loop over observations
    for(size_t obs=0;obs<n;obs++){
      RE_zzt[lev][obs]=z_train[lev].row(obs).transpose()*z_train[lev].row(obs);
      RE_z[lev][obs]=z_train[lev].row(obs).transpose();
    }
  }
  // Draws performed now. Again loop over levels and then over clusters, drawing the random effects and then the variance components.
  // These are initial draws from the prior.
  for(size_t lev=0;lev<id_train_no;lev++){
    Eigen::MatrixXd sum_rereT;
    sum_rereT.resize(z_cols[lev],z_cols[lev]);
    sum_rereT.setZero();
    // Variance component draw
    Eigen::MatrixXd RE_sumzzt(z_cols[lev],z_cols[lev]);;
    Eigen::VectorXd RE_sumzdiff(z_cols[lev]);
    Eigen::MatrixXd RE_varcov(z_cols[lev],z_cols[lev]);
    Eigen::VectorXd RE_mean(z_cols[lev]);
    Eigen::VectorXd RE_tmpre(z_cols[lev]);
    double RE_sumother=0.;
    RE_sumzzt.setZero();
    RE_sumzdiff.setZero();
    for(size_t obs=0;obs<n;obs++){
      if(cur_clust!=id_train[lev][obs]){
	get_reParams(RE_varcov, RE_mean, current_varianceComponents[lev], RE_sumzzt, RE_sumzdiff, sigma);
	RE_tmpre.setZero();
	draw_randomEffects(RE_tmpre, RE_mean, RE_varcov, gen);
	current_randomEffects[lev].row(cur_clust)=RE_tmpre.transpose();
	RE_sumzzt.setZero();
	RE_sumzdiff.setZero();
	cur_clust=id_train[lev][obs];
	sum_rereT+=current_randomEffects[lev].row(id_train[lev][obs]).transpose()*current_randomEffects[lev].row(id_train[lev][id_train[lev][id_train[lev][obs]]]);
      }
      RE_sumzzt+=RE_zzt[lev][obs];
      for(size_t lev2=0;lev2<id_train_no;lev2++){
	if(lev2!=lev) RE_sumother+=RE_z[lev2][obs].transpose()*current_randomEffects[lev2].row(id_train[lev2][obs]).transpose();
      }
      RE_sumzdiff+=RE_z[lev][obs]*(iy[obs]-0-RE_sumother);
      RE_sumother=0;
    }
    get_reParams(RE_varcov, RE_mean, current_varianceComponents[lev], RE_sumzzt, RE_sumzdiff, sigma);
    RE_tmpre.setZero();
    draw_randomEffects(RE_tmpre, RE_mean, RE_varcov, gen);
    current_randomEffects[lev].row(cur_clust)=RE_tmpre.transpose();
    RE_sumzzt.setZero();
    RE_sumzdiff.setZero();
    sum_rereT+=current_randomEffects[lev].row(cur_clust).transpose()*current_randomEffects[lev].row(cur_clust);
    draw_varianceComponents(current_varianceComponents[lev], sum_rereT, priorScale[lev], priorDf[lev], L[lev], priorNo[lev], gen);
    for(size_t obs=0;obs<n;obs++) RE_sumztre(obs)+=z_train[lev].row(obs)*current_randomEffects[lev].row(id_train[lev][obs]).transpose();
    cur_clust=0;
  }

  // for(size_t i=0;i<sd_u_cols;i++){
  //   sigU_prior(sd_u_cols,sigUpr,sigUscale,sigUdf,sd_u[i],prec_u[i],gen);
  // }
  // for(size_t i=0;i<L;i++){
  //   draw_u_prior(u[i],sd_u,gen);
  // }
  
  
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
       << useDart << ',' << theta << ',' << omega << ',' << a << ',' 
       << b << ',' << rho << ',' << aug << endl;
  //printf("*****nkeeptrain,nkeeptest: %zu, %zu\n",nkeeptrain,nkeeptest);
  printf("*****printevery: %zu\n",printevery);

  //create temporaries
  double df=n+nu;
  double *z = new double[n];
  //double *svec = new double[n]; 
  for(size_t i=0; i<n; i++) {
    if(type==1) {
      //svec[i] = iw[i]*sigma; 
      //svec[i] = sigma; 
       z[i]=iy[i]-RE_sumztre(i);
    }
    else {
      //svec[i] = 1.;
      if(iy[i]==0) {
	z[i]=-rtnorm(-RE_sumztre(i), -Offset, 1., gen);
      }
      else{
	z[i]=rtnorm(RE_sumztre(i), Offset, 1., gen);
      }
    }
  }

  //set up BART model
  bm.setprior(alpha,mybeta,tau);
  bm.setdata(p,n,ix,z,numcut);
  bm.setdart(a,b,rho,aug,useDart);

  // dartth iterations
  v1d ivarprb (p,0.);
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
    if(useDart&(i==(size_t)burn/2)) {
      bm.startdart();
    }
    //draw bart
    if(type==1) bm.draw(sigma, gen);
    else bm.draw(1., gen);
    //cout << bm.gettheta() << '\n';;
    //bm.draw(svec,gen);
    if(type==1){
      for(size_t k=0;k<n;k++) {
	z[k]=iy[k]-RE_sumztre(k);
      } 
    }
    else if(type==2) {
      for(size_t k=0; k<n; k++) {
	if(iy[k]==1) z[k]=rtnorm(bm.f(k)+RE_sumztre(k), -Offset, 1.,gen);
	else z[k]=-rtnorm(-bm.f(k)-RE_sumztre(k), Offset, 1., gen);
      }
    }
    if(type1sigest) {
      //draw sigma
      double rss=0.;
      for(size_t k=0;k<n;k++)
	rss += pow(iy[k]-bm.f(k)-RE_sumztre(k), 2.); 
      //rss += pow((iy[k]-bm.f(k)-u[u_train[k]])/(iw[k]), 2.); 
      sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
      }
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

  // Draws performed now. Again loop over levels and then over clusters, drawing the random effects and then the variance components.
  // These are initial draws from the prior.
  RE_sumztre.setZero();
  Eigen::MatrixXd sum_rereT;
  for(size_t lev=0;lev<id_train_no;lev++){
    sum_rereT.resize(z_cols[lev],z_cols[lev]);
    sum_rereT.setZero();
    Eigen::MatrixXd RE_sumzzt;
    Eigen::VectorXd RE_sumzdiff;
    double RE_sumother=0.;
    RE_sumzzt.resize(z_cols[lev],z_cols[lev]);
    RE_sumzdiff.resize(z_cols[lev]);
    RE_sumzzt.setZero();
    RE_sumzdiff.setZero();
    for(size_t obs=0;obs<n;obs++){
      if(cur_clust!=id_train[lev][obs]){
	Eigen::MatrixXd RE_varcov(z_cols[lev],z_cols[lev]);
	Eigen::VectorXd RE_mean(z_cols[lev]);
	Eigen::VectorXd RE_tmpre(z_cols[lev]);
	get_reParams(RE_varcov,RE_mean,current_varianceComponents[lev],RE_sumzzt,RE_sumzdiff,sigma);
	RE_tmpre.setZero();
	draw_randomEffects(RE_tmpre, RE_mean, RE_varcov, gen);
	current_randomEffects[lev].row(cur_clust)=RE_tmpre.transpose();
	RE_sumzzt.setZero();
	RE_sumzdiff.setZero();
	cur_clust=id_train[lev][obs];
	sum_rereT+=current_randomEffects[lev].row(id_train[lev][obs]).transpose()*current_randomEffects[lev].row(id_train[lev][obs]);
      }
      RE_sumzzt+=RE_zzt[lev][obs];
      for(size_t lev2=0;lev2<id_train_no;lev2++){
	if(lev2!=lev) RE_sumother+=RE_z[lev2][obs].transpose()*current_randomEffects[lev2].row(id_train[lev2][obs]).transpose();
      }
      RE_sumzdiff+=RE_z[lev][obs]*(iy[obs]-bm.f(obs)-RE_sumother);
      RE_sumother=0.;
    }
    // so final cluster
    Eigen::MatrixXd RE_varcov(z_cols[lev],z_cols[lev]);
    Eigen::VectorXd RE_mean(z_cols[lev]);
    Eigen::VectorXd RE_tmpre(z_cols[lev]);
    get_reParams(RE_varcov, RE_mean, current_varianceComponents[lev], RE_sumzzt, RE_sumzdiff, sigma);
    RE_tmpre.setZero();
    draw_randomEffects(RE_tmpre, RE_mean, RE_varcov, gen);
    current_randomEffects[lev].row(cur_clust)=RE_tmpre.transpose();
    RE_sumzzt.setZero();
    RE_sumzdiff.setZero();
    cur_clust=0;
    sum_rereT+=current_randomEffects[lev].row(cur_clust).transpose()*current_randomEffects[lev].row(cur_clust);
    draw_varianceComponents(current_varianceComponents[lev], sum_rereT, priorScale[lev], priorDf[lev], L[lev], priorNo[lev], gen);
    for(size_t obs=0;obs<n;obs++) RE_sumztre(obs)+=z_train[lev].row(obs)*current_randomEffects[lev].row(id_train[lev][obs]).transpose();
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
	sdraw[trcnt]=sigma;
	thetadraws[trcnt]=bm.gettheta();
	for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=Offset+bm.f(k);
	for(size_t lev=0;lev<id_train_no;lev++){
	  redraws[lev][trcnt]=current_randomEffects[lev];
	  varcovdraws[lev][trcnt]=current_varianceComponents[lev];
	}
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
  // delete[] u;

#ifndef NoRcpp
  //return list
  Rcpp::List ret;
  //   ret["X"]=X; 
  if(type1sigest) ret["sigma"]=sdraw;
  ret["fhat.train"]=trdraw;
  ret["fhat.test"]=tedraw;
  ret["varcount"]=varcnt;
  ret["varprob"]=varprb;
  ret["re.train"]=redraws;
  ret["re.varcov"]=varcovdraws;
  if(useDart) ret["theta.train"]=thetadraws;

//  ret["u.train"]=udraw;
//  if(sigUpr!=4) ret["sd.u"]=sdudraw;

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
}
#else

#endif

