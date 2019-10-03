/*
 *  DPM: Dirichlet Process Mixtures With Low Information Omnibus Priors
 *  Copyright (C) 2019 Prakash Laud and Rodney Sparapani
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

#ifndef DPMneal7_h
#define DPMneal7_h

#include "DPM.h"

namespace DPM {
  void neal7(SEXP _Y, SEXP _phi, SEXP _C, SEXP _S, 
	     SEXP _prior, SEXP _hyper, DPM::rng &eng,
	     double (*F)(double y, double mu, double tau),
	     double (*G0tau)(double a0, double b0, DPM::rng &eng),
	     double (*G0mu)(double tau, double m0, double k0, 
			    DPM::rng &eng),
	     void (*P0)(size_t c_row, SEXP _Y, SEXP _phi, 
			double m0, double k0, double a0, double b0, 
			DPM::rng &eng)) {

    Rcpp::NumericMatrix Y(_Y);
    //Rcpp::NumericVector Y(_y);
    Rcpp::IntegerVector S(_S), C(_C);
    Rcpp::NumericMatrix phi(_phi);
    //Eigen::MatrixXd phi=Rcpp::as< Eigen::Map<Eigen::MatrixXd> >(_phi);
    Rcpp::List prior(_prior), hyper(_hyper); 

    const int N=Y.nrow(), M=Y.ncol(), p=phi.ncol(); //p=phi.cols();  

    double 
      alpha=Rcpp::as<double>(hyper["alpha"]), 
      alpha_a=Rcpp::as<double>(prior["alpha.a"]), 
      alpha_b=Rcpp::as<double>(prior["alpha.b"]), 
      m0=Rcpp::as<double>(prior["m0"]), 
      k0=Rcpp::as<double>(hyper["k0"]), 
      k0_a=Rcpp::as<double>(prior["k0.a"]), 
      k0_b=Rcpp::as<double>(prior["k0.b"]), 
      a0=Rcpp::as<double>(prior["a0"]), 
      b0=Rcpp::as<double>(hyper["b0"]), 
      b0_a=Rcpp::as<double>(prior["b0.a"]), 
      b0_b=Rcpp::as<double>(prior["b0.b"]); 

    int k=Rcpp::max(C)+1, // number of states
      k0_draw=Rcpp::as<int>(hyper["k0.draw"]), 
      b0_draw=Rcpp::as<int>(hyper["b0.draw"]), 
      alpha_draw=Rcpp::as<int>(hyper["alpha.draw"]); 

    if(M==3) { // right censoring
      Rcpp::NumericVector prob(k);
      //Eigen::VectorXd prob(k);
      for(size_t i=0; i<N; ++i) if(Y(i, 2)==0.) {
	  for(size_t j=0; j<k; ++j) 
	    prob[j]=VEC(S, j)*
	      R::pnorm(Y(i, 1), phi(j, 0), pow(phi(j, 1), -0.5), 0, 0)/N;
	  int h;
	  h=eng.rcat(prob);
	  if(h==-1) h=i;
	  Y(i, 0)=eng.rtnorm(Y(i, 1), phi(h, 0), pow(phi(h, 1), -0.5));
	} 
    }

    if(k0_draw) {
      double k0_a_post=k0_a+0.5*k, k0_b_post=k0_b;
      for(size_t i=0; i<k; ++i) 
	k0_b_post += 0.5*phi(i, 1)*pow(phi(i, 0)-m0, 2.);
      k0=eng.rgamma(k0_a_post, k0_b_post);
      hyper["k0"]=k0;
    }

    if(b0_draw) {
      double b0_b_post=b0_b;
      for(size_t i=0; i<k; ++i) b0_b_post += phi(i, 1);
      b0=eng.rgamma(b0_a+k*a0, b0_b_post);
      hyper["b0"]=b0;
    }

    bool flag=true, accept;

    for(size_t i=0; i<N; ++i) {
      double y, tau, mu, u;
      bool singleton;
      int c, s, cmax; 
      y=Y(i, 0);
      s=C[i]; 
      singleton=(VEC(S, s)==1); //singleton=(S[s]==1);
      cmax=(k-1);
      u=eng.runif();
      if(singleton) {
	Rcpp::NumericVector prob(k);
	//Eigen::VectorXd prob(k);
	S[s]=0; // if accepted, would leave a "gap" in S if s<k-1
	for(size_t j=0; j<k; ++j) prob[j]=VEC(S, j);
	c=eng.rcat(prob);
	if(c==-1) accept=0;
	else accept=(u<(((N-1)*((*F) (y, phi(c, 0), phi(c, 1))))/
			(alpha*((*F) (y, phi(s, 0), phi(s, 1)))))); 
	if(accept) {
	  C[i]=c;
	  S[c]=VEC(S, c)+1; //S[c]=S[c]+1;
	  // swap states s and k-1 if s<k-1 to prevent a "gap"
	  if(s<cmax) {
	    S[s]=VEC(S, cmax); //S[s]=S[cmax];
	    phi.row(s)=phi.row(cmax);
	    for(size_t j=0; j<N; ++j) if(C[j]==cmax) C[j]=s;
	  }
	  VEC(S, cmax)=0;
	  k=cmax;
	}
	else S[s]=1;
      }
      else {
	c=C[i];
	tau= (*G0tau) (a0, b0, eng);
	mu= (*G0mu) (tau, m0, k0, eng);
	accept=(u<((alpha*((*F) (y, mu, tau)))/
		   ((N-1)*((*F) (y, phi(c, 0), phi(c, 1))))));
	if(accept) {
	  phi(k, 0)=mu;
	  phi(k, 1)=tau;
	  C[i]=k;
	  S[c]=VEC(S, c)-1; //S[c]=S[c]-1;
	  VEC(S, k)=1; //S[k]=1;
	  k++;
	}
      }
#ifdef DEBUG
      if(flag && (Rcpp::sum(S)!=N || Rcpp::max(C)!=(k-1))) { 
	flag=false;
	cout << "singleton:" << singleton << ' ' << "s:" << s << ' '
	     << "k:" << k << ' ' << "c:" << c << ' ' 
	     << "i:" << i << ' ' << "cmax:" << cmax << ' '
	     << "accept:" << accept << '\n';
	//      COUT << "prob:" << prob.sum() << '\n';
	//      for(size_t j=0; j<=k; ++j) COUT << j << ':' << prob[j] << ' ';
	//      COUT << '\n';
	cout << "S:" << Rcpp::sum(S) << '\n';
	for(size_t j=0; j<N; ++j) if(S[j]>0) cout << j << ':' << S[j] << ' ';
	cout << '\n';
	cout << "C:" << Rcpp::max(C)+1 << '\n';
	for(size_t j=0; j<N; ++j) cout << j << ':' << C[j] << ' ';
	cout << '\n';
      }
#endif 
    }

    int s, c;
    double y;
    Rcpp::NumericVector prob(k);
    //Eigen::VectorXd prob(k);
    for(size_t i=0; i<N; ++i) {
      s=C[i];
      if(VEC(S, s)>1) {
	y=Y(i, 0);
	S[s]=S[s]-1;
	for(size_t j=0; j<k; ++j) 
	  prob[j]=(VEC(S, j)/(N-1.))*((*F) (y, phi(j, 0), phi(j, 1)));
	c=eng.rcat(prob);
	if(c==-1) c=s;
#ifdef DEBUG
	if(c<0 || c>=N || Rcpp::sum(prob)<=0.) {
	  cout << "N:" << N << '\n';
	  cout << "y:" << y << '\n';
	  cout << "s:" << s << '\n';
	  cout << "c:" << c << '\n';
	  cout << "p:" << prob << '\n';
	  for(size_t j=0; j<k; ++j) cout << S[j] << ' ';
	  cout << '\n';
	  for(size_t j=0; j<k; ++j) 
	    cout << phi(j, 0) << ' ' << phi(j, 1) << '\n';
	  for(size_t j=0; j<k; ++j) 
	    cout << (*F) (y, phi(j, 0), phi(j, 1)) << ' ';
	  cout << '\n';
	  c=s;
	}
#endif 
	C[i]=c;
	S[c]=VEC(S, c)+1;
      }
    }

#ifdef DEBUG
    if(Rcpp::sum(S)!=N || Rcpp::max(C)!=(k-1)) {
      cout << "Sum:" << Rcpp::sum(S) << '\n';
      cout << "Max:" << Rcpp::max(C)+1 << '\n';
      for(size_t i=0; i<N; ++i) { 
	if(S[i]>0)
	  cout << "S:" << i << ' ' << S[i] << '\n';
	if(C[i]>=k)
	  cout << "C:" << i << ' ' << C[i] << '\n';
      }
    }
#endif

    for(int c=0; c<k; c++) {
      int k_c=VEC(S, c);
#ifdef DEBUG
      if(k_c<1) {
	cout << "Sum:" << Rcpp::sum(S) << '\n';
	cout << "Max:" << Rcpp::max(C)+1 << '\n';
	for(size_t i=0; i<N; ++i) { 
	  if(S[i]<1)
	    cout << "S:" << i << ' ' << S[i] << '\n';
	  if(C[i]>=k)
	    cout << "C:" << i << ' ' << C[i] << '\n';
	}
      }
      else 
#endif
	{
	Rcpp::NumericVector y(k_c);
	//Eigen::VectorXd y(k_c);
	//Eigen::MatrixXd A(k_c, M);
  
	for(int i=0, j=0; i<N && j<k_c; ++i) if(C[i] == c) {
	    y[j]=Y(i, 0);
	    //A.row(j)=Y.row(i);
	    j++;
	  }
    
	(*P0) (c, Rcpp::wrap(y), _phi, m0, k0, a0, b0, eng);
      }
    }

    if(alpha_draw) {
      double eta=eng.rbeta(alpha, N);
      alpha=eng.rgamma(alpha_a+k, alpha_b-log(eta)); 
      hyper["alpha"]=alpha;
    }
  }
}

#endif
