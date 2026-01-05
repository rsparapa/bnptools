/*
 *  DPM: Dirichlet Process Mixtures With Low Information Omnibus Priors
 *  Copyright (C) 2019-2026 Prakash Laud and Rodney Sparapani
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

#ifndef DPMjn8_h
#define DPMjn8_h

#include "DPM.h"

namespace DPM {
  void jn8(SEXP _Y, SEXP _phi, SEXP _C, SEXP _S, 
	     SEXP _prior, SEXP _hyper, DPM::rng &eng,
	     double (*F0)(SEXP _Y, SEXP _row),
	     SEXP (*G0tau)(size_t p, double a0, double b0, DPM::rng &eng),
	     SEXP (*G0mu)(SEXP _m0, SEXP _k0, DPM::rng &eng),
	     void (*P0)(size_t row_c, SEXP _Y, SEXP _phi, 
			SEXP _m0, SEXP _k0, double a0, double b0, 
			DPM::rng &eng)) {

    Rcpp::NumericMatrix Y(_Y), phi(_phi);
    Rcpp::IntegerVector S(_S), C(_C);
    Rcpp::List prior(_prior), hyper(_hyper); 
    Rcpp::NumericVector m0(prior["m0"]), k0(hyper["k0"]);

    const int N=Y.nrow(), p=Y.ncol(), m=Rcpp::as<int>(prior["m"]);  

    double
      alpha=Rcpp::as<double>(hyper["alpha"]),
      alpha_a=Rcpp::as<double>(prior["alpha.a"]), 
      alpha_b=Rcpp::as<double>(prior["alpha.b"]), 
      a0=Rcpp::as<double>(prior["a0"]), 
      b0=Rcpp::as<double>(hyper["b0"]);
/*
      m0=Rcpp::as<double>(prior["m0"]), 
      k0=Rcpp::as<double>(hyper["k0"]), 
      k0_a=Rcpp::as<double>(prior["k0.a"]), 
      k0_b=Rcpp::as<double>(prior["k0.b"]), 
      b0_a=Rcpp::as<double>(prior["b0.a"]), 
      b0_b=Rcpp::as<double>(prior["b0.b"]); 
*/

    int k=Rcpp::max(C)+1, // number of states
      alpha_draw=Rcpp::as<int>(hyper["alpha.draw"]); 
/*
      k0_draw=Rcpp::as<int>(hyper["k0.draw"]), 
      b0_draw=Rcpp::as<int>(hyper["b0.draw"]), 

    if(M==3) { // right censoring
      Rcpp::NumericVector prob(k);
      for(size_t i=0; i<N; ++i) if(Y(i, 2)==0.) {
	  int h;
	  for(size_t j=0; j<k; ++j) 
	    prob[j]=VEC(S, j)*
	      R::pnorm(Y(i, 1), phi(j, 0), pow(phi(j, 1), -0.5), 0, 0)/N;
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
*/

    bool singleton, flag=true;

    for(size_t i=0; i<N; ++i) {
      int c, h, s;
      h=k+m;
      Rcpp::NumericVector prob(h);
      c=C[i];
      s=VEC(S, c);
      singleton=(s==1);
      S[c]=s-1;

      double den;
      Rcpp::NumericVector y(p), mu(p), tau(p);
      y=Y(i, Rcpp::_);
      den=N-1.+alpha;
      for(size_t j=0; j<k; ++j) 
	prob[j]=(VEC(S, j)/den)*((*F0) (y, Rcpp::wrap(phi(j, Rcpp::_))));
      for(size_t j=k; j<h; ++j) {
	if(singleton && j==k) phi.row(k)=phi.row(c);
	else {
	  mu = (*G0mu) (m0, k0, eng);
	  tau = (*G0tau) (p, a0, b0, eng);
	  for(size_t l=0; l<p; l++) {
	    phi(j, l)=mu[l];
	    phi(j, l+p)=tau[l];
	  }
	}
	prob[j]=(alpha/(m*den))*((*F0) (y, Rcpp::wrap(phi(j, Rcpp::_))));
      }

      int j;
      j=eng.rcat(prob);

      if(j==-1) j=c; // otherwise C/S will be out of synch, but this is bad
      else if(singleton) {
#ifdef DEBUG 
	if(j==c) {
	  // this CANNOT happen since S[c]=0!?!
	  cout << "j == c:" << c << "\n";
	  cout << prob << '\n';
	  return;
	}
	else 
#endif
	if(j>=k) {
	  if(j>k) phi.row(c)=phi.row(j);
	  j=c;
	}
	else { // singleton cluster c replaced by c+1 unless k-1
	  if(c<(k-1)) {
	    if(j>c) j--;
	    for(size_t l=c; l<(k-1); ++l) {
	      phi.row(l)=phi.row(l+1);
	      VEC(S, l)=VEC(S, l+1);
	    }
	    S[k-1]=0;
	    for(size_t l=0; l<N; ++l) if(C[l]>c) C[l]=C[l]-1;
	  }
	  k--;
	}
      }
      else if(j>=k) {
	if(j>k) {
	  phi.row(k)=phi.row(j);
	  j=k;
	}
	k++;
      }

      C[i]=j;
      S[j]=VEC(S, j)+1;

#ifdef DEBUG
      if(flag && (Rcpp::sum(S)!=N || Rcpp::max(C)!=(k-1))) { 
	flag=false;
	cout << "singleton:" << singleton << ' ' 
	     << "k:" << k << ' ' << "c:" << c << ' ' 
	     << "j:" << j << ' ' << "i:" << i << '\n';
	//      COUT << "prob:" << prob.sum() << '\n';
	//      for(size_t j=0; j<=k; ++j) COUT << j << ':' << prob[j] << ' ';
	//      COUT << '\n';
	cout << "S:" << Rcpp::sum(S) << '\n';
	for(size_t j=0; j<N; ++j) if(S[j]>0) cout << j << ':' << S[j] << ' ';
	cout << '\n';
	cout << "C:" << Rcpp::max(C) << '\n';
	for(size_t j=0; j<N; ++j) cout << j << ':' << C[j] << ' ';
	cout << '\n';
      }
#endif
    }

#ifdef DEBUG 
    if(Rcpp::sum(S)!=N || Rcpp::max(C)!=(k-1)) {
      cout << "Sum:" << Rcpp::sum(S) << '\n';
      cout << "Max:" << Rcpp::max(C) << '\n';
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
	cout << "Max:" << Rcpp::max(C) << '\n';
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
	Rcpp::NumericMatrix y(k_c, p);
  
	for(int i=0, j=0; i<N && j<k_c; ++i) if(C[i] == c) {
	    y.row(j)=Y.row(i);
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
