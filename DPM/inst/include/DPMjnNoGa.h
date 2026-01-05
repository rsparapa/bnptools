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

#ifndef DPMjnNoGa_h
#define DPMjnNoGa_h

#include "DPMjn8.h"

namespace DPM {
  // phi contains mu1, ..., mu6, tau1, ..., tau6
  // prior contains m0, k0, k0_a, k0_b, a0, b0, b0_a, b0_b

  double jnNoGa_F0(SEXP _Y, SEXP _row) {
    Rcpp::NumericVector Y(_Y), phi(_row);
    size_t p=Y.size();
    double result = 1.;
    for(size_t i = 0; i<p; i++) 
      result *= R::dnorm(Y[i], phi[i], sqrt(1./phi[i+p]), 0);
    return result; 
  }

  SEXP jnNoGa_G0tau(size_t p, double a0, double b0, DPM::rng &eng) {
    Rcpp::NumericVector tau(p);
    for(size_t i = 0; i<p; i++) tau[i] = eng.rgamma(a0, b0);
    return Rcpp::wrap(tau);
  }

  SEXP jnNoGa_G0mu(SEXP _m0, SEXP _k0, DPM::rng &eng) {
    Rcpp::NumericVector m0(_m0), k0(_k0); 
    // m0=w
    // precision=k0=B
    size_t p=m0.size();
    Rcpp::NumericVector mu(p);
    for(size_t i = 0; i<p; i++) mu[i] = eng.rnorm(m0[i], sqrt(1./k0[i]));
    return Rcpp::wrap(mu);
  }

  void jnNoGa_P0(size_t row_c, SEXP _Y, SEXP _phi, 
		SEXP _m0, SEXP _k0, double a0, double b0, DPM::rng &eng) {
    Rcpp::NumericMatrix Y(_Y), phi(_phi);
    const int n_c=Y.nrow(), p = Y.ncol();
    Rcpp::NumericVector y, m0(_m0), k0(_k0);
    // a0=r, b0=1/R
    double k1, a1=a0+0.5*n_c, b1, mu_star, tau_star;

    for(size_t i = 0; i<p; i++) {
      y = Y(Rcpp::_, i);
      tau_star = phi(row_c, i+p);
      k1 = k0[i]+n_c*tau_star;
      mu_star=eng.rnorm((k0[i]*m0[i]+tau_star*Rcpp::sum(y))/k1, sqrt(1./k1));
      if(R_finite(mu_star)) phi(row_c, i)=mu_star;
      else mu_star = phi(row_c, i);
      b1 = b0;
      for(size_t j = 0; j<n_c; j++) b1 += 0.5*pow(y[j]-mu_star, 2.); 
      tau_star=eng.rgamma(a1, b1);
      if(R_finite(tau_star)) phi(row_c, i+p)=tau_star;
    }
    return;
  }

  RcppExport SEXP call_jnNoGa(SEXP _Y, SEXP _phi, SEXP _C, SEXP _states, 
			     SEXP _prior, SEXP _hyper, SEXP _mcmc) {
    Rcpp::NumericMatrix Y(_Y);
    Rcpp::IntegerVector states(_states), C(_C);
    Rcpp::NumericMatrix phi(_phi);
    Rcpp::List prior(_prior), hyper(_hyper), mcmc(_mcmc);
    Rcpp::NumericVector m0(prior["m0"]), k0(hyper["k0"]);

    int burn=Rcpp::as<int>(mcmc["burn"]), keep=Rcpp::as<int>(mcmc["keep"]),
      thin=Rcpp::as<int>(mcmc["thin"]), p=phi.ncol(), 
      m=Rcpp::as<int>(prior["m"]);

    double alpha=Rcpp::as<double>(hyper["alpha"]),
      alpha_a=Rcpp::as<double>(prior["alpha.a"]), 
      alpha_b=Rcpp::as<double>(prior["alpha.b"]), 
      a0=Rcpp::as<double>(prior["a0"]), 
      b0=Rcpp::as<double>(hyper["b0"]); 
/*
      m0=Rcpp::as<double>(prior["m0"]), 
      k0=Rcpp::as<double>(hyper["k0"]), 
      k0_a=Rcpp::as<double>(prior["k0.a"]), 
      k0_b=Rcpp::as<double>(prior["k0.b"]), 
      a0=Rcpp::as<double>(prior["a0"]), 
      b0=Rcpp::as<double>(hyper["b0"]), 
      b0_a=Rcpp::as<double>(prior["b0.a"]), 
      b0_b=Rcpp::as<double>(prior["b0.b"]); 
*/

    Rcpp::List GS(keep);

    int l=-burn, h, k, total=keep*thin;

    DPM::arng eng;

    do{
      if(l<0 && (-l*4)%burn == 0) 
	cout << "burn:" << std::abs(-(((-l*4.)/burn)-4.)/4.) << "\n"; 
      else if(l>0 && (l*10)%total == 0) 
	cout << "keep:" << (double)l/total << "\n"; 

/*
      if(m==0) 
	DPM::neal7(_Y, _phi, _C, _states, _prior, _hyper, eng,
		   &DPM::NoGa_F, &DPM::NoGa_G0tau,
		   &DPM::NoGa_G0mu, &DPM::NoGa_P0);
      else */ 
	DPM::jn8(_Y, _phi, _C, _states, _prior, _hyper, eng,
		   &DPM::jnNoGa_F0, &DPM::jnNoGa_G0tau,
		   &DPM::jnNoGa_G0mu, &DPM::jnNoGa_P0);

      if(l>=0 && l%thin == 0) {
	h = (l/thin);

//#define DEBUG
#ifdef DEBUG
	GS[h]=Rcpp::List::create(Rcpp::Named("C")=(Rcpp::clone(C)+1), 
				 Rcpp::Named("phi")=Rcpp::clone(phi), 
				 Rcpp::Named("states")=Rcpp::clone(states),
				 Rcpp::Named("hyper")=Rcpp::clone(hyper));
#else
	k=Rcpp::max(C)+1;

	Rcpp::NumericMatrix Phi(k, p);
	Rcpp::NumericVector States(k);
	for(size_t j=0; j<k; ++j) {
	  States[j]=states[j];
	  Phi.row(j)=phi.row(j);
	}

	GS[h]=Rcpp::List::create(Rcpp::Named("C")=(Rcpp::clone(C)+1), 
				 Rcpp::Named("phi")=Phi,
				 Rcpp::Named("states")=States,
				 Rcpp::Named("hyper")=Rcpp::clone(hyper));
#endif
      }

      l++;

    } while (l<=(keep-1)*thin); 

    return GS;
  }
}

#endif
