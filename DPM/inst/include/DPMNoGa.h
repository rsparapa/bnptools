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

#ifndef DPMNoGa_h
#define DPMNoGa_h

#include "DPMneal7.h"
#include "DPMneal8.h"

namespace DPM {
  // phi contains mu, tau
  // prior contains m0, k0, k0_a, k0_b, a0, b0, b0_a, b0_b
  double NoGa_F(double y, double mu, double tau) {
    // N(mu, 1./tau) 
    return R::dnorm(y, mu, pow(tau, -0.5), 0);
    //return sqrt(tau)*exp(-0.5*tau*pow(y-mu, 2.));
  }

  double NoGa_G0tau(double a0, double b0, DPM::rng &eng) {
    return eng.rgamma(a0, b0);
  }

  double NoGa_G0mu(double tau, double m0, double k0, DPM::rng &eng) {
    return eng.rnorm(m0, 1./sqrt(k0*tau));
  }

  //void NoGa_P0(size_t row_c, Eigen::VectorXd &Y, SEXP _phi, 
  void NoGa_P0(size_t row_c, SEXP _Y, SEXP _phi, 
		double m0, double k0, double a0, double b0, DPM::rng &eng) {
    Rcpp::NumericVector Y(_Y);
    //Eigen::VectorXd Y=Rcpp::as< Eigen::Map<Eigen::VectorXd> >(_Y);
    const int n=Y.size();
    Rcpp::NumericMatrix phi(_phi);

    double /*mu=phi(row_c, 0),*/ 
      k1, tau=phi(row_c, 1), a1=a0+n*0.5, mu_star, tau_star;
    k1=k0+n;

    mu_star=eng.rnorm((k0*m0+Rcpp::sum(Y))/k1, 1./sqrt(k1*tau));

    if(R_finite(mu_star)) phi(row_c, 0)=mu_star;

    double Ybar=Rcpp::mean(Y);
    Rcpp::NumericVector Y0=Y-Ybar;
    double SSY=Rcpp::sum(Y0*Y0);

    tau_star=eng.rgamma(a1, b0+0.5*(SSY+n*k0*pow(Ybar-m0, 2.)/k1));

    if(R_finite(tau_star)) phi(row_c, 1)=tau_star;
/*
    Eigen::VectorXd Ybar;

    Ybar.setConstant(n, Y.sum()/n);

    Y -= Ybar;
    Y = Y.array()*Y.array();

    tau_star=eng.rgamma(a1, b0+0.5*(Y.sum()+n*k0*pow(Ybar[0]-m0, 2.)/k1));

    if(R_finite(tau_star)) phi(row_c, 1)=tau_star;
*/
  }

  RcppExport SEXP call_NoGa(SEXP _Y, SEXP _phi, SEXP _C, SEXP _states, 
			     SEXP _prior, SEXP _hyper, SEXP _mcmc) {
    Rcpp::NumericMatrix Y(_Y);
    //Rcpp::NumericVector y(_y);
    Rcpp::IntegerVector states(_states), C(_C);
    Rcpp::NumericMatrix phi(_phi);
    Rcpp::List prior(_prior), hyper(_hyper), mcmc(_mcmc);

    int burn=Rcpp::as<int>(mcmc["burn"]), keep=Rcpp::as<int>(mcmc["keep"]),
      thin=Rcpp::as<int>(mcmc["thin"]), p=phi.ncol(), 
      m=Rcpp::as<int>(prior["m"]);

    double alpha=Rcpp::as<double>(hyper["alpha"]),
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

    Rcpp::List GS(keep);

    int l=-burn, h, k, total=keep*thin;

    DPM::arng eng;

    do{
      if(l<0 && (-l*4)%burn == 0) 
	cout << "burn:" << std::abs(-(((-l*4.)/burn)-4.)/4.) << "\n"; 
      else if(l>0 && (l*10)%total == 0) 
	cout << "keep:" << (double)l/total << "\n"; 

      if(m==0) 
	DPM::neal7(_Y, _phi, _C, _states, _prior, _hyper, eng,
		   &DPM::NoGa_F, &DPM::NoGa_G0tau,
		   &DPM::NoGa_G0mu, &DPM::NoGa_P0);
      else  
	DPM::neal8(_Y, _phi, _C, _states, _prior, _hyper, eng,
		   &DPM::NoGa_F, &DPM::NoGa_G0tau,
		   &DPM::NoGa_G0mu, &DPM::NoGa_P0);

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
	//Phi=phi(Rcpp::Range(0, k-1), Rcpp::Range(0, phiM-1));
	Rcpp::NumericVector States(k);
	//States=states(Rcpp::Range(0, k-1));
	for(size_t j=0; j<k; ++j) {
	  States[j]=states[j];
	  Phi.row(j)=phi.row(j);
	}

/*
	Eigen::MatrixXd Phi;
Phi=Rcpp::as< Eigen::Map<Eigen::MatrixXd> >(Rcpp::wrap(Rcpp::clone(phi)));
	Phi.conservativeResize(k, phiM);
	Eigen::VectorXd States;
States=Rcpp::as< Eigen::Map<Eigen::VectorXd> >(Rcpp::wrap(Rcpp::clone(states)));
	States.conservativeResize(k);
*/

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
