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

#include "dart.h"

//--------------------------------------------------
//constructor
dart::dart(){}
dart::~dart(){}

//--------------------------------------------------
//operators
dart& dart::operator=(const dart& drt)
{
   return *this;
}

//--------------------------------------------------
//draw variable splitting probabilities from Dirichlet (Linero, 2018)
void dart::draw_s(rn& gen){
// Now draw s, the vector of splitting probabilities  
  std::vector<double> _theta(p);
  for(size_t j=0;j<p;j++) _theta[j]=theta+(double)nv[j];
  //gen.set_alpha(_theta);
  lpv=gen.log_dirichlet(_theta);
  for(size_t j=0;j<p;j++) pv[j]=std::exp(lpv[j]);
}

//--------------------------------------------------
//draw Dirichlet sparsity parameter from posterior using grid
void dart::draw_theta0(rn& gen){
  // Draw sparsity parameter theta_0 (Linero calls it alpha); see Linero, 2018
  // theta / (theta + rho ) ~ Beta(a,b)
  // Set (a=0.5, b=1) for sparsity
  // Set (a=1, b=1) for non-sparsity
  // rho = p usually, but making rho < p increases sparsity
  if(!const_theta){
    if(thetaDraw==1) {
      double sumlpv=0.;
      size_t p=lpv.size();     
      for(size_t j=0;j<p;j++) 
	sumlpv+=lpv[j];
      std::vector<double> lambda_g (1000,0.);
      std::vector<double> theta_g (1000,0.);
      std::vector<double> lwt_g (1000,0.);
      std::vector<double> post_g (1000,0.);
      for(size_t k=0;k<1000;k++){
	lambda_g[k]=(double)(k+1)/1001.;
	theta_g[k]=(lambda_g[k]*rho)/((double)p*(1.-lambda_g[k]));
	double theta_log_lik=lgamma((double)p*theta_g[k])-(double)p*lgamma(theta_g[k])+(theta_g[k]-1.)*sumlpv;
	double theta_log_prior=(a-1.)*std::log((double)p*theta_g[k]/rho)+(-a-b)*std::log(1.+(double)p*theta_g[k]/rho);
//	cout << "SLP: " << sumlpv << "\nTLL: " << theta_log_lik << "\nBLP: " << theta_log_prior << '\n';
	lwt_g[k]=theta_log_lik+theta_log_prior;   
//	cout << theta_log_prior << '\n';
      }
      double mx=lwt_g[0];
      size_t mx_k=0;
      double sum_lwt=0.;
      for(size_t k=1;k<1000;k++) {
	if(lwt_g[k]>mx) {
	  mx=lwt_g[k];
	  mx_k=k;
	}
      }
      for(size_t k=0;k<1000;k++) {
	lwt_g[k]=std::exp(lwt_g[k]-mx);
	sum_lwt+=lwt_g[k];
      }
      for(size_t k=0;k<1000;k++) {
	lwt_g[k]=lwt_g[k]/sum_lwt;
      }
      gen.set_wts(lwt_g);
      theta=theta_g[gen.discrete()];
    }
    else {
      double p=(double)lpv.size();
      double sumlpv=0.;
      double z,y,w,L,R,J,K,m;
      for(size_t j=0;j<(size_t)p;j++) sumlpv+=lpv[j];
      z=lgamma(theta*p)-p*lgamma(theta)+(theta-1)*sumlpv+(-a-b)*log(1+p*theta/rho)+(a-1)*log(p*theta/rho)-::Rf_lbeta(a,b)-gen.gamma(1,1);
      w=.2;
      m=100.;
      L=theta-w*gen.uniform();
      if(L<0) L=0.;
      R=L+w;
      J=std::floor(m*gen.uniform());
      K=(m-1)-J;
      while(J>0&&z<lgamma(L*p)-p*lgamma(L)+(L-1)*sumlpv+(-a-b)*log(1+p*L/rho)+(a-1)*log(p*L/rho)-::Rf_lbeta(a,b)){
	L=L-w;
	if(L<0) L=0.;
	J=J-1;
      }
      while(K>0&&z<lgamma(R*p)-p*lgamma(R)+(R-1)*sumlpv+(-a-b)*log(1+p*R/rho)+(a-1)*log(p*R/rho)-::Rf_lbeta(a,b)){
	R=R+w;
	K=K-1;
      }
      theta=gen.uniform()*(R-L)+L;
      while(theta<0||z>lgamma(theta*p)-p*lgamma(theta)+(theta-1)*sumlpv+(-a-b)*log(1+p*theta/rho)+(a-1)*log(p*theta/rho)-::Rf_lbeta(a,b)) {
	theta=gen.uniform()*(R-L)+L;
      }
    }
  }
}
