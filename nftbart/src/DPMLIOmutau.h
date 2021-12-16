/*
 * Copyright (C) 2021 Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * DPMLIOmutau.h
 *
 * nftbart is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * nftbart is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author contact information
 * Rodney A. Sparapani: rsparapa@mcw.edu
 *
 */

//#include "rn.h"

// phi contains mu, tau
// prior contains m0, k0, k0_a, k0_b, a0, b0, b0_a, b0_b
double DPMLIOmutau_F(double y, double mu, double tau) {
  // N(mu, 1./tau) 
  return R::dnorm(y, mu, pow(tau, -0.5), 0);
  //return sqrt(tau)*exp(-0.5*tau*pow(y-mu, 2.));
}

double DPMLIOmutau_G0tau(double a0, double b0,
			 //double b0_a, double b0_b, 
			 rn &gen) {
//  if(b0==0.) b0=gen.gamma(b0_a, b0_b);
  return gen.gamma(a0, b0);
}

double DPMLIOmutau_G0mu(double tau, double m0, double k0,
			//double k0_a, double k0_b, 
			rn &gen) {
//  if(k0==0.) k0=gen.gamma(k0_a, k0_b);
  return gen.normal(m0, 1./sqrt(k0*tau));
}

//void DPMLIOmutau_P0(size_t row_c, Eigen::VectorXd &Y, SEXP _phi, 
void DPMLIOmutau_P0(size_t row_c, SEXP _Y, SEXP _phi, 
		    double m0, double k0, //double k0_a, double k0_b, 
		    double a0, double b0, //double b0_a, double b0_b, 
		    rn &gen) {
  Rcpp::NumericVector Y(_Y);
  //Eigen::VectorXd Y=Rcpp::as< Eigen::Map<Eigen::VectorXd> >(_Y);
  const int n=Y.size();
  Rcpp::NumericMatrix phi(_phi);

  double /*mu=phi(row_c, 0),*/ 
    k1, tau=phi(row_c, 1), a1=a0+n*0.5, mu_star, tau_star;
  k1=k0+n;

  mu_star=gen.normal((k0*m0+Rcpp::sum(Y))/k1, 1./sqrt(k1*tau));

  if(R_finite(mu_star)) phi(row_c, 0)=mu_star;

  double Ybar=Rcpp::mean(Y);
  Rcpp::NumericVector Y0=Y-Ybar;
  double SSY=Rcpp::sum(Y0*Y0);

  tau_star=gen.gamma(a1, b0+0.5*(SSY+n*k0*pow(Ybar-m0, 2.)/k1));

  if(R_finite(tau_star)) phi(row_c, 1)=tau_star;
/*
  mu_star=gen.normal((k0*m0+Y.sum())/k1, 1./sqrt(k1*tau));

  if(R_finite(mu_star)) phi(row_c, 0)=mu_star;

  Eigen::VectorXd Ybar;

  Ybar.setConstant(n, Y.sum()/n);

  Y -= Ybar;
  Y = Y.array()*Y.array();

  tau_star=gen.gamma(a1, b0+0.5*(Y.sum()+n*k0*pow(Ybar[0]-m0, 2.)/k1));

  if(R_finite(tau_star)) phi(row_c, 1)=tau_star;
*/
}

/*
RcppExport SEXP CallDPMLIOmutau(SEXP _y, SEXP _phi, SEXP _C, SEXP _states, 
				SEXP _prior, SEXP _hyper, SEXP _mcmc) {
  Rcpp::NumericVector y(_y);
  Rcpp::IntegerVector states(_states), C(_C);
  Rcpp::NumericMatrix phi(_phi);
  Rcpp::List prior(_prior), hyper(_hyper), mcmc(_mcmc);

  int burn=Rcpp::as<int>(mcmc["burn"]), keep=Rcpp::as<int>(mcmc["keep"]),
    thin=Rcpp::as<int>(mcmc["thin"]);

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

  const int n=y.size();

  Rcpp::List GS(keep);

  int l=-burn, k;

  rrn gen;

  do{
      DPMLIOneal7(_y, _phi, _C, _states, _prior, _hyper, gen,
		  &DPMLIOmutau_F, &DPMLIOmutau_G0tau,
		  &DPMLIOmutau_G0mu, &DPMLIOmutau_P0);
    
    if(l>=0 && l%thin == 0) {
      int h = (l/thin);

      GS[h]=Rcpp::List::create(Rcpp::Named("C")=(Rcpp::clone(C)+1), 
			       Rcpp::Named("phi")=Rcpp::clone(phi), 
			       Rcpp::Named("states")=Rcpp::clone(states),
			       Rcpp::Named("hyper")=Rcpp::clone(hyper));
    }

    l++;

  } while (l<=(keep-1)*thin); 

  return GS;
}
*/
