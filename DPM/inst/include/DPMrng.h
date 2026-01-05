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

#ifndef DPMrng_h
#define DPMrng_h

//#include <RcppEigen.h>
#include <Rcpp.h>

namespace DPM {
  //pure virtual base class for random number generator
  class rng {
  public:
    virtual double rbeta(double a, double b) = 0;
    virtual int    rbinom(int n, double p) = 0;
    virtual int    rcat(Rcpp::NumericVector prob) = 0;
    //virtual int    rcat(Eigen::VectorXd prob) = 0;
    virtual double rchisq(double df) = 0;
    virtual double rexp(double rate=1.) = 0; 
    virtual double rgamma(double shape, double rate=1., double small=0.1) = 0;
    virtual double rlgamma(double shape) = 0; 
    virtual double rnorm(double mu=0., double sd=1.) = 0; 
    virtual double rtnorm(double tau, double mu=0., double sd=1.) = 0; 
    virtual double runif(double a=0., double b=1.) = 0; 
    virtual ~rng() {}
  };

  //abstract random number generator based on R/Rcpp
  class arng: public rng
  {
  public:
    arng() {}
    ~arng() {}
    virtual double rbeta(double a, double b) {return R::rbeta(a, b);} 
    virtual int    rbinom(int n, double p) {return R::rbinom(n, p);} 

    virtual int    rcat(Rcpp::NumericVector _p) {
      double c=Rcpp::sum(_p), d=Rcpp::min(_p);
      if(c==0. || d<0.) {
#ifdef DEBUG
	cout << "rcat returning -1\n";
	cout << _p << '\n';
#endif
	return -1;
      }
      int K=_p.size();
      Rcpp::NumericVector p = _p/c;
      Rcpp::IntegerVector x(K);
      R::rmultinom(1, &p[0], K, &x[0]);
      for(int j=0; j<K; ++j) if(x[j]==1) return j;
      return -1; // to suppress warning on previous return above
    }

/*
    virtual int    rcat(Eigen::VectorXd prob) {
      // return a draw from Multinomial(1, prob)
      // e.g. c(0:(h-1)) %*% rmultinom(1, 1, prob)
      int j, K=prob.size(), rn[K];
      double p[K];
      Eigen::Map<Eigen::VectorXd> mprob(&p[0], K);
      Eigen::Map<Eigen::VectorXi> x(&rn[0], K);
      Eigen::VectorXi z(K);
      prob /= prob.sum();
      mprob=prob;
      R::rmultinom(1, &p[0], K, &rn[0]);
      //  for(j=0; j<K; ++j) z[j]=j;
      //  return z.dot(x);
      for(j=0; j<K; ++j) if(rn[j]==1) return j;
    }
*/
    virtual double rchisq(double df) {return R::rchisq(df);}
    virtual double rexp(double rate=1.) {return R::rexp(rate);}

    virtual double rgamma(double shape, double rate=1., double small=0.1) {
// in R, the argument order is rgamma(shape, rate)
// and the same is true for rng::rgamma
// E[y]=shape/rate
// V[y]=shape/(rate^2)
// rate=b0=1/R 
// however in Rmath and Rcpp it is defined differently
// y=rgamma(shape, scale) where scale=1/rate
// E[y]=shape*scale 
// V[y]=shape*scale^2
      if(shape<=small) {
	double z;
	do { z=exp(this->rlgamma(shape)-log(rate)); } while (z==0.);
	//do { z=exp(this->rlgamma(shape))/rate; } while (z==0.);
	return z;
      }
      //if(shape<0.01) return ::exp(this->rlgamma(shape))/rate;
      else return R::rgamma(shape, 1.)/rate; 
    } 

    virtual double rlgamma(double shape) {
      double y=log(R::rgamma(shape+1., 1.)), z=log(this->runif())/shape;
      return y+z; 
    }

    virtual double rnorm(double mu=0., double sd=1.) 
    {return R::rnorm(mu, sd);}

    virtual double rtnorm(double tau, double mean=0., double sd=1.) {
      double x, z, lambda;

      tau = (tau - mean)/sd;

      if(tau<=0.) {
	/* draw until we get to the right of tau */
	do { z=this->rnorm(); } while (z < tau);
      }
      else {
	/* optimal exponential rate parameter */
	lambda = 0.5*(tau + sqrt(pow(tau, 2.) + 4.));

	/* rejection sampling */
	do { z = this->rexp()/lambda + tau; } 
	while (this->runif() > exp(-0.5*pow(z - lambda, 2.)));
      }

      /* put x back on the right scale */
      x = z*sd + mean;

      return(x);
    }

    virtual double runif(double a=0., double b=1.) {return R::runif(a, b);}
  
    //virtual double weibull(double shape, double scale) 
    //  {return scale*pow(-R::log(unif_rand()), 1./shape);}
    /*
      virtual double rnormt(double a, double b) {
      double L=R::pnorm(a, 0., 1., 1, 0), R=R::pnorm(b, 0., 1., 1, 0), 
      U=R::unif_rand();
      return R::qnorm(L+U*(R-L), 0., 1., 1, 0);}
      virtual double rnormt(double a, double b, double mu, double sd)
      {return mu+sd*rnormt((a-mu)/sd, (b-mu)/sd);}
    */

    Rcpp::RNGScope get_state(void) { return RNGstate; }
    void set_state(Rcpp::RNGScope RNGstate) { this->RNGstate=RNGstate; }
    
  private:
    Rcpp::RNGScope RNGstate;
    //std::vector<double> wts;
  };
}

RcppExport SEXP call_rcat(SEXP _prob) {
  DPM::arng eng;
  Rcpp::IntegerVector z(1);
  Rcpp::NumericVector prob(_prob);
/*
  size_t K = prob.size();
  Eigen::Map<Eigen::VectorXd> mprob(&prob[0], K);
  z=eng.rcat(mprob);
*/
  z=eng.rcat(prob);
  return Rcpp::wrap(z);
}

RcppExport SEXP call_rgamma(SEXP _n, SEXP _shape, SEXP _rate) {
  DPM::arng eng;
  size_t N = Rcpp::as<int>(_n);
  Rcpp::NumericVector z(N), shape(_shape), rate(_rate);
  size_t S=shape.size(), R=rate.size();
  for(size_t i=0; i<N; ++i) z[i]=eng.rgamma(shape[i%S], rate[i%R]);
  return Rcpp::wrap(z);
}

RcppExport SEXP call_rlgamma(SEXP _n, SEXP _shape) {
  DPM::arng eng;
  size_t N = Rcpp::as<int>(_n);
  Rcpp::NumericVector z(N), shape(_shape);
  size_t S=shape.size();
  for(size_t i=0; i<N; ++i) z[i]=eng.rlgamma(shape[i%S]);
  return Rcpp::wrap(z);
}

RcppExport SEXP call_rtnorm(SEXP _n, SEXP _tau, SEXP _mean, SEXP _sd) {
  DPM::arng eng;
  size_t N = Rcpp::as<int>(_n);
  Rcpp::NumericVector z(N), tau(_tau), mean(_mean), sd(_sd);
  size_t T=tau.size(), M=mean.size(), S=sd.size();
  for(size_t i=0; i<N; ++i) z[i]=eng.rtnorm(tau[i%T], mean[i%M], sd[i%S]);
  return Rcpp::wrap(z);
}

#endif
