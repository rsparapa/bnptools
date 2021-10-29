// rn.h: Random number generator virtual class for using BART in R package.
// Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
//
// This file is part of hbart.
//
// hbart is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// hbart is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author contact information
// Matthew T. Pratola: mpratola@gmail.com
// Robert E. McCulloch: robert.e.mculloch@gmail.com
// Hugh A. Chipman: hughchipman@gmail.com


#ifndef GUARD_rn
#define GUARD_rn

#include <R.h>
#include <Rmath.h>

//pure virtual base class for random numbers
class rn
{
public:
   virtual int bin(int n, double p) = 0;
   virtual double beta(double a, double b) = 0;
   virtual double chi_square(double _df) = 0;
   virtual double chi_square() = 0; //chi-square
   virtual double exp() = 0;
   virtual double gamma(double shape, double rate=1., double small=0.1) = 0;
   virtual double log_gamma(double shape) = 0; 
   virtual double normal() = 0; //standard normal
   virtual double normal(double mu, double sd) = 0; 
   virtual int rcat(Rcpp::NumericVector prob) = 0;
   virtual double rtnorm(double tau, double mu=0., double sd=1.) = 0; 
   virtual void set_df(int df) = 0; //set df for chi-square
   virtual double uniform() = 0; //uniform(0,1)
   virtual ~rn() {}
};

class rrn: public rn
{
public:
//constructor
   rrn():df(1) {}
//virtual
   virtual ~rrn() {}
   virtual double normal() {return norm_rand();}
   virtual double normal(double mu, double sd) 
   {return mu+sd*R::norm_rand();}
   virtual double uniform() { return unif_rand();}
   virtual double chi_square(double _df) {return R::rchisq(_df);}
   virtual double chi_square() {return R::rchisq((double)df);}
   virtual double exp() {return exp_rand();}
   virtual void set_df(int df) {this->df=df;}
   virtual int bin(int n, double p) {return R::rbinom(n, p);}
   virtual double beta(double a, double b) {return R::rbeta(a, b);} 

   virtual double log_gamma(double shape) {
    double y=log(R::rgamma(shape+1., 1.)), z=log(this->uniform())/shape;
    return y+z; 
   }

   virtual double gamma(double shape, double rate=1., double small=0.1) {
      if(shape<=small) {
	double z;
	do { z=::exp(this->log_gamma(shape)-log(rate)); } while (z==0.);
	return z;
      }
      else return R::rgamma(shape, 1.)/rate; 
   } 

   virtual int rcat(Rcpp::NumericVector _p) {
      double c=Rcpp::sum(_p), d=Rcpp::min(_p);
      if(c==0. || d<0.) {
//#ifdef DEBUG
	COUT << "rcat returning -1\n";
	COUT << _p << '\n';
//#endif
	return -1;
      }
      int K=_p.size();
      Rcpp::NumericVector p = _p/c;
      Rcpp::IntegerVector x(K);
      R::rmultinom(1, &p[0], K, &x[0]);
      for(int j=0; j<K; ++j) if(x[j]==1) return j;
      return -2;
    }

    virtual double rtnorm(double tau, double mean=0., double sd=1.) {
      double x, z, lambda;

      tau = (tau - mean)/sd;

      if(tau<=0.) {
	/* draw until we get to the right of tau */
	do { z=this->normal(); } while (z < tau);
      }
      else {
	/* optimal exponential rate parameter */
	lambda = 0.5*(tau + sqrt(pow(tau, 2.) + 4.));

	/* rejection sampling */
	do { z = this->exp()/lambda + tau; } 
	while (this->uniform() > ::exp(-0.5*pow(z - lambda, 2.)));
      }

      /* put x back on the right scale */
      x = z*sd + mean;

      return(x);
    }

//get,set
   int get_df() {return df;}
private:
   int df;
};

#endif
