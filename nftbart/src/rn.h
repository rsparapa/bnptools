/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         Hugh A. Chipman and Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * rn.h
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
 * Matthew T. Pratola: mpratola@gmail.com
 * Robert E. McCulloch: robert.e.mculloch@gmail.com
 * Hugh A. Chipman: hughchipman@gmail.com
 * Rodney A. Sparapani: rsparapa@mcw.edu
 *
 */

#ifndef GUARD_rn
#define GUARD_rn

//#define PROFILER "profiler.log"

#ifdef PROFILER
#include <gperftools/profiler.h>
#endif

#include <sstream>
#include <string>

//#define DEBUG
//#define VEC(x, i) x(i)
#define VEC(x, i) x[i]

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef NotInR
#include <Rcpp.h>
//#include <RcppEigen.h>
#define COUT Rcpp::Rcout 
#define cout Rcpp::Rcout 
#else
#define COUT std::cout 
using std::cout;
#endif

using std::endl;
const double Inf=1./0.;

//pure virtual base class for random numbers
class rn
{
public:
   virtual double normal() = 0; //standard normal
   virtual double normal(double mu, double sd) = 0; 
   virtual double uniform() = 0; //uniform(0,1)
   virtual double chi_square(double _df) = 0;
   virtual double chi_square() = 0; //chi-square
   virtual void set_df(int df) = 0; //set df for chi-square
   virtual double exp() = 0; 
   virtual double log_gamma(double shape) = 0; 
   virtual double gamma(double shape, double rate=1., double small=0.1) = 0;
   virtual double gamma() = 0; 
   virtual void set_alpha(double alpha) = 0; //set df gamma
   virtual int rcat(Rcpp::NumericVector prob) = 0;
   //virtual int rcat(Eigen::VectorXd prob) = 0;
   virtual int bin(int n, double p) = 0;
   virtual double beta(double a, double b) = 0;
   virtual double rtnorm(double tau, double mu=0., double sd=1.) = 0; 
// virtual double weibull(double shape, double scale) = 0;
// virtual double rnormt(double a, double b) = 0;
// virtual double rnormt(double a, double b, double mu, double sd) = 0;
   virtual ~rn() {}
};

//abstract random number generator based on R/Rcpp
class rrn: public rn
{
 public:
  //constructor
 rrn():df(0), alpha(0.) {}
  //virtual
  virtual ~rrn() {}
  virtual double normal() {return R::norm_rand();}
  virtual double normal(double mu, double sd) 
  {return mu+sd*R::norm_rand();}
  virtual double uniform() {return R::unif_rand();}
  virtual double chi_square(double _df) {return R::rchisq(_df);}
  virtual double chi_square() {return R::rchisq((double)this->df);}
  virtual void set_df(int _df) {this->df=_df;}
  virtual double exp() {return R::exp_rand();}
  virtual double log_gamma(double shape) {
    double y=log(R::rgamma(shape+1., 1.)), z=log(this->uniform())/shape;
    return y+z; 
  }

/*
  virtual double gamma(double shape, double rate) {
    if(shape<0.01) return ::exp(this->log_gamma(shape))/rate;
    else return R::rgamma(shape, 1.)/rate; 
  } 
*/

    virtual double gamma(double shape, double rate=1., double small=0.1) {
      if(shape<=small) {
	double z;
	do { z=::exp(this->log_gamma(shape)-log(rate)); } while (z==0.);
	//do { z=::exp(this->log_gamma(shape))/rate; } while (z==0.);
	return z;
      }
      else return R::rgamma(shape, 1.)/rate; 
    } 

  virtual double gamma() {return this->gamma(alpha, 1.);}
  virtual void set_alpha(double _alpha) {this->alpha=_alpha;} 
  int get_df() {return df;}
  double get_alpha() {return alpha;}

    virtual int rcat(Rcpp::NumericVector _p) {
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
#ifdef DEBUG
	cout << "rcat returning -2\n";
	cout << x << '\n';
#endif
      return -2;
    }

/*
  virtual int rcat(Eigen::VectorXd prob) {
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

  virtual int bin(int n, double p) {return R::rbinom(n, p);} 
  virtual double beta(double a, double b) {return R::rbeta(a, b);} 

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
 private:
  //std::vector<double> wts;
  Rcpp::RNGScope RNGstate;
  int df;
  double alpha;
};

/*
void DPMLIOneal7(SEXP _y, SEXP _phi, SEXP _C, SEXP _S, 
		 SEXP _prior, SEXP _hyper, rn &gen,
		 double (*F)(double y, double mu, double tau),
		 double (*G0tau)(double a0, double b0, rn &gen),
		 double (*G0mu)(double tau, double m0, double k0, rn &gen),
		 void (*P0)(size_t c_row, SEXP _Y, SEXP _phi, 
			    double m0, double k0, double a0, double b0, 
			    rn &gen));
*/

void DPMLIOneal8(SEXP _y, SEXP _phi, SEXP _C, SEXP _S, 
		 SEXP _prior, SEXP _hyper, rn &gen,
		 double (*F)(double y, double mu, double tau),
		 double (*G0tau)(double a0, double b0, rn &gen),
		 double (*G0mu)(double tau, double m0, double k0, rn &gen),
		 void (*P0)(size_t c_row, SEXP _Y, SEXP _phi, 
			    double m0, double k0, double a0, double b0, 
			    rn &gen));

double DPMLIOmutau_F(double y, double mu, double tau);

double DPMLIOmutau_G0tau(double a0, double b0, rn &gen);

double DPMLIOmutau_G0mu(double tau, double m0, double k0, rn &gen);

void DPMLIOmutau_P0(size_t row_c, SEXP _Y, SEXP _phi, 
		    double m0, double k0, double a0, double b0, rn &gen);

/*
void DPMtauneal8(SEXP _y, SEXP _phi, SEXP _C, SEXP _S, 
		 SEXP _prior, SEXP _hyper, rn &gen,
		 double (*F)(double y, double tau),
		 double (*G0)(double a0, double b0, rn &gen),
		 void (*P0)(size_t c_row, SEXP _Y, SEXP _phi, 
			    double a0, double b0, rn &gen));

double DPMtau_F(double y, double tau);

double DPMtau_G0(double a0, double b0, rn &gen);

void DPMtau_P0(size_t row_c, SEXP _Y, SEXP _phi, 
		    double a0, double b0, rn &gen);
*/

#endif
