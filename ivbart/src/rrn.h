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

#ifndef RRN_H
#define RRN_H

/*
#ifdef Rcpp_hpp
using R::rchisq;
#else
extern "C" {
#include <R.h>
#include <Rmath.h>
};
#endif
*/
#include <Rcpp.h>
using R::rchisq;
#include "rn.h"

#define printf Rprintf
#define cout Rcpp::Rcout

#ifndef PI
#define PI 3.141592653589793238462643383280
#endif
// log (2 pi)
#define LTPI 1.837877066409345483560659472811 
// sqrt(2*pi)
#define RTPI 2.506628274631000502415765284811

class rrn: public rn
{
public:
//constructor
   rrn():df(1),alpha(1) {}
//virtual
   virtual ~rrn() {}
   virtual double normal() {return norm_rand();}
   virtual double uniform() { return unif_rand();}
   virtual double chi_square() {return rchisq((double)df);}
   virtual void set_df(int df) {this->df=df;}
   virtual double gamma() { return R::rgamma(alpha,1.0);}
   virtual void set_alpha(double alpha) {this->alpha=alpha;} //set alpha shape parameter for gamma
   virtual double exp() {return exp_rand();}
//get,set
   int get_df() {return df;}
private:
   int df;
   double alpha; //shape for gamma
};

#endif
