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

#ifndef DPM_h
#define DPM_h

#ifdef MATHLIB_STANDALONE
#define NoRcpp
#else
#define RNG_Rcpp
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using std::endl;
const double Inf=1./0.;

#ifdef _OPENMP
#include <omp.h>
#endif

#define VEC(x, i) x[i]

#ifdef NoRcpp

#include <stdio.h> // for printf

using std::cout;

#define PI 3.141592653589793238462643383280

#else // YesRcpp

// for debugging only
//#define DEBUG

#ifdef DEBUG
#ifdef NDEBUG
#undef NDEBUG
#endif

#ifdef VEC
#undef VEC
#endif

#define VEC(x, i) x(i)
#endif

#include <Rcpp.h>
//#include <RcppEigen.h>

#define printf Rprintf
#define cout Rcpp::Rcout

#endif

// log(2*pi)
#define LTPI 1.837877066409345483560659472811

#include "DPMrng.h"

namespace DPM {
/*
  template <typename T>
    class floor : public std::unary_function<T,T> {
  public:
    T operator()( T t) const { return std::floor(t) ; }
  } ;

  template <typename T>
    class ceiling : public std::unary_function<T,T> {
  public:
    T operator()( T t) const { return std::ceil(t) ; }
  } ;
*/

/*
  Rcpp::IntegerVector concat(Rcpp::IntegerVector x,
			     Rcpp::IntegerVector y) {
    size_t N=x.size(), M=y.size();
    Rcpp::IntegerVector z(M+N);
    for(size_t i=0; i<N; ++i) z[i]=x[i];
    for(size_t i=0; i<M; ++i) z[N+i]=y[i];
    return z;
  }
*/

  Rcpp::IntegerVector floor(Rcpp::NumericVector x) {
    size_t N=x.size();
    Rcpp::IntegerVector y(N);
    for(size_t i=0; i<N; ++i) y[i]=std::floor(x[i]);
    return y;
  }

  Rcpp::IntegerVector ceiling(Rcpp::NumericVector x) {
    size_t N=x.size();
    Rcpp::IntegerVector y(N);
    for(size_t i=0; i<N; ++i) y[i]=std::ceil(x[i]);
    return y;
  }

  Rcpp::NumericVector sort(Rcpp::NumericVector x) {
   Rcpp::NumericVector y = Rcpp::clone(x);
   std::sort(y.begin(), y.end());
   return y;
  }

  Rcpp::NumericVector quantile(Rcpp::NumericVector x, 
			       Rcpp::NumericVector probs) {
// implementation of type 7
    const int n=x.size(), np=probs.size();
    if(n==0) return x;
    if(np==0) return probs;
/*
    int _np=_probs.size(), np;
    np = _np==0 ? 5 : _np;
    Rcpp::NumericVector probs(np);
    if(_np==0) {
      probs[0]=0.00; probs[1]=0.25; probs[2]=0.50; probs[3]=0.75; probs[5]=1.00;
    }
    Rcpp::IntegerVector lo(Rcpp::sapply(index, floor<double>() )), 
			hi(Rcpp::sapply(index, ceiling<double>() )); 
*/
    Rcpp::NumericVector index=(n-1.)*probs, y=sort(x), x_hi(np), qs(np);
    Rcpp::IntegerVector lo(floor(index)), hi(ceiling(index));

    for(size_t i=0; i<np; ++i) {
      qs[i]=y[lo[i]];
      x_hi[i]=y[hi[i]];
      if((index[i]>lo[i]) && (x_hi[i] != qs[i])) {
	double h;
	h=index[i]-lo[i];
/*
	cout << h << '\n';
	cout << index[i] << '\n';
	cout << lo[i] << '\n';
	cout << x_hi[i] << '\n';
	cout << qs[i] << '\n';
*/
	qs[i]=(1.-h)*qs[i]+h*x_hi[i];
      }
    }

    return qs;
  }

  double Q95(Rcpp::NumericVector x) {
    Rcpp::NumericVector probs(1);
    probs[0]=0.95;
    Rcpp::NumericVector Q95=quantile(x, probs);
    return Q95[0];
  }
}

RcppExport SEXP call_quantile(SEXP _x, SEXP _probs) {
  Rcpp::NumericVector x(_x), probs(_probs);
  return Rcpp::wrap(DPM::quantile(x, probs));
}

RcppExport SEXP call_floor(SEXP _x) {
  Rcpp::NumericVector x(_x);
  return Rcpp::wrap(DPM::floor(x));
}

RcppExport SEXP call_ceiling(SEXP _x) {
  Rcpp::NumericVector x(_x);
  return Rcpp::wrap(DPM::ceiling(x));
}

RcppExport SEXP call_sort(SEXP _x) {
  Rcpp::NumericVector x(_x);
  return Rcpp::wrap(DPM::sort(x));
}

#endif
