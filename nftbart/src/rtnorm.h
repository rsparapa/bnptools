/*
 * Copyright (C) 2021 Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * rtnorm.h
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

#ifndef GUARD_rtnorm
#define GUARD_rtnorm

double rtnorm(double mean, double tau, double sd, rn& gen);

#ifndef NotInR

RcppExport SEXP crtnorm(SEXP n, SEXP mean, SEXP tau, SEXP sd) {
  rrn gen;
  size_t N = Rcpp::as<int>(n);
//double a=Rcpp::as<double>(mean), b=Rcpp::as<double>(tau), 
//c=Rcpp::as<double>(sd);
  Rcpp::NumericVector z(N), a(mean), b(tau), c(sd);
  size_t A=a.size(), B=b.size(), C=c.size();
  for(size_t i=0; i<N; ++i) z[i]=rtnorm(a[i%A], b[i%B], c[i%C], gen);
  //for(size_t i=0; i<N; ++i) z[i]=rtnorm(a, b, c, gen);
  return Rcpp::wrap(z);
}

#endif

double rtnorm(double mean, double tau, double sd, rn& gen)
{
  double x, z, lambda;

  /* Christian Robert's way */
  //assert(mean < tau); //assert unnecessary: Rodney's way
  tau = (tau - mean)/sd;

  /* originally, the function did not draw this case */
  /* added case: Rodney's way */
  if(tau<=0.) {
    /* draw until we get one in the right half-plane */
    do { z=gen.normal(); } while (z < tau);
  }
  else {
    /* optimal exponential rate parameter */
    lambda = 0.5*(tau + sqrt(tau*tau + 4.0));

    /* do the rejection sampling */
    do {
      z = gen.exp()/lambda + tau;
      //z = lambda*gen.exp() + tau;
    } while (gen.uniform() > exp(-0.5*pow(z - lambda, 2.)));
  }

  /* put x back on the right scale */
  x = z*sd + mean;

  //assert(x > 0); //assert unnecessary: Rodney's way
  return(x);

}

#endif
