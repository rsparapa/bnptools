#ifndef GUARD_rmvnorm
#define GUARD_rmvnorm


// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace Rcpp;

RNGScope scope;

VectorXd rnormXd(int size, double mu=0., double sd=1.);

VectorXd rnormXd(VectorXd mu, MatrixXd Sigma);

//#ifndef NoRcpp

//RcppExport SEXP crnormXd(SEXP, SEXP);

//#endif

//#ifndef NoRcpp
//needs more work
//RcppExport SEXP crnormXd(SEXP n, SEXP mean, SEXP tau, SEXP sd) {
//  arn gen;
//  size_t N = Rcpp::as<int>(n);
//double a=Rcpp::as<double>(mean), b=Rcpp::as<double>(tau), 
//c=Rcpp::as<double>(sd);
//  Rcpp::NumericVector z(N), a(mean), b(tau), c(sd);
//  size_t A=a.size(), B=b.size(), C=c.size();
//  for(size_t i=0; i<N; ++i) z[i]=rtnorm(a[i%A], b[i%B], c[i%C], gen);
//for(size_t i=0; i<N; ++i) z[i]=rtnorm(a, b, c, gen);
//  return Rcpp::wrap(z);
//}

//#endif

VectorXd rnormXd(int size, double mu, double sd) {
  return as<VectorXd>(rnorm(size, mu, sd));
}

VectorXd rnormXd(VectorXd &mu, MatrixXd &Sigma)
{
  MatrixXd L=Sigma.llt().matrixL();

  return L*rnormXd(mu.size())+mu;

}

#endif
