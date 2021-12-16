/*
 * Copyright (C) 2021 Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * DPMLIOneal8.h
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

void DPMLIOneal8(SEXP _Y, SEXP _phi, SEXP _C, SEXP _S, 
		 SEXP _prior, SEXP _hyper, rn &eng,
		 double (*F)(double y, double mu, double tau),
		 double (*G0tau)(double a0, double b0, rn &eng),
		 double (*G0mu)(double tau, double m0, double k0, 
				rn &eng),
		 void (*P0)(size_t c_row, SEXP _Y, SEXP _phi, 
			    double m0, double k0, double a0, double b0, 
			    rn &eng)) {

  Rcpp::NumericMatrix Y(_Y);
  //Rcpp::NumericVector Y(_y);
  Rcpp::IntegerVector S(_S), C(_C);
  Rcpp::NumericMatrix phi(_phi);
  //Eigen::MatrixXd phi=Rcpp::as< Eigen::Map<Eigen::MatrixXd> >(_phi);
  Rcpp::List prior(_prior), hyper(_hyper); 

  const int N=Y.nrow(), M=Y.ncol(), p=phi.ncol(), //p=phi.cols(), 
    m=Rcpp::as<int>(prior["m"]);  

  double 
    alpha=Rcpp::as<double>(hyper["alpha"]), 
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

  int k=Rcpp::max(C)+1, // number of states
    k0_draw=Rcpp::as<int>(hyper["k0.draw"]), 
    b0_draw=Rcpp::as<int>(hyper["b0.draw"]), 
    alpha_draw=Rcpp::as<int>(hyper["alpha.draw"]); 

  if(M==3) { // right censoring
    Rcpp::NumericVector prob(k);
    //Eigen::VectorXd prob(k);
    for(size_t i=0; i<N; ++i) if(Y(i, 2)==0.) {
	for(size_t j=0; j<k; ++j) 
	  prob[j]=VEC(S, j)*
	    R::pnorm(Y(i, 1), phi(j, 0), pow(phi(j, 1), -0.5), 0, 0)/N;
	int h;
	h=eng.rcat(prob);
	if(h==-1) h=i;
	Y(i, 0)=eng.rtnorm(Y(i, 1), phi(h, 0), pow(phi(h, 1), -0.5));
      } 
  }

// use LIO non-standardized parameterization
//double Q2, /*_95,*/ scale, scale2;
/*
  Rcpp::NumericVector Y0(Y.column(0));
  Q2=Rcpp::median(Y0);
  //_95=Q95(Y0);
  //scale=0.5*(_95-Q2);
  scale=0.5*(Rcpp::max(Y0)-Q2);
  scale2=pow(scale, 2.);
  m0 += Q2/scale;
  k0_b /= scale2;
  a0 += 0.5;
  b0_b /= scale2;
*/
  
  if(k0_draw) {
    double k0_a_post=k0_a+0.5*k, k0_b_post=k0_b;
    for(size_t i=0; i<k; ++i) 
      k0_b_post += 0.5*phi(i, 1)*pow(phi(i, 0)-m0, 2.);
    k0=eng.gamma(k0_a_post, k0_b_post);
    hyper["k0"]=k0;
  }

  if(b0_draw) {
    double b0_b_post=b0_b;
    for(size_t i=0; i<k; ++i) b0_b_post += phi(i, 1);
    b0=eng.gamma(b0_a+k*a0, b0_b_post);
    hyper["b0"]=b0;
  }

  bool singleton, flag=true;

  for(size_t i=0; i<N; ++i) {
    int c, h, s;
    h=k+m;
    Rcpp::NumericVector prob(h);
    c=C[i];
    s=VEC(S, c);
    singleton=(s==1);
    S[c]=s-1;

    double den, tau, y;
    y=Y(i, 0);
    den=N-1.+alpha;
    for(size_t j=0; j<k; ++j) 
      prob[j]=(VEC(S, j)/den)*((*F) (y, phi(j, 0), phi(j, 1)));
    for(size_t j=k; j<h; ++j) {
      if(singleton && j==k) phi.row(k)=phi.row(c);
      else {
	tau= (*G0tau) (a0, b0, eng);
	phi(j, 0)= (*G0mu) (tau, m0, k0, eng);
	phi(j, 1)=tau;
      }
      prob[j]=(alpha/(m*den))*((*F) (y, phi(j, 0), phi(j, 1)));
    }

    int j;
    j=eng.rcat(prob);

    if(j==-1) j=c; // otherwise C/S will be out of synch, but this is bad
    else if(singleton) {
#ifdef DEBUG
      if(j==c) {
	// this CANNOT happen since S[c]=0!?!
	cout << "j == c:" << c << "\n";
	cout << prob << '\n';
	return;
      }
      else 
#endif
	if(j>=k) {
	  if(j>k) phi.row(c)=phi.row(j);
	  j=c;
	}
	else { // singleton cluster c replaced by c+1 unless k-1
	  if(c<(k-1)) {
	    if(j>c) j--;
	    for(size_t l=c; l<(k-1); ++l) {
	      phi.row(l)=phi.row(l+1);
	      VEC(S, l)=VEC(S, l+1);
	    }
	    S[k-1]=0;
	    for(size_t l=0; l<N; ++l) if(C[l]>c) C[l]=C[l]-1;
	  }
	  k--;
	}
    }
    else if(j>=k) {
      if(j>k) {
	phi.row(k)=phi.row(j);
	j=k;
      }
      k++;
    }

    C[i]=j;
    S[j]=VEC(S, j)+1;

#ifdef DEBUG
    if(flag && (Rcpp::sum(S)!=N || Rcpp::max(C)!=(k-1))) { 
      flag=false;
      cout << "singleton:" << singleton << ' ' 
	   << "k:" << k << ' ' << "c:" << c << ' ' 
	   << "j:" << j << ' ' << "i:" << i << '\n';
      //      COUT << "prob:" << prob.sum() << '\n';
      //      for(size_t j=0; j<=k; ++j) COUT << j << ':' << prob[j] << ' ';
      //      COUT << '\n';
      cout << "S:" << Rcpp::sum(S) << '\n';
      for(size_t j=0; j<N; ++j) if(S[j]>0) cout << j << ':' << S[j] << ' ';
      cout << '\n';
      cout << "C:" << Rcpp::max(C) << '\n';
      for(size_t j=0; j<N; ++j) cout << j << ':' << C[j] << ' ';
      cout << '\n';
    }
#endif
  }

#ifdef DEBUG 
  if(Rcpp::sum(S)!=N || Rcpp::max(C)!=(k-1)) {
    cout << "Sum:" << Rcpp::sum(S) << '\n';
    cout << "Max:" << Rcpp::max(C) << '\n';
    for(size_t i=0; i<N; ++i) { 
      if(S[i]>0)
	cout << "S:" << i << ' ' << S[i] << '\n';
      if(C[i]>=k)
	cout << "C:" << i << ' ' << C[i] << '\n';
    }
  }
#endif

  for(int c=0; c<k; c++) {
    int k_c=VEC(S, c);
#ifdef DEBUG 
    if(k_c<1) {
      cout << "Sum:" << Rcpp::sum(S) << '\n';
      cout << "Max:" << Rcpp::max(C) << '\n';
      for(size_t i=0; i<N; ++i) { 
	if(S[i]<1)
	  cout << "S:" << i << ' ' << S[i] << '\n';
	if(C[i]>=k)
	  cout << "C:" << i << ' ' << C[i] << '\n';
      }
    }
    else {
#endif
	Rcpp::NumericVector y(k_c);
	//Eigen::VectorXd y(k_c);
	//Eigen::MatrixXd A(k_c, M);
  
	for(int i=0, j=0; i<N && j<k_c; ++i) if(C[i] == c) {
	    y[j]=Y(i, 0);
	    //A.row(j)=Y.row(i);
	    j++;
	  }
    
	(*P0) (c, Rcpp::wrap(y), _phi, m0, k0, a0, b0, eng);
#ifdef DEBUG 
    }
#endif
  }

  if(alpha_draw) {
    double eta=eng.beta(alpha, N);
    alpha=eng.gamma(alpha_a+k, alpha_b-log(eta)); 
    hyper["alpha"]=alpha;
  }
}
