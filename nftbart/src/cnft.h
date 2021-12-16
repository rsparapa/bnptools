/*
 * Copyright (C) 2012-2021 Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * cnft.h
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

/*
#define ARNG
#include "rn.h"

#include <iostream>
#include <string>
#include <ctime>
#include <sstream>

#include <fstream>
#include <vector>
#include <limits>

#include "ambrt.h"
#include "psbrt.h"
#include "brt.h"
#include "brtfuns.h"
#include "dinfo.h"
#include "mbrt.h"
#include "treefuns.h"
#include "tree.h"
*/

//#include "rrn.h"
//#include "rtnorm.h"

RcppExport SEXP cnft(
		     SEXP _ix, 		//train x
		     SEXP _iy, 		//train y
		     SEXP _idelta,         //delta
		     //  SEXP _ievents,        //events
		     SEXP _ixp,		//test x
		     SEXP _im, 		//num trees mean and var
		     SEXP _ind,		//number of draws kept  (adapt, then burn, then draw)
		     SEXP _iburn, 	        //number of burn draws
		     SEXP _inadapt, 	//number of adapt draws 
		     SEXP _iadaptevery,	//adapt every adaptevery draws
		     SEXP _itau,		//mu_i ~N(0,\tau^2)
		     SEXP _ioveralllambda, //sigma_i ~ prod (nu'*lambda')/chisquared^_nu', nu',lambda' from overallnu, overalllambda
		     SEXP _ioverallnu,
		     SEXP _ibase,		//tree prior for mean and var
		     SEXP _ipower,	        //tree prior for mean and var
		     SEXP _itc,		//thread count
		     SEXP _isigmav,	//intial values for sigma_i 
		     SEXP _ichv,		//change variable (correlation) matrix
		     SEXP _ipbd,		//prob of birth/death step, mean and var trees
		     SEXP _ipb,		//pb prob of birth given birth/death
		     SEXP _istepwpert,	//(in (0,1)) shrinkage factor for range of acceptable cutpoints for MH, mean and var, (adapts)
		     SEXP _iprobchv,	//prob of trying change variable move, mean and var, (prob of perturb is 1 minus this)
		     SEXP _iminnumbot,	//minimum number of observations in a bottom node, mean and var
		     SEXP _iprintevery,    //how often you print progress
		     SEXP _ixicuts,        //variable cutpoints for each predictor
		     //SEXP _isummarystats,  //boolean, do you want summary stats
//		     SEXP _ialphao, 
//		     SEXP _ibetao,
		     //  SEXP _mstart, //  SEXP _sstart,
		     SEXP _hyper, 
		     SEXP _C,
		     SEXP _states, 
		     SEXP _phi,
		     SEXP _prior,
		     //SEXP _idraws,
		     SEXP _idrawMuTau,
		     SEXP _impute_bin, 
		     SEXP _impute_prior 
		     // SEXP _impute_mult, // integer vector of column indicators for missing covariates
		     // SEXP _impute_miss, // integer vector of row indicators for missing values
		     // SEXP _impute_prior // matrix of prior missing imputation probability
		     )
{
  //random number generation

#if defined(RRNG)
  GetRNGstate();
#endif

  rrn gen;

  //--------------------------------------------------
  //process args
  //x
  Rcpp::NumericMatrix Xt(_ix);
  size_t p = Xt.nrow();
  size_t n = Xt.ncol();
  double *x = &Xt[0];

  int impute_bin = Rcpp::as<int>(_impute_bin), impute_flag=(impute_bin>=0);
  Rcpp::NumericVector impute_prior(_impute_prior); 
  double impute_post0, impute_post1, *impute_Xrow_ptr = 0;
  /*
    Rcpp::IntegerVector impute_mult(_impute_mult); // integer vector of column indicators for missing covariates
    size_t _K = impute_mult.size(); // number of columns to impute
    Rcpp::IntegerVector impute_miss(_impute_miss); // length n: integer vector of row indicators for missing values
    Rcpp::NumericMatrix impute_prior(_impute_prior); // n X K: matrix of prior missing imputation probability
    Rcpp::NumericVector impute_post(_K); // length K: double vector of posterior missing imputation probability
    Rcpp::NumericVector impute_fhat(_K); 
    double *impute_fhat_ptr = 0, *impute_shat_ptr = 0, *impute_Xrow_ptr = 0;
    if(_K>0) impute_fhat_ptr = &impute_fhat[0];
    Rcpp::NumericVector impute_shat(_K); 
  */   
  //y
  Rcpp::NumericVector yv(_iy);
  double *y = &yv[0], ymax=Rcpp::max(yv);

  //delta
  Rcpp::IntegerVector delta(_idelta), censor(n);
/*
  Rcpp::NumericVector deltav(_idelta);
  double *delta = &deltav[0];
*/

  //events
  //  Rcpp::NumericVector events(_ievents);
  //size_t K=events.size();

  //z and w
  Rcpp::NumericVector zv(n), wv(n);
  double *z = &zv[0], *w = &wv[0];

  //x out-of-sample
  Rcpp::NumericMatrix xpm(_ixp);
  size_t np = xpm.ncol();
  double *xp = nullptr;
  if(np)  xp = &xpm[0];

  //number of trees
  Rcpp::IntegerVector im(_im);
  size_t m = im[0], mh = im[1];

  //mu prior (tau, ambrt) and sigma prior (lambda,nu, psbrt)
  double tau = Rcpp::as<double>(_itau);
  double overalllambda = Rcpp::as<double>(_ioveralllambda);
  double overallnu = Rcpp::as<double>(_ioverallnu);

  //nd and burn
  size_t nd = Rcpp::as<int>(_ind);
  size_t burn = Rcpp::as<int>(_iburn);
  size_t nadapt = Rcpp::as<int>(_inadapt);
  size_t adaptevery = Rcpp::as<int>(_iadaptevery);

  //tree strings
  /*
    std::stringstream mtrees, strees;  
    mtrees.precision(10);
    mtrees << nd << " " << m << " " << p << endl;
    strees.precision(10);
    strees << nd << " " << mh << " " << p << endl;
  */

  int drawMuTau = Rcpp::as<int>(_idrawMuTau);
    //draws = Rcpp::as<int>(_idraws), 
    //drawsd=0; //(1-draws)*(1-drawMuTau);

  Rcpp::NumericMatrix mdraws(nd, n), //sdraws(nd, pow(n, draws)),
    sdraws(nd, n), mpred(nd, np), spred(nd, np), zdraws(nd, n);
  //Rcpp::NumericVector sddraws((nd+burn)*drawsd);

  // for varcounts
  Rcpp::IntegerMatrix fvc(nd, p), svc(nd, p);
    Rcpp::IntegerVector varcount(p);
    unsigned int tmaxd=0;
    unsigned int tmind=0;
    double tavgd=0.0;

  //tree prior
  Rcpp::NumericVector ibase(_ibase), ipower(_ipower);
  double falpha=ibase[0], salpha=ibase[1], 
    mybeta=ipower[0], mybetah=ipower[1];

  //thread count
  int tc = Rcpp::as<int>(_itc);

  //sigma vector
  Rcpp::NumericVector sigmav(_isigmav);
  double *sig = &sigmav[0];
  dinfo disig;
  disig.n=n; disig.p=p; disig.x = x; disig.y = sig; disig.tc=tc;

  //f(x) function to be used like sig vector
  double *fun = new double[n]; 

  //change variable
  Rcpp::NumericMatrix chvm(_ichv);
  COUT << "row, cols chvm: " << chvm.nrow() << ", " << chvm.ncol() << endl;
  std::vector<std::vector<double> > chgv(chvm.nrow());
  for(int i=0;i<chvm.nrow();i++) {
    chgv[i].resize(chvm.ncol());
    for(int j=0;j<chvm.ncol();j++) chgv[i][j]=chvm(i,j);
  }

  //control
  size_t printevery = Rcpp::as<int>(_iprintevery);
  Rcpp::NumericVector ipbd(_ipbd), ipb(_ipb), istepwpert(_istepwpert),
    iprobchv(_iprobchv);
  double pbd=ipbd[0], pbdh=ipbd[1], pb=ipb[0], pbh=ipb[1],
    stepwpert = istepwpert[0], stepwperth = istepwpert[1],
    probchv = iprobchv[0], probchvh = iprobchv[1];
  Rcpp::IntegerVector iminnumbot(_iminnumbot);
  size_t minnumbot = iminnumbot[0], minnumboth = iminnumbot[1];
  bool dopert=true, doperth=true;
  if(probchv<0) dopert=false;
  if(probchvh<0) doperth=false;

  //summary statistics yes/no
  bool summarystats = true; //Rcpp::as<bool>(_isummarystats);

  //error variance prior
/*double alphao = Rcpp::as<double>(_ialphao);
  double betao = Rcpp::as<double>(_ibetao);
  double lambdao=betao/alphao, nuo=2.*alphao;*/

  // DPM LIO
  Rcpp::IntegerVector C(_C), states(_states);
  Rcpp::NumericMatrix phi(_phi);
  Rcpp::List prior(_prior), hyper(_hyper);
  const int neal_m=Rcpp::as<int>(prior["m"]), 
    constrain=Rcpp::as<int>(prior["constrain"]);

  double alpha=Rcpp::as<double>(hyper["alpha"]); // inital value
  Rcpp::NumericMatrix Y(n, 1);

  // return data structures
  int ndMT=0, ndbMT=0, nMT=0;
  if(drawMuTau>0) {
    ndMT=nd;
    ndbMT=(nd+burn);
    nMT=n;
  }
  Rcpp::NumericVector dalpha(ndMT); 
  for(int i=0;i<ndMT;i++) dalpha[i]=0.;
  Rcpp::IntegerVector dnpart(ndbMT);
  Rcpp::NumericMatrix dmu(ndMT, nMT);
  Rcpp::NumericMatrix dsig(ndMT, nMT);
  Rcpp::IntegerMatrix dpC(ndMT, nMT);
  Rcpp::NumericMatrix dpMU(ndMT, nMT);
  Rcpp::NumericMatrix dpSD(ndMT, nMT);
  //Rcpp::NumericMatrix dpWT(ndMT, nMT);

  //--------------------------------------------------
  //print args
  //  Rprintf("**********************\n");
  Rprintf("n: %ld\n",n);
  Rprintf("p: %ld\n",p);
  Rprintf("first row: %lf, %lf\n",x[0],x[p-1]);
  Rprintf("second row: %lf, %lf\n",x[p],x[2*p-1]);
  Rprintf("last row: %lf, %lf\n",x[(n-1)*p],x[n*p-1]);
  Rprintf("first and last y: %lf, %lf\n",y[0],y[n-1]);
  if(np) {
    Rprintf("np: %d\n",np);
    Rprintf("first row xp: %lf, %lf\n",xp[0],xp[p-1]);
    Rprintf("second row xp: %lf, %lf\n",xp[p],xp[2*p-1]);
    Rprintf("last row xp : %lf, %lf\n",xp[(np-1)*p],xp[np*p-1]);
  } else {
    Rprintf("no test observations\n");
  }
  Rprintf("number of trees mean: %ld\n",m);
  Rprintf("number of trees stan dev: %ld\n",mh);
  Rprintf("tau: %lf\n",tau);
  Rprintf("overalllambda: %lf\n",overalllambda);
  Rprintf("overallnu: %lf\n",overallnu);
  Rprintf("burn (nskip): %ld\n",burn);
  Rprintf("nd (ndpost): %ld\n",nd);
  Rprintf("nadapt: %ld\n",nadapt);
  Rprintf("adaptevery: %ld\n",adaptevery);
  Rprintf("mean tree prior base: %lf\n",falpha);
  Rprintf("mean tree prior power: %lf\n",mybeta);
  Rprintf("variance tree prior base: %lf\n",salpha);
  Rprintf("variance tree prior power: %lf\n",mybetah);
  Rprintf("thread count: %ld\n",tc);
  Rprintf("first and last sigmav: %lf, %lf\n",sigmav[0],sigmav[n-1]);
  Rprintf("chgv first row: %lf, %lf\n",chgv[0][0],chgv[0][p-1]);
  Rprintf("chgv last row: %lf, %lf\n",chgv[p-1][0],chgv[p-1][p-1]);
  Rprintf("mean trees prob birth/death: %lf\n",pbd);
  Rprintf("mean trees prob birth: %lf\n",pb);
  Rprintf("variance trees prob birth/death: %lf\n",pbdh);
  Rprintf("variance trees prob birth: %lf\n",pbh);
  Rprintf("mean trees initial step width pert move: %lf\n",stepwpert);
  Rprintf("variance trees initial step width pert move: %lf\n",stepwperth);
  Rprintf("mean trees prob of a change var move : %lf\n",probchv);
  Rprintf("variance trees prob of a change var move : %lf\n",probchvh);
  Rprintf("mean trees min num obs in bottom node: %ld\n",minnumbot);
  Rprintf("variance trees min num obs in bottom node: %ld\n",minnumboth);
//  Rprintf("sigma prior : alphao=%lf, betao=%lf\n", alphao, betao);
  //Rprintf("base prior : muo=%lf, ko=%lf, alphao=%lf, betao=%lf\n",
  // muo, ko, alphao, betao);
  //Rprintf("init values : mu=%lf, tau=%lf, sd=%lf, alpha=%lf\n",
  // mstart, tauinit, sstart, alpha);
  COUT << "m:" << Rcpp::as<int>(prior["m"]) << '\n';
  COUT << "states:\n"; for(size_t u=0; u<R::imin2(5, states.size()); ++u) 
			 COUT << states[u] << ' ';
  COUT << '\n';
  COUT << "C:\n"; for(size_t u=0; u<R::imin2(5, C.size()); ++u) 
		    COUT << C[u] << ' ';
  COUT << '\n';
  COUT << "phi:\n"; for(size_t u=0; u<R::imin2(5, phi.rows()); ++u) {
    for(size_t v=0; v<phi.cols(); ++v)
      COUT << phi(u, v) << ' ';
    COUT << '\n'; }
  Rprintf("alpha draw: alpha_a=%lf, alpha_b=%lf\n", 
	  Rcpp::as<double>(prior["alpha.a"]),
	  Rcpp::as<double>(prior["alpha.b"]));
  Rprintf("draw: MuTau=%ld\n", drawMuTau);
  //Rprintf("draw: s=%ld, MuTau=%ld\n", draws, drawMuTau);
  if(impute_flag) {
    cout << "Missing imputation column index:\n" << impute_bin << endl;
    cout << "Missing imputation probabilities:\n index 0=" << impute_prior[0]
	 << ',' << "index n-1=" << impute_prior[n-1] << endl;
  }
  Rprintf("printevery: %d\n",printevery);

  //--------------------------------------------------
  //make xinfo
  xinfo xi;
  Rcpp::List ixi(_ixicuts);
  if(ixi.length()!= (int)p) {
    Rprintf("Cutpoint definition does not match number of predictor variables!\n");
    return 0;
  }
  xi.resize(p);

  for(size_t i=0;i<p;i++)
    xi[i]=Rcpp::as< std::vector<double> >(ixi[i]);

  //summarize input variables:
  for(size_t i=0;i<p;i++)
    {
      COUT << "Variable " << i << " has numcuts=" << xi[i].size() << " : ";
      COUT << xi[i][0] << " ... " << xi[i][xi[i].size()-1] << endl;
    }

  for(size_t i=0;i<n;i++) {
    z[i]=y[i]; // initialize z
    w[i]=1.;   // initialize w
    censor[i] = 1-delta[i]; // -1 left, 0 event, 1 right
    Y(i, 0)=y[i];
    if(impute_flag) 
      Xt(impute_bin, i)=gen.bin(1, impute_prior[i]);
    /*
      if(_K>0) {
      if(impute_miss[i]==1) {
      for(size_t j=0; j<_K; j++) {
      XV(impute_mult[j], i)=0;
      }
      size_t k;
      Rcpp::NumericVector impute_prob(impute_prior.row(i));
      k=gen.rcat(impute_prob); // use prior prob only
      XV(impute_mult[k], i)=1;
      }
      }
    */
  }

  //--------------------------------------------------
  //dinfo
  dinfo di;
  di.n=n; di.p=p; di.x = x; di.y = z; di.tc=tc;
  //--------------------------------------------------
  // set up ambrt object
  ambrt ambm(m);

  //cutpoints
  ambm.setxi(&xi);    //set the cutpoints for this model object
  //data objects
  ambm.setdata(&di);  //set the data
  //thread count
  ambm.settc(tc);      //set the number of threads when using OpenMP, etc.
  //tree prior
  ambm.settp(falpha, //the alpha parameter in the tree depth penalty prior
	     mybeta     //the beta parameter in the tree depth penalty prior
	     );
  //MCMC info
  ambm.setmi(
	     pbd,  //probability of birth/death
	     pb,  //probability of birth
	     minnumbot,    //minimum number of observations in a bottom node
	     dopert, //do perturb/change variable proposal?
	     stepwpert,  //initialize stepwidth for perturb proposal.  If no adaptation it is always this.
	     probchv,  //probability of doing a change of variable proposal.  perturb prob=1-this.
	     &chgv  //initialize the change of variable correlation matrix.
	     );
  ambm.setci(tau,sig);


  //--------------------------------------------------
  //setup psbrt object
  psbrt psbm(mh,overalllambda);

  //make di for psbrt object
  dinfo dips;
  double *r = new double[n]; 
  for(size_t i=0;i<n;i++) r[i]=sigmav[i];
  dips.x = x; dips.y=r; dips.p=p; dips.n=n; dips.tc=tc;

  double opm=1.0/((double)mh);
  double nu=2.0*pow(overallnu,opm)/(pow(overallnu,opm)-pow(overallnu-2.0,opm));
  double lambda=pow(overalllambda,opm);

  //cutpoints
  psbm.setxi(&xi);    //set the cutpoints for this model object
  //data objects
  psbm.setdata(&dips);  //set the data
  //thread count
  psbm.settc(tc); 
  //tree prior
  psbm.settp(salpha, //the alpha parameter in the tree depth penalty prior
	     mybetah     //the beta parameter in the tree depth penalty prior
	     );
  psbm.setmi(
	     pbdh,  //probability of birth/death
	     pbh,  //probability of birth
	     minnumboth,    //minimum number of observations in a bottom node
	     doperth, //do perturb/change variable proposal?
	     stepwperth,  //initialize stepwidth for perturb proposal.  If no adaptation it is always this.
	     probchvh,  //probability of doing a change of variable proposal.  perturb prob=1-this.
	     &chgv  //initialize the change of variable correlation matrix.
	     );
  psbm.setci(nu,lambda);

  // x.test predictions
  double *_f = new double[np], *_s = new double[np];
  dinfo f, s;
  f.x = xp; f.y=_f; f.p = p; f.n=np; f.tc=tc;
  s.x = xp; s.y=_s; s.p = p; s.n=np; s.tc=tc;

  // DPM LIO
  double *mvec = new double[n];
  double *svec = new double[n];
  //Rcpp::NumericVector mvec(n), svec(n);
  for(size_t i=0;i<n;i++) {
    mvec[i]=phi(0, 0);
    svec[i]=pow(phi(0, 1), -0.5);
  }
  //--------------------------------------------------
  //run mcmc
  std::vector<int> onn(nd*m,1);
  std::vector<std::vector<int> > oid(nd*m, std::vector<int>(1));
  std::vector<std::vector<int> > ovar(nd*m, std::vector<int>(1));
  std::vector<std::vector<int> > oc(nd*m, std::vector<int>(1));
  std::vector<std::vector<double> > otheta(nd*m, std::vector<double>(1));
  std::vector<int> snn(nd*mh,1);
  std::vector<std::vector<int> > sid(nd*mh, std::vector<int>(1));
  std::vector<std::vector<int> > svar(nd*mh, std::vector<int>(1));
  std::vector<std::vector<int> > sc(nd*mh, std::vector<int>(1));
  std::vector<std::vector<double> > stheta(nd*mh, std::vector<double>(1));
  //brtMethodWrapper fambm(&brt::f,ambm);
  //brtMethodWrapper fpsbm(&brt::f,psbm);

  // x predictions
  double fhat0, fhat1, shat0=1., shat1=1.;
  dinfo f0, f1, s0, s1;
  f0.y=&fhat0; f0.p = p; f0.n=1; f0.tc=tc;
  f1.y=&fhat1; f1.p = p; f1.n=1; f1.tc=tc;
  s0.y=&shat0; s0.p = p; s0.n=1; s0.tc=tc;
  s1.y=&shat1; s1.p = p; s1.n=1; s1.tc=tc;

  Rprintf("Starting MCMC...\n");

  bool adapting=true, burning=false, keeping=false, adapting_every, drawDP=false;
//, draw_s=true; 
#ifdef PROFILER 
  ProfilerStart(PROFILER);
#endif

  for(size_t h, i=0, j, L=nadapt+burn, M=L+nd; i<M; i++) {
    adapting=(i<nadapt);
    h=i-nadapt;
    j=i-L;

    if(!adapting) {
      adapting_every=false;
      burning=(i<L); 
      if(!burning) keeping=true;
      if(i==nadapt) {
	drawDP=(drawMuTau>0);
      }

      if(burning && (h % printevery)==0) COUT << "draw burn " << h << endl;
      else if(keeping) {
	if ((j % printevery)==0) COUT << "draw keep " << j << endl;
	//if(summarystats && j==0) {
	if(j==0) {
	  ambm.setstats(true);
	  psbm.setstats(true);
	}
      }
    } 
    else adapting_every=(i>0 && (i%adaptevery)==0);

    if(adapting_every) {
      COUT << "adapting  " << i << endl;
    }

    size_t K=Rcpp::max(C)+1;
    for(size_t k=0;k<n;k++) {
      if(delta[k]==1) z[k] = y[k];
      else {
	if(drawDP) {
	  Rcpp::NumericVector prob(K);
	  for(size_t g=0;g<K;g++) 
	    prob[g]=states[g]*
	      R::pnorm(y[k], phi(g, 0)*sig[k]+ambm.f(k), 
		       pow(phi(g, 1), -0.5)*sig[k], 
		       (int)(delta[k]==2), 0)/n; // left and right
		       //pow(phi(g, 1), -0.5)*sig[k], 0, 0)/n; right only
	  size_t g;
	  g=gen.rcat(prob);
	  if(g==-1) g=C[k];
	  // right and left censoring
	  z[k]=censor[k]*gen.rtnorm(censor[k]*y[k], 
				    censor[k]*(phi(g, 0)*sig[k]+ambm.f(k)), 
			  pow(phi(g, 1), -0.5)*sig[k]);
	}
	else z[k]=censor[k]*gen.rtnorm(censor[k]*y[k], 
				       censor[k]*ambm.f(k), sig[k]);
/* right censoring only
	  z[k]=gen.rtnorm(y[k], phi(g, 0)*sig[k]+ambm.f(k), 
			  pow(phi(g, 1), -0.5)*sig[k]);
	}
	else z[k]=gen.rtnorm(y[k], ambm.f(k), sig[k]);
*/
      }

      if(keeping) zdraws(j, k)=z[k];
    }

    for(size_t k=0;k<n;k++) {
      if(impute_flag) {
	if(impute_prior[k]>0. && impute_prior[k]<1.) {
	  impute_Xrow_ptr=&Xt(0, k);
	  impute_post0=1.-impute_prior[k];
	  impute_post1=impute_prior[k];
	  Xt(impute_bin, k)=0;
	  f0.x=impute_Xrow_ptr;
	  ambm.predict(&f0); // result in fhat0
/*
	  if(draws) {
	    s0.x=impute_Xrow_ptr;
	    psbm.predict(&s0); // result in shat0
	  }
*/
	  Xt(impute_bin, k)=1;
	  f1.x=impute_Xrow_ptr;
	  ambm.predict(&f1); // result in fhat1
// imputation too simplistic if drawDP: need mvec/svec
	  // if(draws) {
	  //   s1.x=impute_Xrow_ptr;
	  //   psbm.predict(&s1); // result in shat1
	  //   impute_post0 *= R::dnorm(z[k], fhat0, shat0, 0);
	  //   impute_post1 *= R::dnorm(z[k], fhat1, shat1, 0);
	  // }
	  // else
	    {
	    impute_post0 *= R::dnorm(z[k], fhat0, sig[k], 0);
	    impute_post1 *= R::dnorm(z[k], fhat1, sig[k], 0);
	  }
	  Xt(impute_bin, k)=
	    gen.bin(1, impute_post1/(impute_post0+impute_post1));
	}
      }

      if(drawDP) {
	z[k] = z[k]-mvec[k]*sig[k];
	sig[k] = sig[k]*svec[k];
	//if(draws) sig[k] = sig[k]*svec[k];
	//else sig[k] = svec[k];
      }
    }

    ambm.draw(gen);
    for(size_t k=0;k<n;k++) fun[k]=ambm.f(k);
    if(adapting_every) ambm.adapt();

    //if(draws) 
    {
      for(size_t k=0;k<n;k++) { 
	r[k]=z[k]-fun[k]; 
	//r[k]=z[k]-ambm.f(k); 
	if(drawDP) r[k]=r[k]/svec[k];
      }
      //dips = di;    
      //dips -= fambm;
      psbm.draw(gen);

      for(size_t k=0;k<n;k++) { 
	sig[k] = psbm.f(k);
	//if(drawDP) sig[k] *= svec[k]; 
      }
      //disig = fpsbm;
      if(adapting_every) psbm.adapt(); 
    }
/*
    else if(!drawDP) {
      double rss, sigma;
      rss=nuo*lambdao;
      for(size_t k=0;k<n;k++) {
	//rss += pow(z[k]-ambm.f(k), 2.);
	rss += pow(z[k]-fun[k], 2.);
      }
      sigma = sqrt(rss/gen.chi_square(n+nuo));
      for(size_t k=0;k<n;k++) sig[k]=sigma;
      sddraws[h] = sigma; 
    }
*/

    if(drawDP) {
      for(size_t k=0;k<n;k++) {
	//z[k] = mvec[k]+(z[k]-ambm.f(k))/sig[k];
	z[k] = mvec[k]+(z[k]-fun[k])/sig[k];
	Y(k, 0) = z[k];
	//Y[k] = z[k];
	//Y[k] = mvec[k]+(z[k]-ambm.f(k))/sig[k];
	//sig[k] /= svec[k];
	//Y[k] = (z[k] + mvec[k]-ambm.f(k))/sig[k];
	//w[k] = 1./sig[k];
      }

      //if(drawMuTau==1)
	DPMLIOneal8(Rcpp::wrap(Y), _phi, _C, _states, _prior, _hyper,
		  gen, &DPMLIOmutau_F, &DPMLIOmutau_G0tau, 
		  &DPMLIOmutau_G0mu, &DPMLIOmutau_P0);
      // else
      // 	DPMtauneal8(Rcpp::wrap(Y), _phi, _C, _states, _prior, _hyper,
      // 		  gen, &DPMtau_F, &DPMtau_G0, &DPMtau_P0);

      // testing mu.=0
      //mvec[k]=0.;
      double mu0, var0; 
      mu0=0.; var0=0.; 
      for(size_t k=0;k<n;k++) {
	mvec[k] = phi(C[k], 0);
	mu0 += mvec[k]/n;
	svec[k] = pow(phi(C[k], 1), -0.5);
	var0 += pow(svec[k], 2.)/n;
      }

      if(constrain) {
	double sd0=sqrt(var0); 
	for(size_t k=0;k<n;k++) {
	  mvec[k] = mvec[k]-mu0;
	  phi(C[k], 0)=mvec[k];
	  //if(draws) {
	    svec[k] = svec[k]/sd0;
	    phi(C[k], 1)=pow(svec[k], -2.);
	  //}
	}
      }

      dnpart[h]=Rcpp::max(C)+1;
    }

    /*
      if(impute_flag) {
      for(size_t k=0; k<n; ++k) {
      if(impute_prior[k]>0. && impute_prior[k]<1.) {
      impute_Xrow_ptr=&Xt(0, k);
      impute_post0=1.-impute_prior[k];
      impute_post1=impute_prior[k];
      Xt(impute_bin, k)=0;
      f0.x=impute_Xrow_ptr;
      ambm.predict(&f0); // result in fhat0
      if(draws) {
      s0.x=impute_Xrow_ptr;
      psbm.predict(&s0); // result in shat0
      }
      Xt(impute_bin, k)=1;
      f1.x=impute_Xrow_ptr;
      ambm.predict(&f1); // result in fhat1
      if(draws) {
      s1.x=impute_Xrow_ptr;
      psbm.predict(&s1); // result in shat1
      if(drawDP) {
      impute_post0 *= 
      R::dnorm(z[k], mvec[k]+fhat0, svec[k]*shat0, 0);
      impute_post1 *= 
      R::dnorm(z[k], mvec[k]+fhat1, svec[k]*shat1, 0);
      } 
      else {
      impute_post0 *= R::dnorm(z[k], fhat0, shat0, 0);
      impute_post1 *= R::dnorm(z[k], fhat1, shat1, 0);
      }
      }
      else {
      if(drawDP) {
      impute_post0 *= 
      R::dnorm(z[k], mvec[k]+fhat0, svec[k]*sig[k], 0);
      impute_post1 *= 
      R::dnorm(z[k], mvec[k]+fhat1, svec[k]*sig[k], 0);
      } else {
      impute_post0 *= R::dnorm(z[k], fhat0, sig[k], 0);
      impute_post1 *= R::dnorm(z[k], fhat1, sig[k], 0);
      }
      }
      Xt(impute_bin, k)=
      gen.bin(1, impute_post0/(impute_post0+impute_post1));
      }
      }
      }
    */
  
    /* cannot use j or h in this loop
       if(_K>0) {
       for(size_t k=0; k<n; ++k) {
       if(impute_miss[k]==1) {
       impute_Xrow_ptr=&XV(0, k);
       impute_post=impute_prior.row(k);
       for(size_t j=0; j<_K; ++j) {
       for(size_t h=0; h<_K; ++h) XV(impute_mult[h], k)=0;
       XV(impute_mult[j], k)=1;
       _X.x=impute_Xrow_ptr;
       ambm.predict(&_X);
       impute_fhat_ptr[j]=_x;
       //ambm.predict(p, 1, impute_Xrow_ptr, &impute_fhat_ptr[j]);
       if(draws) {
       psbm.predict(&_X);
       impute_shat_ptr[j]=_x;
       //psbm.predict(p, 1, impute_Xrow_ptr, &impute_shat_ptr[j]);
       if(drawDP) {
       impute_post[j] *= R::dnorm(z[k], mvec[k]+impute_fhat_ptr[j],
       svec[k]*impute_shat_ptr[j], 0);
       } else {
       impute_post[j] *= R::dnorm(z[k], impute_fhat_ptr[j],
       impute_shat_ptr[j], 0);
       }
       }
       else {
       if(drawDP) {
       impute_post[j] *= R::dnorm(z[k], mvec[k]+impute_fhat_ptr[j],
       svec[k]*sig[k], 0);
       } else {
       impute_post[j] *= R::dnorm(z[k], impute_fhat_ptr[j],
       sig[k], 0);
       }
       }
       }
       for(size_t j=0; j<_K; j++) { 
       XV(impute_mult[j], k)=0;
       //prevXV[impute_mult[j]]=0;
       }
       size_t h;
       h=gen.rcat(impute_post); 
       XV(impute_mult[h], k)=1;
       //prevXV[impute_mult[h]]=1;
       }
       }
       }
    */

    if(keeping) {
      //save tree to vec format

      ambm.savetree(j,m,onn,oid,ovar,oc,otheta);
      //if(draws) 
      psbm.savetree(j,mh,snn,sid,svar,sc,stheta);

      if(drawMuTau>0) dalpha[j]=hyper["alpha"];

      ambm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
      for(size_t k=0;k<p;k++) fvc(j, k)=varcount[k];

      //if(draws) 
      {
	psbm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
	for(size_t k=0;k<p;k++) svc(j, k)=varcount[k];
      }

      double wt=1./n;
      for(size_t k=0;k<n;k++) {
	//mdraws(j, k) = ambm.f(k);
	mdraws(j, k) = fun[k];
	//if(draws) {
	  //sdraws(j, k) = psbm.f(k);
	  sdraws(j, k) = sig[k];
	//}
	if(drawMuTau>0) {
	  dmu(j, k) = mvec[k];
	  dsig(j, k) = svec[k];
	  dpMU(j, k)=phi(k, 0);
	  dpSD(j, k)=pow(phi(k, 1), -0.5);
	  dpC(j, k) = C[k]+1;
	  //dpWT(j, C[k])=dpWT(j, C[k])+wt;
	}
      }

      if(np) {
	ambm.predict(&f);
	//if(draws) 
	psbm.predict(&s);

	for(size_t k=0;k<np;k++) {
	  mpred(j, k) = f.y[k];
	  //if(draws) 
	  spred(j, k) = s.y[k];
	}
      }
    }
  }

#ifdef PROFILER
  ProfilerStop();
#endif

  std::stringstream mtrees, strees;  
  mtrees=ambm.gettrees(nd,m,onn,oid,ovar,oc,otheta,0.);
  //if(draws) 
  strees=psbm.gettrees(nd,mh,snn,sid,svar,sc,stheta,-1.);

  /*
    Flatten posterior trees to a few XXL vectors so we can just pass pointers
    to these vectors back to R (which is much faster than copying all the data back)
  */

  std::vector<int>* e_ots=new std::vector<int>(nd*m);
  std::vector<int>* e_oid=new std::vector<int>;
  std::vector<int>* e_ovar=new std::vector<int>;
  std::vector<int>* e_oc=new std::vector<int>;
  std::vector<double>* e_otheta=new std::vector<double>;
  for(size_t i=0;i<nd;i++)
    for(size_t j=0;j<m;j++) {
      e_ots->at(i*m+j)=static_cast<int>(oid[i*m+j].size());
      e_oid->insert(e_oid->end(),oid[i*m+j].begin(),oid[i*m+j].end());
      e_ovar->insert(e_ovar->end(),ovar[i*m+j].begin(),ovar[i*m+j].end());
      e_oc->insert(e_oc->end(),oc[i*m+j].begin(),oc[i*m+j].end());
      e_otheta->insert(e_otheta->end(),otheta[i*m+j].begin(),otheta[i*m+j].end());
    }
  Rcpp::XPtr< std::vector<int> > extern_ots(e_ots,true);
  Rcpp::XPtr< std::vector<int> > extern_oid(e_oid,true);
  Rcpp::XPtr< std::vector<int> > extern_ovar(e_ovar,true);
  Rcpp::XPtr< std::vector<int> > extern_oc(e_oc,true);
  Rcpp::XPtr< std::vector<double> > extern_otheta(e_otheta,true);

  Rcpp::List ret = Rcpp::List::create(
				      Rcpp::Named("ots")=extern_ots,
				      Rcpp::Named("oid")=extern_oid,
				      Rcpp::Named("ovar")=extern_ovar,
				      Rcpp::Named("oc")=extern_oc,
				      Rcpp::Named("otheta")=extern_otheta,
				      /*
					Rcpp::Named("sts")=extern_sts,
					Rcpp::Named("sid")=extern_sid,
					Rcpp::Named("svar")=extern_svar,
					Rcpp::Named("sc")=extern_sc,
					Rcpp::Named("stheta")=extern_stheta,
				      */
				      Rcpp::Named("f.train")=mdraws,
				      Rcpp::Named("f.trees")=Rcpp::CharacterVector(mtrees.str()));

  //if(draws) 
  {
    std::vector<int>* e_sts=new std::vector<int>(nd*mh);
    std::vector<int>* e_sid=new std::vector<int>;
    std::vector<int>* e_svar=new std::vector<int>;
    std::vector<int>* e_sc=new std::vector<int>;
    std::vector<double>* e_stheta=new std::vector<double>;
    for(size_t i=0;i<nd;i++)
      for(size_t j=0;j<mh;j++) {
	e_sts->at(i*mh+j)=static_cast<int>(sid[i*mh+j].size());
	e_sid->insert(e_sid->end(),sid[i*mh+j].begin(),sid[i*mh+j].end());
	e_svar->insert(e_svar->end(),svar[i*mh+j].begin(),svar[i*mh+j].end());
	e_sc->insert(e_sc->end(),sc[i*mh+j].begin(),sc[i*mh+j].end());
	e_stheta->insert(e_stheta->end(),stheta[i*mh+j].begin(),stheta[i*mh+j].end());
      }

    Rcpp::XPtr< std::vector<int> > extern_sts(e_sts,true);
    Rcpp::XPtr< std::vector<int> > extern_sid(e_sid,true);
    Rcpp::XPtr< std::vector<int> > extern_svar(e_svar,true);
    Rcpp::XPtr< std::vector<int> > extern_sc(e_sc,true);
    Rcpp::XPtr< std::vector<double> > extern_stheta(e_stheta,true);

    ret["sts"]=extern_sts;
    ret["sid"]=extern_sid;
    ret["svar"]=extern_svar;
    ret["sc"]=extern_sc;
    ret["stheta"]=extern_stheta;
    ret["s.train"]=sdraws;
    ret["s.trees"]=Rcpp::CharacterVector(strees.str());
  }
//  else if(drawMuTau==0) ret["sigma"]=sddraws;

  if(drawMuTau>0) {
    ret["dpalpha"]=dalpha;
    ret["dpn"]=dnpart;
    ret["dpmu"]=dmu;
    ret["dpsd"]=dsig;
    ret["dpC"]=dpC;
    ret["dpmu."]=dpMU;
    ret["dpsd."]=dpSD;
    //ret["dpwt."]=dpWT;
  }

  if(np) {
    ret["f.test"]=mpred;
    //if(draws) 
    ret["s.test"]=spred;
  }

  ret["z.train"]=zdraws;

  // summary statistics
  if(summarystats) {
    //unsigned int varcount[p];
    //    std::vector<unsigned int> varcount(p);
/*
    Rcpp::IntegerVector varcount(p);
    for(size_t i=0;i<p;i++) varcount[i]=0;
    unsigned int tmaxd=0;
    unsigned int tmind=0;
    double tavgd=0.0;
*/

    ambm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
    tavgd/=(double)(nd*m);
    ret["f.tavgd"]=tavgd;
    ret["f.tmaxd"]=tmaxd;
    ret["f.tmind"]=tmind;
    ret["f.varcount"]=fvc;
    //ret["f.varcount"]=Rcpp::clone(varcount);
    //Rcpp::NumericVector vc(p);
    //Rcpp::IntegerVector vc(p);
    //for(size_t i=0;i<p;i++) vc[i]=(int)varcount[i];
    //ret["f.varcount"]=vc;

    //if(draws) 
    {
      //for(size_t i=0;i<p;i++) varcount[i]=0;
      //tmaxd=0; tmind=0; tavgd=0.0;
      psbm.getstats(&varcount[0],&tavgd,&tmaxd,&tmind);
      tavgd/=(double)(nd*mh);
      ret["s.tavgd"]=tavgd;
      ret["s.tmaxd"]=tmaxd;
      ret["s.tmind"]=tmind;
      ret["s.varcount"]=svc;
      //ret["s.varcount"]=varcount;
      //Rcpp::NumericVector sdvc(p);
      //Rcpp::IntegerVector sdvc(p);
      //for(size_t i=0;i<p;i++) sdvc[i]=(int)varcount[i];
      //ret["s.varcount"]=sdvc;
    }
  }

#if defined(RRNG)
  PutRNGstate();
#endif

  if(_f) delete [] _f;
  if(_s) delete [] _s;
  if(r) delete [] r;
  if(fun) delete [] fun;
  if(mvec) delete [] mvec;
  if(svec) delete [] svec;

  return ret;
}


