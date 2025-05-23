/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2018-2025 Robert McCulloch and Rodney Sparapani
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

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)
#define XV(a, b) xv(a, b)
#define IMPUTE_DRAW1(a, b) impute_draw1(a, b)
#define IMPUTE_DRAW2(a, b) impute_draw2(a, b)

RcppExport SEXP cgbart(
   SEXP _type,          //1:wbart, 2:pbart, 3:lbart
   SEXP _in,            //number of observations in training data
   SEXP _ip,            //dimension of x
   SEXP _inp,           //number of observations in test data
   SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,            //y, train,  nx1
   SEXP _ixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
   //SEXP _z_train,
   SEXP _im,            //number of trees
   SEXP _inc,           //number of cut points
   SEXP _ind,           //number of kept draws (except for thinnning ..)
   SEXP _iburn,         //number of burn-in draws skipped
   SEXP _ithin,         //thinning
   SEXP _ipower,
   SEXP _ibase,
   SEXP _Offset,
   SEXP _itau,
   SEXP _inu,
   SEXP _ilambda,
   SEXP _isigest,
   SEXP _iw,
//   SEXP _treeinit,
//   SEXP _itrees,
   SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,            //param a for sparsity prior
   SEXP _ib,            //param b for sparsity prior
   SEXP _irho,          //param rho for sparsity prior (default to p)
   SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP _varprob,
   SEXP _inprintevery,
   SEXP _Xinfo,
   SEXP _shards,
   SEXP _verbose,
   SEXP _impute_mult, // integer vector of column indicators for missing covariates
   SEXP _impute_miss, // integer vector of row indicators for missing values
   SEXP _impute_prior // matrix of prior missing imputation probability
)
{
   //process args
   int type = Rcpp::as<int>(_type), 
     shards = Rcpp::as<int>(_shards),
     verbose = Rcpp::as<int>(_verbose);
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
 
/*
  // the pointer from a vector should be equivalent to the pointer from a matrix
  // since, as far as R is concerned, a matrix is just a fancy vector
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
*/
   
   Rcpp::NumericMatrix xv(_ix); // transposed:  p rows, n columns
   double *ix = &XV(0, 0);

/*
   size_t K;
   Rcpp::List impute(_impute);
   Rcpp::IntegerVector impute_mult, impute_miss;
//   Rcpp::NumericVector impute_y;
//   Rcpp::NumericMatrix impute_prior;
   if(impute.size()==0) K=0;
   else {
     K=impute_mult.size();
     impute_mult=impute["mult"];
     impute_miss=impute["miss"];
//     impute_prior=impute["prior"];
//     impute_y=impute["y"];
   }
*/
   Rcpp::IntegerVector impute_mult(_impute_mult); // integer vector of column indicators for missing covariates
   size_t K = impute_mult.size(); // number of columns to impute
   Rcpp::IntegerVector impute_miss(_impute_miss); // length n: integer vector of row indicators for missing values
   Rcpp::NumericMatrix impute_prior(_impute_prior); // n X K: matrix of prior missing imputation probability
   Rcpp::NumericVector impute_post(K); // length K: double vector of posterior missing imputation probability
   Rcpp::NumericVector impute_fhat(K); 
   double *impute_fhat_ptr = 0, *impute_Xrow_ptr = 0;
   if(K>0) impute_fhat_ptr = &impute_fhat[0];
   
   Rcpp::NumericVector  yv(_iy); 
   double *iy = &yv[0];
   //Rcpp::NumericVector  xpv(_ixp);
   //double *ixp = &xpv[0];
   Rcpp::NumericMatrix xpv(_ixp);
   double *ixp = nullptr;
   if(np>0) ixp = &xpv[0];
   //Rcpp::NumericVector z_train(_z_train);
   //double *z = &z_train[0];
   //Rcpp::IntegerVector z_draw(_z_draw);
   size_t m = Rcpp::as<int>(_im);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   int burn = Rcpp::as<int>(_iburn);
   size_t thin = Rcpp::as<int>(_ithin);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   double Offset = Rcpp::as<double>(_Offset);
   double tau = Rcpp::as<double>(_itau);
   double nu = Rcpp::as<double>(_inu);
   double lambda = Rcpp::as<double>(_ilambda);
   double sigma=Rcpp::as<double>(_isigest);
   Rcpp::NumericVector  wv(_iw); 
   double *iw = &wv[0];
//   int treeinit=Rcpp::as<int>(_treeinit); 
   //Rcpp::CharacterVector itrees(_itrees); 
   //std::string itv(itrees[0]);
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
//   Rcpp::NumericVector irho(_irho);
//   double *rho = &irho[0];
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   Rcpp::NumericVector varprob(_varprob);
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
   int *grp = &_grp[0];
//   if(rho==0.) for(size_t i=0; i<p; ++i) rho += (1./grp[i]);
   size_t nkeeptrain = nd/thin;     //Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = nd/thin;      //Rcpp::as<int>(_inkeeptest);
   size_t nkeeptreedraws = nd/thin; //Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
   Rcpp::NumericVector sdraw(burn+nkeeptrain);
   //Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericVector accept(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
   Rcpp::NumericMatrix impute_draw1(nkeeptrain, n);
   Rcpp::NumericMatrix impute_draw2(nkeeptrain, n);
   //Rcpp::NumericVector impute_draw(Rcpp::Dimension(nkeeptrain, n, K));

   //random number generation
   arn gen;

   heterbart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       for(int j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
       //for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void cgbart(
   int type,            //1:wbart, 2:pbart, 3:lbart
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int *numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   int burn,		//number of burn-in draws skipped
   size_t thin,		//thinning
   double mybeta,
   double alpha,
   double Offset,
   double tau,
   double nu,
   double lambda,
   double sigma,
   double* iw,
   bool dart,           //dart prior: true(1)=yes, false(0)=no   
   double theta,
   double omega, 
   int* grp,
   double a,		//param a for sparsity prior                                          
   double b,		//param b for sparsity prior                                          
//   double rho,		//param rho for sparsity prior (default to p)                         
   bool aug,		//categorical strategy: true(1)=data augment false(0)=degenerate trees
//   size_t nkeeptrain, //   size_t nkeeptest, //   size_t nkeeptreedraws,
   size_t printevery,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
   double* sdraw,
   double* _trdraw,
   double* _tedraw
)
{
   //return data structures (using C++)
   size_t nkeeptrain=nd/thin, nkeeptest=nd/thin, nkeeptreedraws=nd/thin;
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];

   //matrix to return dart posteriors (counts and probs)
   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2);
   heterbart bm(m);
#endif

   /* multiple imputation hot deck implementation
      this will cause trouble for multiple threads
      so it must be done prior to calling C++ */
/*
   bool hotdeck=false;

   std::vector<int> _missing(n*p);
   std::vector<int*> missing(n);
   std::vector<double*> x(n);
   Rcpp::NumericMatrix X(n, p); 

   for(size_t i=0; i<n; ++i) {
     missing[i]=&_missing[i*p];
     x[i]=&ix[i*p];
     for(size_t j=0; j<p; ++j) 
       if(x[i][j]!=x[i][j]) {
	 hotdeck=true;
	 missing[i][j]=1;
       }
       else missing[i][j]=0;
   }     

   if(hotdeck) {
     for(size_t i=0; i<n; ++i) {
       for(size_t j=0; j<p; ++j) { 
	 if(missing[i][j]==1) {
	   while(x[i][j]!=x[i][j]) {
	     size_t k=n*gen.uniform();
	     x[i][j]=x[k][j];
	   }
	 }
	 X(i, j)=x[i][j]; 
       }
     }    
   } 
*/

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;

   size_t skiptr=thin, skipte=thin, skiptreedraws=thin;
/*
   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;
*/

   //--------------------------------------------------
   //print args

   if(verbose==1) {
     printf("*****Calling gbart: type=%d\n", type);
     printf("*****Data:\n");
     printf("n,p,np: %zu, %zu, %zu\n",n,p,np);
     printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
     printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
     //   if(hotdeck) 
 //printf("warning: missing elements in x multiply imputed with hot decking\n");
     if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
     printf("*****Number of Trees: %zu\n",m);
     printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
     printf("*****burn,nd,thin: %d,%zu,%zu\n",burn,nd,thin);
     //   printf("*****Value of treeinit: %zu\n", treeinit);
 // printf("Prior:\nbeta,alpha,tau,nu,lambda,offset: %lf,%lf,%lf,%lf,%lf,%lf\n",
     //                    mybeta,alpha,tau,nu,lambda,Offset);
     cout << "*****Prior:beta,alpha,tau,nu,lambda,offset,shards:\n" 
	  << mybeta << ',' << alpha << ',' << tau << ',' << nu << ',' 
          << lambda << ',' << Offset << ',' << shards << endl;
     if(type==1) {
       printf("*****sigma: %lf\n",sigma);
       printf("*****w (weights): %lf ... %lf\n",iw[0],iw[n-1]);
     }
     if(K>0) {
       cout << "*****Missing imputation row indices:\n index 0=" << impute_miss[0] << ','
	    << "index n-1=" << impute_miss[n-1] << endl;
       cout << "*****Missing imputation column indices:\n index 0=" << impute_mult[0] << ','
	    << "index K-1=" << impute_mult[K-1] << endl;
// cout << "*****Missing imputation probability: prob[0]=" << impute_prior[0] 
// << ',' << "prob[K-1]=" << impute_prior[K-1] << endl;
     }
     //printf("*****nkeeptrain,nkeeptest: %zu, %zu\n",nkeeptrain,nkeeptest);
     //printf("*****printevery: %zu\n",printevery);
 cout << "*****Dirichlet:sparse,theta,omega,rho,a,b,augment,grp[0],grp[p-1]:\n" 
      << dart << ',' << theta << ',' << omega << ',' << rho << ',' << a << ',' 
      << b << ',' << aug << ',' << grp[0] << ',' << grp[p-1] << endl;
   }
   //--------------------------------------------------
   //create temporaries
   double df=n+nu;
   double *z = new double[n]; 
   double *svec = new double[n]; 
   double *sign = NULL;
   if(type!=1) sign = new double[n]; 
   Rcpp::IntegerVector prevXV(p);

   for(size_t i=0; i<n; i++) {
     if(type==1) {
       svec[i] = iw[i]*sigma; 
       z[i] = iy[i];
     }
     else {
       svec[i] = iw[i];
       //svec[i] = 1.;
       if(iy[i]==0) sign[i] = -1.;
       else sign[i] = 1.;
       z[i] = sign[i];
       //z[i] = sign[i]*iw[i];
     }
     if(K>0) {
       if(impute_miss[i]==1) {
       //if(impute_miss[i]>0) {
	 for(size_t j=0; j<K; j++) {
	   XV(impute_mult[j], i)=0;
	   prevXV[impute_mult[j]]=0;
	 }
	 size_t k;
	 Rcpp::NumericVector impute_prob(impute_prior.row(i));
	 k=gen.rcat(impute_prob); // use prior prob only
	 XV(impute_mult[k], i)=1;
	 prevXV[impute_mult[k]]=1;
       }
     else if(impute_miss[i]==2) 
       for(size_t j=0; j<K; j++) 
	 XV(impute_mult[j], i)=prevXV[impute_mult[j]];
     }
   }
   //--------------------------------------------------
   //set up BART model
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,z,numcut);
/*
   if(treeinit==1) {
     Rcpp::CharacterVector itrees(_itrees); 
     std::string itv(itrees[0]);
     bm.settree(itv);
   }
*/
   bm.setdart(a,b,grp,aug,dart,rho,theta,omega);
   //bm.setdart(a,b,rho,aug,dart);
   bm.setpv(&varprob[0]);

   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);
   ivarprb=bm.getpv();
      if(verbose==1) {
	cout << "*****Variable selection probability pv[0],pv[p-1]:\n"
        << ivarprb[0] << ',' << ivarprb[p-1] << endl;
	printf("\nMCMC\n");
      }

   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; 
   if(np) { fhattest = new double[np]; }

   //--------------------------------------------------
   //mcmc
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   bool keeptest,/*keeptestme*/keeptreedraw,
     type1sigest=(type==1 && lambda!=0.);

   time_t tp;
   int time1 = time(&tp), total=nd+burn;
   xinfo& xi = bm.getxinfo();

   for(int i=0;i<total;i++) {
   //for(size_t i=0;i<total;i++) {
      if(verbose==1 && i%printevery==0) 
	printf("done %d (out of %d)\n",i,total);
      //if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(svec,gen,shards);
      accept[i]=bm.getaccept();

      if(type1sigest) {
      //draw sigma
	double rss=0.;
	for(size_t k=0;k<n;k++) rss += pow((iy[k]-bm.f(k))/(iw[k]), 2.); 
	sigma = sqrt((nu*lambda + rss)/gen.chi_square(df));
	if(i<burn) sdraw[i]=sigma;
      }

      for(size_t k=0; k<n; k++) {
	if(type==1) svec[k]=iw[k]*sigma;
	else {
	  z[k]=sign[k]*rtnorm(sign[k]*bm.f(k), -sign[k]*Offset, svec[k], gen);
	  if(type==3) 
	    svec[k]=sqrt(draw_lambda_i(pow(svec[k], 2.), sign[k]*bm.f(k), 
				       1000, 1, gen));
	  }
      }

      if(K>0) {
	for(size_t k=0; k<n; ++k) {
	  if(impute_miss[k]==1) {
	  //if(impute_miss[k]>0) {
	    impute_Xrow_ptr=&XV(0, k);
	    impute_post=impute_prior.row(k);
	    for(size_t j=0; j<K; ++j) {
	      for(size_t h=0; h<K; ++h) XV(impute_mult[h], k)=0;
	      XV(impute_mult[j], k)=1;
	      bm.predict(p, 1, impute_Xrow_ptr, &impute_fhat_ptr[j]);
	      impute_post[j] *= R::dnorm(z[k], impute_fhat_ptr[j], svec[k], 0); 
	    }
	    for(size_t j=0; j<K; j++) { 
	      XV(impute_mult[j], k)=0;
	      prevXV[impute_mult[j]]=0;
	    }
	    size_t h;
	    h=gen.rcat(impute_post); 
	    XV(impute_mult[h], k)=1;
	    prevXV[impute_mult[h]]=1;
	  }
	  else if(impute_miss[k]==2) 
	    for(size_t j=0; j<K; j++) 
	      XV(impute_mult[j], k)=prevXV[impute_mult[j]];
	}
      }

/*
      if(hotdeck) {
	//draw x
	for(size_t h=0; h<n; ++h) {
	  for(size_t j=0; j<p; ++j) {
	    if(missing[h][j]==1) {
	      size_t k=n*gen.uniform();
	      x[h][j]=x[k][j];
	    }
	  }
	}    
      } 
*/

      if(i>=burn) {
	//size_t idcnt=0;
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
	   if(type1sigest) sdraw[trcnt+burn]=sigma;
            for(size_t k=0;k<n;k++) {
	      TRDRAW(trcnt,k)=Offset+bm.f(k);
	      if(K>0) {
		if(impute_miss[k]==1) {
		  IMPUTE_DRAW1(trcnt,k)=XV(impute_mult[0], k);
		  IMPUTE_DRAW2(trcnt,k)=XV(impute_mult[1], k);
		} else {
		  IMPUTE_DRAW1(trcnt,k)=2;
		  IMPUTE_DRAW2(trcnt,k)=2;
		}
	      }
/*
	      for(size_t j=0; j<K; ++j) {
		impute_draw[idcnt]=XV(impute_mult[j], k);
		idcnt+=1;
		//impute_draw(trcnt, k, j)=XV(impute_mult[j], k);
	      }
*/
	    }
	    trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=Offset+fhattest[k];
            tecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
	    size_t k=(i-burn)/skiptreedraws;
            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);

	      #ifndef NoRcpp
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
//	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t h=0;h<p;h++){
	      varcnt(k,h)=ivarcnt[h];
	      varprb(k,h)=ivarprb[h];
	    }
/* dangerous re-use of j
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	    }
*/
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif
	    }
         }
      }
   }
   if(verbose==1) {
     int time2 = time(&tp);
     printf("time: %ds\n",time2-time1);
     printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);
   }
   if(fhattest) delete[] fhattest;
   delete[] z;
   delete[] svec;
   if(type!=1) delete[] sign;

#ifndef NoRcpp
   //return list
   Rcpp::List ret;
//   ret["X"]=X; 
   if(type1sigest) ret["sigma"]=sdraw;
   ret["yhat.train"]=trdraw;
   ret["yhat.test"]=tedraw;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;
   ret["accept"]=accept;

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   ret["treedraws"] = treesL;

   if(K>0) {
     ret["impute.draw1"]=impute_draw1;
     ret["impute.draw2"]=impute_draw2;
   }
   //if(K>0) ret["impute.mult"]=impute_draw;
   
   return ret;
#else

#endif

}
