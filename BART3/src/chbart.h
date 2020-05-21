/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2018 Robert McCulloch and Rodney Sparapani
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
#define TRDRAW2(a, b) trdraw2(a, b)
#define TEDRAW2(a, b) tedraw2(a, b)

RcppExport SEXP chbart(
//   SEXP _type,          //1:wbart, 2:pbart, 3:lbart
   SEXP _in,            //number of observations in training data
   SEXP _ip,            //dimension of x
   SEXP _inp,           //number of observations in test data
   SEXP _ix,            //x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _iy,            //y, train,  nx1
   SEXP _ixp,           //x, test, pxnp (transposed so rows are contiguous in memory)
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
   SEXP _idart,         //dart prior: true(1)=yes, false(0)=no
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,            //param a for sparsity prior
   SEXP _ib,            //param b for sparsity prior
   SEXP _irho,          //param rho for sparsity prior (default to p)
   SEXP _iaug,          //categorical strategy: true(1)=data augment false(0)=degenerate trees
   SEXP _inprintevery,
   SEXP _Xinfo
)
{
   //process args
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::NumericVector  yv(_iy); 
   double *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   Rcpp::IntegerVector m(_im);
//   size_t m = Rcpp::as<int>(_im);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   size_t thin = Rcpp::as<int>(_ithin);
/*
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   double Offset = Rcpp::as<double>(_Offset);
   double tau = Rcpp::as<double>(_itau);
   double nu = Rcpp::as<double>(_inu);
   double lambda = Rcpp::as<double>(_ilambda);
   double sigma=Rcpp::as<double>(_isigest);
*/
   Rcpp::NumericVector mybeta(_ipower);
   Rcpp::NumericVector alpha(_ibase);
   Rcpp::NumericVector Offset(_Offset);
   Rcpp::NumericVector tau(_itau);
   Rcpp::NumericVector nu(_inu);
   Rcpp::NumericVector lambda(_ilambda);
   Rcpp::NumericVector sigma(_isigest);

   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
//   int *grp = &_grp[0];
   size_t nkeeptrain = nd/thin;     //Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = nd/thin;      //Rcpp::as<int>(_inkeeptest);
   size_t nkeeptreedraws = nd/thin; //Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
   Rcpp::NumericMatrix varprb2(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt2(nkeeptreedraws,p);
   Rcpp::NumericMatrix trdraw2(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw2(nkeeptest,np);

   //random number generation
   arn gen;

   heterbart bm(m[0]); // f(x)
   bart bm2(m[1]);// s(x)

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
     bm2.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]

void chbart(
//   int type,            //1:wbart, 2:pbart, 3:lbart
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int *numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   size_t thin,		//thinning
   double mybeta,
   double alpha,
   double Offset,
   double tau,
   double nu,
   double lambda,
   double sigma,
   bool dart,           //dart prior: true(1)=yes, false(0)=no   
   double theta,
   double omega, 
   int* grp,
   double a,		//param a for sparsity prior                                          
   double b,		//param b for sparsity prior                                          
   double rho,		//param rho for sparsity prior (default to p)                         
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
   std::vector< std::vector<size_t> > varcnt2;
   std::vector< std::vector<double> > varprb2;

   //random number generation
   arn gen(n1, n2);
   heterbart bm(m[0]);
   bart bm2(m[1]);
#endif

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m[0] << " " << p << endl;
   std::stringstream treess2;  //string stream to write trees to
   treess2.precision(10);
   treess2 << nkeeptreedraws << " " << m[1] << " " << p << endl;

   printf("*****Calling hbart\n");

   size_t skiptr=thin, skipte=thin, skiptreedraws=thin;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
   printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu, %zu\n",m[0],m[1]);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn,nd,thin: %zu,%zu,%zu\n",burn,nd,thin);
   cout << "*****Prior:beta,alpha,tau,nu,lambda,offset: " 
	<< mybeta[0] << ',' << alpha[0] << ',' << tau[0] << ',' 
        << nu[0] << ',' << lambda[0] << ',' << Offset[0] << endl;
   cout << "*****Prior2:beta,alpha,tau,nu,lambda,offset: " 
	<< mybeta[1] << ',' << alpha[1] << ',' << tau[1] << ',' 
        << nu[1] << ',' << lambda[1] << ',' << Offset[1] << endl;
   printf("*****sigma: %lf\n",sigma[1]);
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   printf("*****printevery: %zu\n",printevery);

   //--------------------------------------------------
   //create temporaries
   double df=n+nu[1];
   double *z = new double[n]; 
   double *svec = new double[n]; 

   for(size_t i=0; i<n; i++) {
       svec[i] = sigma[0]; 
   }
   //--------------------------------------------------
   //set up BART model
   bm.setprior(alpha[0],mybeta[0],tau[0]);
   bm.setdata(p,n,ix,iy,numcut);
   bm.setdart(a,b,rho,aug,dart);
   bm2.setprior(alpha[1],mybeta[1],tau[1]);
   bm2.setdata(p,n,ix,z,numcut);
   bm2.setdart(a,b,rho,aug,dart);

   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);
   std::vector<double> ivarprb2(p,0.);
   std::vector<size_t> ivarcnt2(p,0);
   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; 
   double* fhattest2=0; 
   if(np) { 
     fhattest = new double[np]; 
     fhattest2 = new double[np]; 
   }

   //--------------------------------------------------
   //mcmc
   printf("\nMCMC\n");
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   bool keeptest,/*keeptestme*/keeptreedraw;

   time_t tp;
   int time1 = time(&tp), total=nd+burn;
   xinfo& xi = bm.getxinfo();

   for(size_t i=0;i<total;i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      //if(i%printevery==0) printf("%22zu/%zu\r",i,total);
      if(i==(burn/2)&&dart) {
	bm.startdart();
	bm2.startdart();
      }

      //draw bart
      bm.draw(svec,gen);

      for(size_t k=0;k<n;k++) {
	z[k] = log(pow(iy[k]-bm.f(k), 2.)); 
      }

      //draw hbart
      bm2.draw(sigma[1],gen);

      //draw sigma
	double rss=0.;
	for(size_t k=0;k<n;k++) {
	  rss += pow(z[k]-bm2.f(k), 2.);
	  svec[k] = exp(bm2.f(k));
	}
	sigma[1] = sqrt((nu[1]*lambda[1] + rss)/gen.chi_square(df));
	sdraw[i]=sigma[1];

      if(i>=burn) {
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
            for(size_t k=0;k<n;k++) {
	      TRDRAW(trcnt,k)=Offset[0]+bm.f(k);
	      TRDRAW2(trcnt,k)=exp(Offset[1]+bm2.f(k));
	    }
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
	   bm2.predict(p,np,ixp,fhattest2);
            for(size_t k=0;k<np;k++) {
	      TEDRAW(tecnt,k)=Offset[0]+fhattest[k];
	      TEDRAW2(tecnt,k)=exp(Offset[0]+fhattest2[k]);
	    }
            tecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {
	    size_t k=(i-burn)/skiptreedraws;
            for(size_t j=0;j<m[0];j++) {
	      treess << bm.gettree(j);
	      #ifndef NoRcpp
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    for(size_t h=0;h<p;h++){
	      varcnt(k,h)=ivarcnt[h];
	      varprb(k,h)=ivarprb[h];
	    }
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif
	    }
	    for(size_t j=0;j<m[1];j++) {
	      treess2 << bm2.gettree(j);
	      #ifndef NoRcpp
	    ivarcnt2=bm2.getnv();
	    ivarprb2=bm2.getpv();
	    for(size_t h=0;h<p;h++){
	      varcnt2(k,h)=ivarcnt2[h];
	      varprb2(k,h)=ivarprb2[h];
	    }
            #else
	    varcnt2.push_back(bm2.getnv());
	    varprb2.push_back(bm2.getpv());
	    #endif
	    }
         }
      }
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);

   if(fhattest) delete[] fhattest;
   if(fhattest2) delete[] fhattest2;
   delete[] z;
   delete[] svec;

#ifndef NoRcpp
   //return list
   Rcpp::List ret;
//   ret["X"]=X; 
   ret["sigma"]=sdraw;
   ret["yhat.train"]=trdraw;
   ret["yhat.test"]=tedraw;
   ret["shat.train"]=trdraw2;
   ret["shat.test"]=tedraw2;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;
   ret["varcount2"]=varcnt2;
   ret["varprob2"]=varprb2;

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   treesL["trees2"]=Rcpp::CharacterVector(treess2.str());
   ret["treedraws"] = treesL;

   return ret;
#else

#endif

}
