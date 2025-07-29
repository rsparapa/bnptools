#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>

#include "rng.h"
#include "tree.h"
#include "info.h"
#include "funs.h"
#include "bd.h"

using std::cout;
using std::endl;

void xitofile(std::ofstream& fnm, xinfo xi)
{
   fnm << xi.size() << endl;
   for(uint i=0;i< xi.size();i++) {
      fnm << xi[i].size() << endl;
      for(uint j=0;j<xi[i].size();j++) fnm << xi[i][j] << "  ";
      fnm << endl;
   }
   fnm << endl;
}

int main(int argc, char** argv)
{
   cout << "\n*****Into bart main\n";
   if(argc<3) {cout << "error need at least 2 arguments\n"; return 1;}
   double xtemp,ytemp; //used to read in x and y without knowing how many there are.

   //--------------------------------------------------
   //random number generation
   uint seed=99;
   RNG gen(seed); //this one random number generator is used in all draws
   
   //--------------------------------------------------
   //read in data
   //read y NOTE: we assume y is already centered at 0.
   std::ifstream yf(argv[2]); //file to read y from
   std::vector<double> y; //storage for y
   double miny = INFINITY; //use range of y to calibrate prior for bottom node mu's
   double maxy = -INFINITY;
   sinfo allys;       //sufficient stats for all of y, use to initialize the bart trees.
   while(yf >> ytemp) {
      y.push_back(ytemp);
      if(ytemp<miny) miny=ytemp;
      if(ytemp>maxy) maxy=ytemp;
      allys.sy += ytemp; // sum of y
      allys.sy2 += ytemp*ytemp; // sum of y^2
   }
   size_t n = y.size();
   if(n<1) {
      cout << "error n<1\n";
      return 1;
   }
   allys.n = n;
   cout << "\ny read in:\n";
   cout << "n: " << n << endl;
   cout << "y first and last:\n";
   cout << y[0] << ", " << y[n-1] << endl;
   double ybar = allys.sy/n; //sample mean
   double shat = sqrt((allys.sy2-n*ybar*ybar)/(n-1)); //sample standard deviation
   cout << "ybar,shat: " << ybar << ", " << shat <<  endl;

   //read x
   //the n*p numbers for x are stored as the p for first obs, then p for second, and so on.
   std::ifstream xf(argv[1]);  //file to read x from
   std::vector<double> x;
   while(xf >> xtemp) 
      x.push_back(xtemp);
   size_t p = x.size()/n;
   if(x.size() != n*p) {
      cout << "error: input x file has wrong number of values\n"; 
      return 1;
   }
   cout << "\nx read in:\n";
   cout << "p: " << p << endl;
   cout << "first row: " <<  x[0] << " ...  " << x[p-1] << endl;
   cout << "last row: " << x[(n-1)*p] << " ...  " << x[n*p-1] << endl;

   //x for predictions
   dinfo dip; //data information for prediction
   dip.n=0;
   std::vector<double> xp;     //prediction observations, stored like x
   if(argc>3) {
      std::ifstream xpf(argv[3]); //file to read x to predict at from
      while(xpf >> xtemp) 
         xp.push_back(xtemp);
      size_t np = xp.size()/p;
      if(xp.size() != np*p) {
         cout << "error, wrong number of elements in prediction data set\n";
         return 1;
      }
      if(np) 
         dip.n=np; dip.p=p; dip.x = &xp[0]; dip.y=0; //there are no y's!
   }
   cout << "\nx for prediction read in:\n";
   cout << "n: " << dip.n << endl;
   if(dip.n) {
      cout << "first row: " <<  dip.x[0] << " ...  " << dip.x[p-1] << endl;
      cout << "last row: " << dip.x[(dip.n-1)*p] << " ...  " << dip.x[dip.n*p-1] << endl;
   }

   //--------------------------------------------------
   //optionally read in additional arguments
   size_t burn = 100; //number of mcmc iterations called burn-in
   size_t nd = 1000; //number of mcmc iterations
   size_t m=200; //number of trees in BART sum
   double lambda = 1.0; //this one really needs to be set
   double nu = 3.0;
   double kfac=2.0;
   if(argc>4) nd = atoi(argv[4]);
   if(argc>5) lambda = atof(argv[5]);
   if(argc>6) burn = atoi(argv[6]);
   if(argc>7) m = atoi(argv[7]);
   if(argc>8) nu = atof(argv[8]);
   if(argc>9) kfac = atof(argv[9]);

   cout <<"\nburn,nd,number of trees: " << burn << ", " << nd << ", " << m << endl;
   cout <<"\nlambda,nu,kfac: " << lambda << ", " << nu << ", " << kfac << endl;

   //--------------------------------------------------
   //x cutpoints
   xinfo xi;
   size_t nc=100; //100 equally spaced cutpoints from min to max.
   makexinfo(p,n,&x[0],xi,nc);

   //--------------------------------------------------
   //trees
   std::vector<tree> t(m);
   for(size_t i=0;i<m;i++) t[i].setm(ybar/m); //if you sum the fit over the trees you get the fit.

   //--------------------------------------------------
  //prior and mcmc
   pinfo pi;
   pi.pbd=1.0; //prob of birth/death move
   pi.pb=.5; //prob of birth given  birth/death

   pi.alpha=.95; //prior prob a bot node splits is alpha/(1+d)^beta, d is depth of node
   pi.beta=2.0; //2 for bart means it is harder to build big trees.
   pi.tau=(maxy-miny)/(2*kfac*sqrt((double)m));
   pi.sigma=shat;

   cout << "\nalpha, beta: " << pi.alpha << ", " << pi.beta << endl;
   cout << "sigma, tau: " << pi.sigma << ", " << pi.tau << endl;

   //--------------------------------------------------
   //dinfo
   double* allfit = new double[n]; //sum of fit of all trees
   for(size_t i=0;i<n;i++) allfit[i]=ybar;
   double* r = new double[n]; //y-(allfit-ftemp) = y-allfit+ftemp
   double* ftemp = new double[n]; //fit of current tree
   dinfo di;
   di.n=n; di.p=p; di.x = &x[0]; di.y=r; //the y for each draw will be the residual 

   //--------------------------------------------------
   //storage for ouput
   //in sample fit
   double* pmean = new double[n]; //posterior mean of in-sample fit, sum draws,then divide
   for(size_t i=0;i<n;i++) pmean[i]=0.0;

   //out of sample fit
   double* ppredmean=0; //posterior mean for prediction
   double* fpredtemp=0; //temporary fit vector to compute prediction
   if(dip.n) {
      ppredmean = new double[dip.n];
      fpredtemp = new double[dip.n];
      for(size_t i=0;i<dip.n;i++) ppredmean[i]=0.0;
   }
   //for sigma draw
   double rss;  //residual sum of squares
   double restemp; //a residual
   std::ofstream bsdf("bart-sd.txt"); //note that we write all burn+nd draws to this file.

   std::ofstream btsf("bart-ave-tree-size.txt"); //note that we write all burn+nd draws to this file.
   double ats; //place for average tree size

   std::ofstream bnbf("bart-ave-num-bots.txt"); //note that we write all burn+nd draws to this file.
   double anb; //place for average number of bottom nodes

   //--------------------------------------------------
   //mcmc

   cout << "\nMCMC:\n";
   clock_t tp;
   tp = clock();
   for(size_t i=0;i<(nd+burn);i++) {
      if(i%100==0) cout << "i: " << i << endl;
      //draw trees
      for(size_t j=0;j<m;j++) {
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) {
            allfit[k] = allfit[k]-ftemp[k];
            r[k] = y[k]-allfit[k];
         }
         bd(t[j],xi,di,pi,gen);
         drmu(t[j],xi,di,pi,gen);
         fit(t[j],xi,di,ftemp);
         for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
      }
      //draw sigma
      rss=0.0;
      for(size_t k=0;k<n;k++) {restemp=y[k]-allfit[k]; rss += restemp*restemp;}
      pi.sigma = sqrt((nu*lambda + rss)/gen.chi_square(nu+n));
      bsdf << pi.sigma << endl;
      ats = 0.0; anb=0.0;
      for(size_t k=0;k<m;k++) {
          ats += t[k].treesize();
          anb += t[k].nbots();
      }
      btsf << ats/m << endl;
      bnbf << anb/m << endl;
      if(i>=burn) {
         for(size_t k=0;k<n;k++) pmean[k] += allfit[k];
         if(dip.n) {
            for(size_t j=0;j<m;j++) {
               fit(t[j],xi,dip,fpredtemp);
               for(size_t k=0;k<dip.n;k++) ppredmean[k] += fpredtemp[k];
            }
         }
      }
   }
   tp=clock()-tp;
   double thetime = (double)(tp)/(double)(CLOCKS_PER_SEC);
   cout << "time for loop: " << thetime << endl;
   std::ofstream timef("time.txt");
   timef << thetime << endl;

   for(size_t k=0;k<n;k++) pmean[k] /= (double)nd;
   std::ofstream bfitf("bart-fit.txt");
   for(size_t i=0;i<n;i++) bfitf << pmean[i] << endl;

   if(dip.n) {
      for(size_t k=0;k<dip.n;k++) ppredmean[k] /= (double)nd;
      std::ofstream bpredfitf("bart-pred-fit.txt");
      for(size_t i=0;i<dip.n;i++) bpredfitf << ppredmean[i] << endl;
   }

   std::ofstream treef("trees.txt");
   //treef << xi << endl; //the cutpoints
   xitofile(treef,xi);
   treef << m << endl;  //number of trees
   treef << p << endl;  //dimension of x's
   for(size_t j=0;j<m;j++) treef << t[j] << endl;  //all the trees

   return 0;
}
