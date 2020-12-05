#ifndef GUARD_info_h
#define GUARD_info_h

#include "rrn.h"

//data
class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;}
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j)
   double *y; // ith y is *(y+i) or y[i]
};
//prior and mcmc
class pinfo
{
public:
   pinfo(): pbd(1.0),pb(.5),alpha(.95),mybeta(2.0),tau(1.0) {}
//mcmc info
   double pbd; //prob of birth/death
   double pb;  //prob of birth
//prior info
   double alpha;
   double mybeta;
   double tau;
   void pr() {
      cout << "pbd,pb: " << pbd << ", " << pb << std::endl;
      cout << "alpha,beta,tau: " << alpha << 
             ", " << mybeta << ", " << tau << std::endl;
   }
};

#endif
