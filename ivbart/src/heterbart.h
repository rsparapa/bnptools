#ifndef GUARD_heterbart_h
#define GUARD_heterbart_h
#include "bart.h"

class heterbart : public bart
{
public:
   heterbart():bart() {}
   heterbart(size_t m):bart(m) {}
   void pr();
   void draw(double *sigma, rn& gen);
};

//#include "heterbart.h"
#include "heterbartfuns.h"
#include "heterbd.h"

//--------------------------------------------------
void heterbart::pr()
{
   cout << "+++++heterbart object:\n";
   bart::pr();
}
//--------------------------------------------------
void heterbart::draw(double *sigma, rn& gen)
{
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      heterbd(t[j],xi,di,pi,sigma,gen);
      heterdrmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }
}

#endif
