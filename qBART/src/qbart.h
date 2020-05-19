#ifndef GUARD_QBART_h
#define GUARD_QBART_h

#include "BART3/bart.h"
#include "qbartfuns.h"
#include "qbd.h"

class qbart : public bart
{
public:
  qbart():bart(),q(0) {}
  qbart(size_t m):bart(m),q(0) {}
  void setdata(size_t p, size_t n, double *x, double *y, int *q, size_t nc=100);
  void setdata(size_t p, size_t n, double *x, double *y, int *q, int* nc);
  void draw(double sigma, rn& gen);
  void pr();

protected:
  int *q; //individual cure status
};


void qbart::setdata(size_t p, size_t n, double *x, double *y, int *q, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  this->setdata(p, n, x, y, q, nc);
  delete [] nc;
}

void qbart::setdata(size_t p, size_t n, double *x, double *y, int *q, int *nc)
{
   this->p=p; this->n=n; this->x=x; this->y=y; this->q=q;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

   if(allfit) delete[] allfit;
   allfit = new double[n];
   predict(p,n,x,allfit);

   if(r) delete[] r;
   r = new double[n];

   if(ftemp) delete[] ftemp;
   ftemp = new double[n];

   di.n=n; di.p=p; di.x = &x[0]; di.y=r; di.q=q;
   for(size_t j=0;j<p;j++){
     nv.push_back(0);
     pv.push_back(1/(double)p);
   }
}

void qbart::draw(double sigma, rn& gen)
{
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      qbd(t[j],xi,di,pi,sigma,nv,pv,aug,gen);
      qdrmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }
   if(dartOn) {
     draw_s(nv,lpv,theta,gen);
     draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
     for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
   }
}

void qbart::pr(){
  cout << "+++++++latent cure status\n";
  bart::pr();
}

#endif
