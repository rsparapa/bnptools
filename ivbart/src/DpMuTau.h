#ifndef GUARD_DPMUTAU
#define GUARD_DPMUTAU

#include "Dp.h"

class DpMuTau: public Dp {
//Dp for Normal/gamma conjugate prior for y~N(mu,1/tau), theta is (mu,tau) and eta is empty

public:
   //--------------------------------------------------
   //const/dest
   DpMuTau(): Dp(), muo(0.0), ko(1.0), alphao(5.0),betao(5.0) {setlc();}
   DpMuTau(size_t _n, Dp::dv _theta): Dp(_n,_theta), muo(0.0), ko(1.0), alphao(5.0),betao(5.0) {setlc();}
   ~DpMuTau() {}

   //get/set
   void setData(double *y) {Dp::setData(1,y);} //call setData from Dp with p=1
   void setPrior(double muo, double ko, double alphao, double betao) {
      cout << "into setPrior:\n";
      this->muo = muo;
      this->ko = ko;
      this->alphao = alphao;
      this->betao = betao;
      setlc();
   }

   //-------------------------------
   //the virtuals
   void toscreen();
   double f(double* y, dv& theta, double* eta);
   double qo(double* y, double* eta);
   dv draw_one_theta(std::list<size_t>& ind, rn& gen);

   //utility
   double setlc();

   //--------------------------------------------------
   //private

   //------------------------------
   //prior parameters
   double muo,ko,alphao,betao;

   
   //------------------------------
   //working parameters
   double lc;  //all of log(qo) that does not need y, just depends on prior

};


//#include "DpMuTau.h"

//################################################################################
//### toscreen
void DpMuTau::toscreen() {
   cout << "########DpMuTau object, y~N(mu,1/tau); theta is (mu,tau) and eta is empty:\n";
   cout <<"prior:\n";
   cout << "\tmuo,ko,alphao,betao: " << muo << ", " << ko << ", " << alphao << ", " << betao << std::endl;
   cout << "\tlc: " << lc << std::endl;
   Dp::toscreen();
}
//################################################################################
//### draw_one_theta
Dp::dv DpMuTau::draw_one_theta(std::list<size_t>& ind, rn& gen)
{
   //cout << "&&&&&&&&&&&&&& into draw_one_theta\n";

   //--------------------------------------------------
   //return value
   Dp::dv rtheta{0.0,1.0};

   //--------------------------------------------------
   //sufficient statistics
   int nn = ind.size();
   //cout << "nn: " << nn << std::endl;

   //mean
   double ybar=0.0;
   double yy;
   for(thetaGP::iiter i=ind.begin();i!=ind.end();i++) {
      yy = y[*i];
      ybar += yy;
   }
   ybar /= nn;

   //sum of squares
   double S=0.0;
   for(thetaGP::iiter i=ind.begin();i!=ind.end();i++) {
      yy = y[*i];
      S += (yy-ybar)*(yy-ybar);
   }
   //cout << "ybar: " << ybar << std::endl;
   //cout << "S: " << S << std::endl;

   //--------------------------------------------------
   //posterior parameters
   double alpha1 = alphao + nn/2.0;
   double k1 = ko+nn;
   double mu1 = (nn*ybar+ko*muo)/(nn+ko);
   double beta1 = betao + .5*(S+(muo-ybar)*(muo-ybar)*nn*ko/(nn+ko));
   //cout << "alpha1,k1,mu1,beta1: " << alpha1 << ", " << k1 << ", " << mu1 << ", " << beta1 << std::endl; 
   
   //--------------------------------------------------
   //draw from posterior
   gen.set_alpha(alpha1);
   double tau = gen.gamma()/beta1;
   double sigma = sqrt(1.0/(k1*tau));
   double mu = mu1 + sigma*gen.normal();
   rtheta[0]=mu; rtheta[1]=tau;

   //cout << "&&&&&&&&&&&&&& out draw_one_theta\n";
   return rtheta;
}
//################################################################################
//### f(y,theta), theta=(mu,tau), no eta
double DpMuTau::f(double* y, Dp::dv& theta, double* eta)
{
   //recall: (mu,tau) = theta
   double sigma = sqrt(1.0/theta[1]);
   double mu = theta[0];
   double r = (*y-mu)/sigma;
   return exp(-.5*r*r)/(sigma*RTPI);
}
//################################################################################
//### qo
double DpMuTau::qo(double* y, double* eta)
{

   double beta1=betao;
   double d2 = (muo - *y);
   d2 = d2*d2;

   beta1 += .5*ko*d2/(1+ko);
   double logqo = lc - (alphao+.5)*log(beta1);

   return exp(logqo);
}
//--------------------------------------------------
double DpMuTau::setlc()
{
   //lc = 14.0;
   lc = alphao*log(betao) + .5*log(ko) - logam(alphao);
   lc = lc + logam(alphao+.5) -.5*log(1+ko);
   lc = lc - log(RTPI);
}



#endif

