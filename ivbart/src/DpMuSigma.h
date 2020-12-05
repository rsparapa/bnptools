#ifndef GUARD_DPMUSIGMA
#define GUARD_DPMUSIGMA

#include "Dp.h"

class DpMuSigma: public Dp {
//Dp for Normal/IWishart conjugate prior for y~N(mu,Sigma), theta is (mu,Chol(Sigma)) and eta is empty

public:
   //--------------------------------------------------
   //const/dest
   DpMuSigma(): Dp(), v(1.0),nu(2), a(1.0) {}
   DpMuSigma(size_t _n, Dp::dv _theta): Dp(_n,_theta), v(1.0),nu(2), a(1.0) {}
   ~DpMuSigma() {}

   //-------------------------------
   // get/set
   void setData(double *y) {Dp::setData(2,y);} //call setData from Dp with p=2
   void setPrior(double _v, double _nu, double _a) {v=_v; nu=_nu; a=_a;}

   //the virtuals
   void toscreen();
   double f(double* y, dv& theta, double* eta);
   double qo(double* y, double* eta);
   dv draw_one_theta(std::list<size_t>& ind, rn& gen);

   //public functions
   void center();


private:

   //------------------------------
   //prior parameters
   //Sigma^{-1} ~ Wishart_nu(V^{-1}), V = vI
   double v; 
   double nu;
   // mu | Sigma ~ N(0,Sigma/a)
   double a;
};


//#include "DpMuSigma.h"

//################################################################################
//### toscreen
void DpMuSigma::toscreen() {
   cout << "########DpMuSigma object, y~N(mu,Sigma),\n";
   cout << "\t  conjugate prior: mu~N(mubar,Sigma/a), Sigma^{-1} ~ Wishart(V^{-1})\n";
   cout <<"\t theta is (mu,Chol(Sigma)) and eta is empty:\n";
   cout <<"prior:\n";
   cout << "\tv: " << v << std::endl;
   cout << "\tnu: " << nu << std::endl;
   cout << "\ta: " << a << std::endl;
   Dp::toscreen();
}
//################################################################################
//### draw_one_theta
Dp::dv DpMuSigma::draw_one_theta(std::list<size_t>& ind, rn& gen)
{
   //cout << "&&&&&&&&&&&&&& into draw_one_theta\n";

   //--------------------------------------------------
   //return value
   Dp::dv rtheta{0.0,0.0,1.0,0.0,1.0};

   //--------------------------------------------------
   //sufficient statistics
   int nn = ind.size();
   //cout << "nn: " << nn << std::endl;

   //ybar
   double ybar1=0.0,ybar2=0.0;
   size_t ii;
   for(thetaGP::iiter i=ind.begin();i!=ind.end();i++) {
      ii = *i;
      ybar1 += *(y+ii*p);
      ybar2 += *(y+ii*p+1);
   }
   ybar1 /= nn;
   ybar2 /= nn;
   //cout << "yb: " << ybar1 << ", " << ybar2 << std::endl;

   double s11=0.0,s21=0.0,s22=0.0;
   double y1,y2;
   for(thetaGP::iiter i=ind.begin();i!=ind.end();i++) {
      ii = *i;
      y1 = *(y+ii*p);
      y2 = *(y+ii*p+1);
      s11 += (y1-ybar1)*(y1-ybar1);
      s22 += (y2-ybar2)*(y2-ybar2);
      s21 += (y2-ybar2)*(y1-ybar1);
   }
   //printf("s: %lf, %lf, %lf\n",s11,s21,s22);


   //--------------------------------------------------
   //posterior parameters
   double a1 = a + nn;
   double nu1 = nu + nn;
   double mu11 = (nn*ybar1)/(nn+a);
   double mu21 = (nn*ybar2)/(nn+a);
   double V11 = v + s11 + (nn*a*ybar1*ybar1)/(nn+a);
   double V22 = v + s22 + (nn*a*ybar2*ybar2)/(nn+a);
   double V21 = 0.0 + s21 + (nn*a*ybar1*ybar2)/(nn+a);
   //printf("a1,nu1,mu11,mu21: %lf,%lf,%lf,%lf\n",a1,nu1,mu11,mu21);
   //printf("V11,V21,V22: %lf,%lf,%lf\n",V11,V21,V22);

   //--------------------------------------------------
   //draw from posterior
   // draw G and then get Sigma and the chol of Sigma
   // for G, need the chol of V^{-1} 
   //V inverse
   double dt = V11*V22 - V21*V21; //determinant of V
   double iV11,iV22,iV21;
   iV11 = V22/dt;
   iV22 = V11/dt;
   iV21 = -V21/dt;
   //printf("Vinv: %lf, %lf, %lf\n",iV11,iV21,iV22);
   //chol of V inverse
   double L11,L21,L22;
   L11 = sqrt(iV11);
   L21 = iV21/L11;
   L22 = sqrt(iV22 - L21*L21);
   //printf("chol of Vinv: %lf, %lf, %lf\n",L11,L21,L22);

   //virtual double normal() = 0; //standard normal
   //virtual double chi_square() = 0; //chi-square
   //virtual void set_df(int df) = 0; //set df for chi-square
   //virtual void set_alpha(double alpha) = 0; //set df for chi-square

   //A : choleski root of G = Sigma^{-1}
   double Z = gen.normal();
   gen.set_df(nu1);
   double chi1 = sqrt(gen.chi_square());
   gen.set_df(nu1-1);
   double chi2 =  sqrt(gen.chi_square());
   double A11,A21,A22;
   A11 = L11 * chi1;
   A21 =  L21 * chi1 + L22*Z;
   A22 = L22*chi2;

   //G 
   double G11 = A11*A11;
   double G21 = A21*A11;
   double G22 = A21*A21 + A22*A22;

   //Sigma
   double Sig11,Sig21,Sig22;
   dt = G11*G22 - G21*G21;
   Sig11 = G22/dt;
   Sig21 = -G21/dt;
   Sig22 = G11/dt;
   //printf("Sigma 11,21,22: %lf, %lf %lf\n",Sig11,Sig21,Sig22);

   //LS : chol of Sigma
   double LS11 = sqrt(Sig11);
   double LS21 = Sig21/LS11;
   double LS22 = sqrt(Sig22 - LS21*LS21);
   //printf("chol of Sigma 11,21,22: %lf, %lf %lf\n",LS11,LS21,LS22);

   //draw mu
   double Z1 = gen.normal();
   double Z2 = gen.normal();
   double drmu1 = LS11*Z1;
   double drmu2 = LS21*Z1 + LS22*Z2;
   drmu1 = mu11 + sqrt(1/a1) * drmu1;
   drmu2 = mu21 + sqrt(1/a1) * drmu2;
   //printf("drmu1, drmu2: %lf, %lf\n",drmu1,drmu2);

   //fill in rvtheta with draw
   rtheta[0] = drmu1;
   rtheta[1] = drmu2;
   rtheta[2] = LS11;
   rtheta[3] = LS21;
   rtheta[4] = LS22;

   //cout << "&&&&&&&&&&&&&& out draw_one_theta\n";
   //printdv(rtheta);
   return rtheta;
}
//################################################################################
//### f(y,theta), theta=(mu,tau), no eta
double DpMuSigma::f(double* y, Dp::dv& theta, double* eta)
{
   //cout << "in f, theta is\n";
   //printdv(theta);
   //cout << "in f, y is \n";
   //cout << y[0] << ", " << y[1] << std::endl;

   // recall: theta = (mu1,mu2,L11,L21,L22) where L is chol(Sigma)
   // extract names params from theta
   double m1 = theta[0];
   double m2 = theta[1];
   double L11 = theta[2];
   double L21 = theta[3];
   double L22 = theta[4];
   //printf("L: %lf, %lf, %lf\n",L11,L21,L22);

   //inverse of L = IL
   double IL11 = 1.0/L11;
   double IL21 = -L21*IL11/L22;
   double IL22 = 1.0/L22;
   //printf("IL: %lf, %lf, %lf\n",IL11,IL21,IL22);

   double retval = -2.0*log(RTPI);  // log(1/2*pi)
   //cout << "log(p1): " << retval << std::endl;
   retval -= log(L11); // det(Sigma)^{-1/2}
   retval -= log(L22);
   //cout << "log(p2): " << retval << std::endl;

   double r1 = IL11*(y[0] - m1);
   double r2 = IL21*(y[0] - m1) + IL22*(y[1]-m2);
   retval += -.5*(r1*r1 + r2*r2);
   //cout << "log(p3): " << retval << std::endl;

   return(exp(retval));
}
//################################################################################
//### qo
double DpMuSigma::qo(double* y, double* eta)
{
   double retval = -2.0*log(RTPI);  // log(1/2*pi)
   double temp = (a*(nu-1))/((a+1)*v);
   retval += log(temp);
   temp = (y[0]*y[0]+y[1]*y[1])/v;
   temp *= a/(a+1.0);
   temp += 1.0;
   retval -= .5*(nu+1)*log(temp);

   return exp(retval);
}
//################################################################################
//### center
void DpMuSigma::center() 
{
   double m0=0.0;
   double m1=0.0;
   double p=0.0;
   //compute means
   for(auto c = theta.begin(); c != theta.end(); c++) {
      p = c->ind.size();
      p /= n;
      //cout << "p: " << p << std::endl;
      m0 += p*(c->theta[0]);
      m1 += p*(c->theta[1]);
   }
   //cout << "m0,m1: " << m0 << ", " << m1 << std::endl;
   //subtract off means
   for(auto c = theta.begin(); c != theta.end(); c++) {
      (c->theta[0]) -=m0;;
      (c->theta[1]) -=m1;;
   }
}


#endif
