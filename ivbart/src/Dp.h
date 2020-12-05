#ifndef GUARD_Dp
#define GUARD_Dp

#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>

#include "rrn.h"

//--------------------------------------------------------------------------------
// "theta cluster" class, has a theta value and indices of which obs use this theta
struct thetaGP {
public:
   //--------------------------------------------------
   //typedef
   //index iterator
   typedef std::list<size_t>::iterator iiter;

   //--------------------------------------------------
   //constructors/dest
   thetaGP():theta(),ind() {}
   thetaGP(const std::vector<double>& _theta):theta(_theta),ind() {}
   thetaGP(const std::vector<double>& _theta, const std::list<size_t>& _ind):theta(_theta),ind(_ind) {}
   ~thetaGP() {}

   //--------------------------------------------------
   // public functions
   void toscreen() {
      int p = theta.size(); 
      cout << "\n**********************\n";
      cout << "thetaGP object:\n" << "\tcluster size: " << ind.size() << 
            "\n\ttheta dim: "  << p << std::endl;
      if(p) cout << "\tfirst and last theta: " << theta[0] << ", " << theta[p-1] << std::endl;
   }

   //--------------------------------------------------
   //data members
   std::vector<double> theta; //the "theta" value
   std::list<size_t> ind;  //which observations have this theta value
};

class Dp {

public:
   //typedef
   //iterate over theta clusers (thetaGP instances)
   typedef std::list<thetaGP>::iterator titer;
   typedef std::vector<double> dv;
   typedef std::vector<dv> dvv;

   //const-dest
   //note: you either use (q,eta) for eta_i or etaconst for constant eta
   //   y is pxn
   Dp():theta(),n(0),p(0),y(nullptr),alpha(1.0),q(0),eta(nullptr),etaconst(),ag(),priag(),lag(),lgr() {}
   Dp(size_t _n, dv& _theta);
   virtual ~Dp() {}

   //public functions
   size_t npart() {return theta.size();}
   void setalpha(double _alpha) { alpha = _alpha;}
   double getalpha() { return alpha;}
   virtual void toscreen();
   dvv thetaMatrix();
   dv thetaVector();
   dvv thetastar();
   bool check();
   titer findTheta(dv& thetaval);
   void set(size_t vind, dv& v);
   std::vector<size_t> counts();
   void setData(size_t _p, double *_y) { p = _p; y = _y;} //assumes y has length p*n
   //note: three "eta states"
   //  (i) q=0, etaconst.size()=0 : no eta
   //  (ii) q=0, etaconst.size()!=0 : constant eta
   //  (iii) q=!0, eta!=nullptr, then eta has to point to q*n doubles
   void setetaconst(dv& _eta) {q=0; etaconst=_eta;}
   void seteta(size_t _q, double* _eta) {q=_q; eta=_eta;}
   void setAlphaPrior(dv& _ag, dv& _priag) {ag=_ag; priag=_priag; setAg();} //check them? at least lengths?

   //the key virtuals
   virtual dv draw_one_theta(std::list<size_t>& ind, rn& gen); //theta | y_i in ind
   //note that the "y and eta" in these calls will differ by the "y and eta" of the class
   //  in that they will refer to a single observation, e.g. y+i*p for observation i.
   virtual double f(double* y, dv& theta, double* eta); //f(y|theta,eta)
   virtual double qo(double* y, double* eta); // \int f(y|theta,eta) p(theta) dtheta

   //key posterior draws
   void drawTheta(rn& gen); //the central Escobar algorithm
   void remix(rn& gen);
   void drawAlpha(size_t k, rn& gen);
   void draw(rn& gen) { drawTheta(gen); remix(gen); drawAlpha(npart(),gen);}

protected:
   //parameter
   std::list<thetaGP> theta;
   //data
   size_t n;
   size_t p;
   double *y;                  //observations, should be pxn
   //alpha
   double alpha;

   //--------------------
   // eta, the other parameter taken as given f(y|theta,eta)
   // if q=0 (or eta==nullptr) then use etaconst
   size_t q; //dimension of eta parameter, if 0 then we are using etaconst
   double *eta; //if q!=0 should be q*n eta values, each block of q is an eta for an obsn.
   dv etaconst;  //if eta not set just plug in etaconst

   //---------------------
   // alpha prior
   dv ag; // discrete possible values for alpha
   dv priag; // prior on grid

   //---------------------
   //working 
   //for alpha prior draw on grid
   dv lag; // log(alpha) for alpha in ag
   dv lgr; // log(Gamma(a)/Gamma(n+a)) for a in ag
   void setAg();
};

void printdvv(Dp::dvv& m);
void printdv(Dp::dv& v);
size_t drind(double sum, Dp::dv& pv, rn& gen);
double logam(double x);


//#include <Dp.h>

//###--------------------------------------------------
//### Dp to screen
void Dp::toscreen()
{
   cout << "****************************\n";
   cout << "*****Dp object:\n";
   cout << "***alpha: " << alpha << std::endl;
   cout << "***alpha prior grid size: " << ag.size() << std::endl;
   if(ag.size()) {
      cout << "\tgrid first and last: " << ag[0] << ", " << *(--ag.end()) << std::endl;
      cout << "\tprior on grid first and last: " << *(priag.begin()) << ", " << *(--priag.end()) << std::endl;
   }
   cout << "***Number of Clusters: " << theta.size() << std::endl;
   if(theta.size()) {
      auto firstpart = theta.begin();
      cout << "***Dimension of theta: " << firstpart->theta.size() << std::endl;      
   } else {
      cout << "***Partition not set\n";
   }
   cout << "***eta:\n";
   if(q) {
      cout << "\t eta dimension (q): " << q << std::endl;
      if(n)
         cout << "\tfirst and last eta: " << eta[0] << ", " << eta[q*n-1] << std::endl;
   } else {
      cout << "\tUsing etaconst (q=0)\n";
      size_t ecd = etaconst.size();
      cout << "\tetaconst dimension: " << ecd << std::endl;
      if(ecd) cout << "\tfirst and last etaconst: " << etaconst[0] << ", " << etaconst[ecd-1] << std::endl;
   }
   if(y==nullptr) {
      cout << "***data not set\n";
   } else {
      cout << "***data:\n";
      cout << "\tn: " << n << std::endl;
      cout << "\tp: " << p << std::endl;
      cout << "\tfirst and last y: " << y[0] << ", " << y[p*n-1] << std::endl;
   }
   if(theta.size()) {
      cout << "***Clusters: " << std::endl;
      int clid=0; //cluster id
      for(titer i=theta.begin();i!=theta.end();i++) { 
         cout << "\tOn cluster: " << clid << ", size: " << i->ind.size() << std::endl;
         //printdv(i->theta);
         cout << "\t\tCluster theta: ";
         for(auto ii=i->theta.begin();ii!=i->theta.end();ii++) cout << *ii << ", ";
         cout << std::endl;
         clid+=1;
      }
   }
}
//inititialize to one cluster with the input theta value
Dp::Dp(size_t _n, std::vector<double>& _theta) : theta(),n(0),p(0),y(nullptr),alpha(1.0),q(0),eta(nullptr),etaconst(),ag(),priag(),lag(),lgr()
{
   n = _n;

   //make a single theta cluster
   thetaGP th;
   th.theta = _theta;
   th.ind.clear();

   //assign all observations to theta cluster 
   for(size_t i=0;i<n;i++) {
      th.ind.push_back(i);
   }

   //add the single cluster to the theta list 
   theta.push_back(th);
}
//###--------------------------------------------------
//### get matrix with value for each observation (so some may be repeated)
Dp::dvv Dp::thetaMatrix() 
{
   Dp::dvv m(n);
   for(titer i=theta.begin(); i!=theta.end();i++) {
      for(thetaGP::iiter j=i->ind.begin();j!=i->ind.end();j++) m[*j] = i->theta;
   }
   return m;
}
//###--------------------------------------------------
//### get vector of stacked that values, one for each observation, vec of thetaMatrix
Dp::dv Dp::thetaVector() 
{
   auto ii = theta.begin();
   size_t p = ii->theta.size();
   //cout << "in thetaVector, p is: " << p << std::endl;
   //cout << "in thetaVector, n is: " << n << std::endl;
   Dp::dv v(n*p);
   for(titer ip=theta.begin(); ip!=theta.end();ip++) { // loop over partitions
      for(thetaGP::iiter j=ip->ind.begin();j!=ip->ind.end();j++)  { //loop over obs in partition
         size_t obsn = *j;
         for(size_t k=0;k<p;k++) { //loop over theta elements for theta of partition
            v[obsn*p+k]=ip->theta[k];
         }
      }
   }
   return v; 
}
//###--------------------------------------------------
//### get theta star vector 
Dp::dvv Dp::thetastar() 
{
   Dp::dvv tstar;
   for(titer i= theta.begin();i!=theta.end();i++) tstar.push_back(i->theta);
   return tstar;
}
//###--------------------------------------------------
//###check that the thetastar are unique, and the the union of indices is 0,1,2,,n-1
bool Dp::check() 
{
   Dp::dvv tst = thetastar();
   //std::sort(tst.begin(),tst.end());
   //for(std::vector<double>::size_type i=1;i<tst.size();i++) if(tst[i]==tst[i-1]) return false;
   size_t nts = tst.size();
   cout << "in check nts (num thetastar): " << nts << std::endl;
   for(size_t i=0;i<nts;i++) {
      Dp::dv curts = tst[i];
      for(size_t j=0;j<i;j++) {
         if (tst[j]==curts) return false;
      }
      for(size_t j=(i+1);j<nts;j++) {
         if (tst[j]==curts) return false;
      }
   }

   std::list<size_t> allind;
   for(titer i=theta.begin();i!=theta.end();i++) copy(i->ind.begin(),i->ind.end(),back_inserter(allind));
   allind.sort();

   size_t cnt=0;
   for(thetaGP::iiter i=allind.begin();i!=allind.end();i++) {
      if(*i != cnt) return false;
      cout << "checking observation idex: " << *i << std::endl;
      cnt++;
   }
   return true;
}
//###--------------------------------------------------
//### find iterator of cluster with theta = val
Dp::titer Dp::findTheta(Dp::dv& thetaval)
{
   titer i = theta.begin();
   bool not_found = true;
   while(not_found && i!=theta.end()) {
      if(i->theta == thetaval) not_found = false;
      else i++;
   }
   //if(i == theta.end()) cout << "\n*********theta not found in Dp::findTheta\n";
   //else cout << "\n**** theta found in Dp::findTheta\n";
   return i;
}
//###--------------------------------------------------
//### set the theta value for a particular observation
// ### set theta[vind] to v
void Dp::set(size_t vind, Dp::dv& v)
{
   //should check that n>0
   cout << "\n\n**********************into set\n\n";

   // first get a pointer the cluster, and a pointer to the index in that cluster's ind
   titer i  = theta.begin(); //which cluster is vind in, has to be in one of them
   thetaGP::iiter j; //which index in the cluster is vind
   bool not_found = true;
   while(not_found && i != theta.end()) {
       j = find(i->ind.begin(),i->ind.end(),vind);
       if(j != i->ind.end())
          not_found = false;
        else
          i++;
   }
   if(not_found) {cout << "error: index not found in Dp::set\n"; 
     //exit(EXIT_FAILURE);
   }
   if(!not_found) cout << "index found in set check: " << *j << ", " << vind << std::endl;

   //now do the updating
   if(i->theta != v) { //if it already has the value, nothing to do
      i->ind.erase(j); // drop it from old cluster
      if(i->ind.size() == 0) theta.erase(i); //if cluster now empty, drop it from theta
      titer k = findTheta(v); //see if you can find a theta with value v
      if(k == theta.end()) { //if you can't make a new cluster and add it to theta
         thetaGP tg;
         tg.theta = v;
         tg.ind.push_back(vind);
         theta.push_back(tg);
      } else { //if you can, just add the index to the list for the cluster
         k->ind.push_back(vind);
      }
   }
}
//###--------------------------------------------------
//### key virtual, f
double Dp::f(double* y, dv& theta, double* eta)
{
   return 1.0;
}
//###--------------------------------------------------
//### key virtual, draw_one_theta
Dp::dv Dp::draw_one_theta(std::list<size_t>& ind, rn& gen)
{
   size_t k = theta.begin()->theta.size();
   Dp::dv ret;
   for(size_t i=0;i<k;i++) ret.push_back(gen.uniform());
   return ret;
}
//###--------------------------------------------------
//### key virtual, qo
double Dp::qo(double* y, double* eta)
{
   return 1.0;
}
//###--------------------------------------------------
//### get the size of each theta cluster
std::vector<size_t> Dp::counts()
{
   std::vector<size_t> cnts;
   for(titer i= theta.begin();i!=theta.end();i++) cnts.push_back(i->ind.size());
   return cnts;
}
//###--------------------------------------------------
//###  basic DP draw of \theta_i | the rest
void Dp::drawTheta(rn& gen)
{
   //cout << "&&&& Into drawTheta\n\n";

    //get counts and thetastar for initial partition
   Dp::dvv tst = thetastar();
   std::vector<size_t> cnt = counts();
   std::vector<size_t> icnt = cnt; //need to know the initial number
   size_t np = cnt.size(); //initial number of partitions or clusters

   /*
   cout << "size of tst: " << tst.size() << std::endl;
   printdvv(tst);
   cout << "cnt:\n";
   for(auto i=cnt.begin();i!=cnt.end();i++) cout << "count: " << *i << ", ";
   cout << "\n\n";
   */

   //The idea is that we loop over all theta_i (all observations) by looping over the clusters and then
   //   the observations in a cluster.
   //We update the custers as we go.  If we get a new one we put it on the end of the theta list.
   //  After we are done we delete any empties.

   //this is tricky code because dimension changes, new partitions get put on the end
   //i loops over the clusters and j loops over observations in a cluster
   //  but, since these change as we go along we mangage:
   //   par: pointer to current cluster in updated partition
   //   ob: pointer to current observation index with in a cluster

   //cout << "\n\n&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& Start Loop:\n\n";
   titer par = theta.begin();
   for(size_t i=0;i<np;i++) {
      thetaGP::iiter ob = par->ind.begin();
      for(size_t j=0;j<icnt[i];j++) { 
         size_t nts = tst.size();
         std::vector<double> pv(nts + 1);
         double sum=0.0;
         for(size_t k=0;k<nts;k++) {
            if(cnt[k] == 0) { //tricky, because you drop empties at the end
               pv[k] = 0.0;
            } else {
               size_t thecnt = cnt[k];
               if(tst[k]==par->theta) thecnt -=1; //par->theta will be the theta of the obs you are on
               if(eta == nullptr) {
                  pv[k] = thecnt*f(y+(*ob)*p,tst[k],&etaconst[0]);
               } else {
                  pv[k] = thecnt*f(y+(*ob)*p,tst[k],eta+(*ob)*q);
               }
            }
            sum += pv[k];
         }
         if(eta == nullptr) {
            pv[nts] = alpha * qo(y+(*ob)*p,&etaconst[0]);
         } else {
            pv[nts] = alpha * qo(y+(*ob)*p,eta+(*ob)*q);
         }
         sum += pv[nts];
         size_t ii = drind(sum,pv,gen);

         //cout << "\n****&&&& in drawTheta, ii: " << ii << std::endl;
         //for(size_t jj = 0;jj<nts;jj++) cout << "i,j,jj: "
          //             << i << " " << j << " " << jj << std::endl;
   /*
   cout << "size of tst: " << tst.size() << std::endl;
   cout << "thetastar\n";
   printdvv(tst);
   cout << "cnt:\n";
   for(auto i=cnt.begin();i!=cnt.end();i++) cout << *i << ", ";
   cout << "\npv\n";
   printdv(pv);
   */

         if(ii==nts) {  // birth
            //draw new theta
            std::list<size_t> tind; tind.push_back(*ob);
            Dp::dv ntheta = draw_one_theta(tind,gen); 

            //update tst and cnt
            tst.push_back(ntheta); cnt.push_back(1); cnt[i] -=1;

            //add new cluster (of size 1) to theta
            thetaGP ngp; ngp.theta = ntheta; ngp.ind.push_back(*ob);
            theta.push_back(ngp);

            //manage ob
            ob = par->ind.erase(ob);
         } else if(tst[ii] == par->theta) { //drew yourself, just go on to next observation
            ob++;
         } else {  //drew one of the other old thetas
            //update cnt
            cnt[i] -=1; cnt[ii] += 1;
            // add to other old
            titer tp = findTheta(tst[ii]); tp->ind.push_back(*ob);
            // delete from current old
            ob = par->ind.erase(ob);
         }
      }
      par++;
   }
   //drop all empty partitions
   par = theta.begin();
   while(par!=theta.end()) {
      if(par->ind.size() == 0) par = theta.erase(par);
      else  par++;
   }

}
//###--------------------------------------------------
//### Dp::remix, redraw theta in each cluster.
void Dp::remix(rn& gen)
{
   for(titer i = theta.begin(); i!= theta.end(); i++) {
      i->theta = draw_one_theta(i->ind,gen);
   }
}
//###--------------------------------------------------
//### Dp::setAg, set up stuff for p(a|k) which does not depend on k
void Dp::setAg()
{
   size_t ng = ag.size();  //size of alpha grid
   lag.resize(ng);
   lgr.resize(ng);

   //the working vectors lag and lgr store precomputed parts of p(a|k) 
   //   which depend on alpha but not k= num unique theta
   for(size_t i=0;i<ng;i++) {
      lag[i] = log(ag[i]);
      lgr[i] = logam(ag[i]) - logam(ag[i]+n);
   }
}
//###--------------------------------------------------
void Dp::drawAlpha(size_t k, rn& gen)
{
   size_t ng = ag.size();  //size of alpha grid
   Dp::dv postv(ng); //log post then post
   double maxll = - std::numeric_limits<double>::infinity();
   for(size_t i=0;i<ng;i++) {
      postv[i] = k*lag[i] + lgr[i] + log(priag[i]);
      if(postv[i] > maxll) maxll = postv[i];
   }

   double sum=0.0;
   for(size_t i=0;i<ng;i++) {
      postv[i] =  exp(postv[i]-maxll);
      sum += postv[i];
   }
   for(size_t i=0;i<ng;i++) postv[i] /= sum;

   size_t ii = drind(1.0,postv,gen);
   alpha = ag[ii];
}

//###--------------------------------------------------
//###--------------------------------------------------
//### utility functions
void printdvv(Dp::dvv& m) {
   size_t nr = m.size();
   for(int i=0;i<nr;i++) {
      int p = m[i].size();
      for(int j=0;j<p;j++) cout << m[i][j] << " ";
      cout << std::endl;
   }

}
void printdv(Dp::dv& v) {
   size_t nv = v.size();
   for(size_t i=0;i<nv;i++) cout << v[i] << std::endl;
}
//###--------------------------------------------------
//### draw from a discrete, note that sum is sum of pv, so pv may be unnormalized
size_t drind(double sum,Dp::dv& pv, rn& gen) 
{
   double u = gen.uniform();
   double psum=pv[0]/sum;
   size_t ii=0;
   while(psum<u) {
      ii++;
      psum += pv[ii]/sum;
   }
   return ii;
}

//--------------------------------------------------
// compute the logarithm of the Gamma function
// people.sc.fsu.edu/~jburkardt/cpp_src/toms291/toms291.html
double logam(double x)
{
  double f;
  double value;
  double y;
  double z;

  if ( x <= 0.0 )
  {
    value = 0.0;
    return value;
  }

  y = x;

  if ( x < 7.0 )
  {
    f = 1.0;
    z = y;

    while ( z < 7.0 )
    {
      f = f * z;
      z = z + 1.0;
    }
    y = z;
    f = - log ( f );
  }
  else
  {
    f = 0.0;
  }

  z = 1.0 / y / y;

  value = f + ( y - 0.5 ) * log ( y ) - y 
    + 0.918938533204673 + 
    ((( 
    - 0.000595238095238   * z 
    + 0.000793650793651 ) * z 
    - 0.002777777777778 ) * z 
    + 0.083333333333333 ) / y;

  return value;
}

#endif

