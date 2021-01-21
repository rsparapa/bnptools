//     brt.h: Base BART model class definition.
//     Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
//
//     This file is part of BART.
//
//     BART is free software: you can redistribute it and/or modify
//     it under the terms of the GNU Affero General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
//
//     BART is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU Affero General Public License for more details.
//
//     You should have received a copy of the GNU Affero General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//     Author contact information
//     Matthew T. Pratola: mpratola@gmail.com
//     Robert E. McCulloch: robert.e.mculloch@gmail.com
//     Hugh A. Chipman: hughchipman@gmail.com


#ifndef GUARD_brt_h
#define GUARD_brt_h

/*
#include "tree.h"
#include "treefuns.h"
#include "dinfo.h"

#ifdef _OPENMP
#   include <omp.h>
#endif

#include "brt.h"
#include "brtfuns.h"
#include <iostream>
#include <map>
#include <vector>

#ifndef NotInR
#include <Rcpp.h>
#define COUT Rcpp::Rcout 
#else
#define COUT std::cout  << "std::cout "
#endif

using std::cout;
using std::endl;
*/

class sinfo { //sufficient statistics (will depend on end node model)
public:
   sinfo(): n(0) {}
   sinfo(const sinfo& is):n(is.n) {}
   virtual ~sinfo() {}  //need this so memory is properly freed in derived classes.

   size_t n;
   // compound addition operator needed when adding suff stats
   virtual sinfo& operator+=(const sinfo& rhs) {
      n=n+rhs.n;
      return *this;
   }
   // assignment operator for suff stats
   virtual sinfo& operator=(const sinfo& rhs)
   {
      if(&rhs != this) {
         this->n = rhs.n;
      }
      return *this; 
   }
   // addition opertor is defined in terms of compound addition
   const sinfo operator+(const sinfo& other) const {
      sinfo result = *this; //copy of myself.
      result += other;
      return result;
   }
};

class brt {
public:
   //--------------------
   //classes
   class tprior { //prior on the tree structure
   public:
      tprior(): alpha(.95),beta(1.0) {}
      //prob(split) = alpha/(1+d)^beta, where d is node depth.
      double alpha;
      double beta;
   };
   class mcmcinfo { //algorithm parameters (eg move probabilities)
   public:
      mcmcinfo(): pbd(1.0),pb(.5),minperbot(5),dopert(true),pertalpha(0.1),pertproposal(1),
                  pertaccept(0),rotproposal(0),rotaccept(1),bproposal(0),baccept(1),dproposal(0),
                  daccept(1),pchgv(0.2),chgvproposal(1),chgvaccept(0),corv(0),
                  dostats(false),varcount(0),tavgd(0.0),tmaxd(0) {} 
      double pbd;
      double pb;
      size_t minperbot;
      bool dopert;
      double pertalpha;
      size_t pertproposal;  //number of perturb proposals
      size_t pertaccept;    //number of accepted perturb proposals
      size_t rotproposal;   //number of rotate proposals
      size_t rotaccept;     //number of accepted rotate proposals
      size_t bproposal;     //number of birth proposals
      size_t baccept;       //number of birth accepts
      size_t dproposal;     //number of death proposals
      size_t daccept;       //number of death accepts
      double pchgv;         //probability of change of variable proposal.  Probability of perturb proposal is 1-pchgv
      size_t chgvproposal;  //number of change of variable proposals
      size_t chgvaccept;    //number of accepted chnage of varialbe proposals
      std::vector<std::vector<double> >* corv; //initial proposal distribution for changing variables in pert
      //statistics
      bool dostats;    //keep track of statistics yes or no?
      unsigned int* varcount;//count number of splits in tree on each variable
      double tavgd;         //average tree depth
      unsigned int tmaxd;   //maximum tree depth
      unsigned int tmind;   //minimum tree depth
   };
   class cinfo { //parameters for end node model prior
   };
   //--------------------
   //constructors/destructors
   brt():t(0.0),tp(),xi(0),ci(),di(0),mi(),tc(1) {}
   virtual ~brt() { if(mi.varcount) delete[] mi.varcount; }
   //--------------------
   //methods
   void settc(int tc) {this->tc = tc;}
   void setxi(xinfo *xi) {this->xi=xi; this->ncp1=2.0;
                          for(size_t i=0;i<(*xi).size();i++) 
                            if(this->ncp1<(double)((*xi)[i].size()+1.0))
                              this->ncp1=(double)((*xi)[i].size()+1.0);
                         }
   void setdata(dinfo *di) {this->di=di; resid.resize(di->n); yhat.resize(di->n); setf(); setr(); }
   void pr();
   void settp(double alpha, double beta) {tp.alpha=alpha;tp.beta=beta;}
   void setmi(double pbd, double pb, size_t minperbot, bool dopert, double pertalpha, double pchgv, std::vector<std::vector<double> >* chgv)
             {mi.pbd=pbd; mi.pb=pb; mi.minperbot=minperbot; mi.dopert=dopert;
              mi.pertalpha=pertalpha; mi.pchgv=pchgv; mi.corv=chgv; }
   void setstats(bool dostats) { mi.dostats=dostats; if(dostats) mi.varcount=new unsigned int[xi->size()]; }
   void getstats(unsigned int* vc, double* tad, unsigned int* tmd, unsigned int* tid) { *tad=mi.tavgd; *tmd=mi.tmaxd; *tid=mi.tmind; for(size_t i=0;i<xi->size();i++) vc[i]=mi.varcount[i]; }
   void addstats(unsigned int* vc, double* tad, unsigned int* tmd, unsigned int* tid) { *tad+=mi.tavgd; *tmd=std::max(*tmd,mi.tmaxd); *tid=std::min(*tid,mi.tmind); for(size_t i=0;i<xi->size();i++) vc[i]+=mi.varcount[i]; }
   void resetstats() { mi.tavgd=0.0; mi.tmaxd=0; mi.tmind=0; for(size_t i=0;i<xi->size();i++) mi.varcount[i]=0; }
   void setci() {}
   void draw(rn& gen);
   virtual sinfo* newsinfo() { return new sinfo; }
   virtual std::vector<sinfo*>& newsinfovec() { std::vector<sinfo*>* si= new std::vector<sinfo*>; return *si; }
   virtual std::vector<sinfo*>& newsinfovec(size_t dim) { std::vector<sinfo*>* si = new std::vector<sinfo*>; si->resize(dim); for(size_t i=0;i<dim;i++) si->push_back(new sinfo); return *si; }
   void setf();  //set the vector of predicted values
   void setr();  //set the vector of residuals
   double f(size_t i) { return yhat[i]; }  //get the i'th predicted value
   double r(size_t i) { return resid[i]; }  //get the i'th residual
   std::vector<double>* getf() { return &yhat; }
   std::vector<double>* getr() { return &resid; }
   void predict(dinfo* dipred); // predict y at the (npred x p) settings *di.x
//   void savetree(int* id, int* v, int* c, double* theta);  //save tree to vector output format
   void savetree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
//   void loadtree(size_t nn, int* id, int* v, int* c, double* theta);  //load tree from vector input format
   void loadtree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
   //--------------------
   //data
   tree t;
   //--------------------------------------------------
   //stuff that maybe should be protected
   void bd(rn& gen);      //uses getsuff
   void pertcv(rn& gen);  //uses getpertsuff which in turn uses subsuff
   void drawtheta(rn& gen);
   void allsuff(tree::npv& bnv,std::vector<sinfo*>& siv);  //assumes brt.t is the root node
   void subsuff(tree::tree_p nx, tree::npv& bnv, std::vector<sinfo*>& siv); //does NOT assume brt.t is the root node.
                                                                           //Instead, uses the path from nx that it constructs.
   bool rot(tree::tree_p tnew, tree& x, rn& gen);  //uses subsuff
   void adapt();
protected:
   //--------------------
   //model information
   tprior tp; //prior on tree (alpha and beta)
   xinfo *xi; //cutpoints for each x.
   double ncp1; //number of cutpoints for each x plus 1.
   cinfo ci; //conditioning info (e.g. other parameters and prior and end node models)
   //--------------------
   //data
   dinfo *di; //n,p,x,y
   std::vector<double> yhat; //the predicted vector
   std::vector<double> resid; //the actual residual vector
   //--------------------
   //mcmc info
   mcmcinfo mi;
   //thread count
   int tc;
   //number of trees -- is always 1 in the base class.
//   size_t m;
   //--------------------
   //methods
   virtual void add_observation_to_suff(diterator& diter, sinfo& si); //add in observation i (from di) into si (possibly using ci)
   void getsuff(tree::tree_p nx, size_t v, size_t c, sinfo& sil, sinfo& sir);  //assumes brt.t is the root node
   void getsuff(tree::tree_p l, tree::tree_p r, sinfo& sil, sinfo& sir);       //assumes brt.t is the root node
   void getchgvsuff(tree::tree_p pertnode, tree::npv& bnv, size_t oldc, size_t oldv, bool didswap, 
                  std::vector<sinfo*>& sivold, std::vector<sinfo*>& sivnew);     //uses subsuff
   void getpertsuff(tree::tree_p pertnode, tree::npv& bnv, size_t oldc,        //uses subsuff
                  std::vector<sinfo*>& sivold, std::vector<sinfo*>& sivnew);
   void local_getsuff(diterator& diter, tree::tree_p nx, size_t v, size_t c, sinfo& sil, sinfo& sir); 
   void local_getsuff(diterator& diter, tree::tree_p l, tree::tree_p r, sinfo& sil, sinfo& sir);
   virtual double lm(sinfo& si); //uses pi. 
   virtual double drawnodetheta(sinfo& si, rn& gen);
   void local_allsuff(diterator& diter, tree::npv& bnv,std::vector<sinfo*>& siv);
   void local_subsuff(diterator& diter, tree::tree_p nx, tree::npv& path, tree::npv& bnv, std::vector<sinfo*>& siv);
   virtual void local_setf(diterator& diter);
   virtual void local_setr(diterator& diter);
   virtual void local_predict(diterator& diter);
//#  ifdef _OPENMP
   void local_ompgetsuff(tree::tree_p nx, size_t v, size_t c, dinfo di, sinfo& sil, sinfo& sir);
   void local_ompgetsuff(tree::tree_p l, tree::tree_p r, dinfo di, sinfo& sil, sinfo& sir);
   void local_ompallsuff(dinfo di, tree::npv bnv,std::vector<sinfo*>& siv);
   void local_ompsubsuff(dinfo di, tree::tree_p nx, tree::npv& path, tree::npv bnv,std::vector<sinfo*>& siv);
   void local_ompsetf(dinfo di);
   void local_ompsetr(dinfo di);
   void local_omppredict(dinfo dipred);
//   void local_ompsavetree(int* id, int* v, int* c, double* theta);
//   void local_savetree(int* id, int* v, int* c, double* theta);
//   void local_omploadtree(size_t nn, int* id, int* v, int* c, double* theta);
//   void local_loadtree(size_t nn, int* id, int* v, int* c, double* theta);
   void local_ompsavetree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
   virtual void local_savetree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
   void local_omploadtree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
   virtual void local_loadtree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
//#  endif
};


//--------------------------------------------------
//a single iteration of the MCMC for brt model
void brt::draw(rn& gen)
{
   // Structural/topological proposal(s)
   if(gen.uniform()<mi.pbd)
      bd(gen);
   else
   {
      tree::tree_p tnew;
      tnew=new tree(t); //copy of current to make life easier upon rejection
      rot(tnew,t,gen);
      delete tnew;
   }

   // Perturbation Proposal
   if(mi.dopert)
      pertcv(gen);
   // Gibbs Step
    drawtheta(gen);

   //update statistics
   if(mi.dostats) {
      tree::npv bnv; //all the bottom nodes
      for(size_t k=0;k< xi->size();k++) mi.varcount[k]+=t.nuse(k);
      t.getbots(bnv);
      //unsigned int tempdepth[bnv.size()];
      std::vector<unsigned int> tempdepth(bnv.size());
      unsigned int tempavgdepth=0;
      for(size_t i=0;i!=bnv.size();i++)
         tempdepth[i]=(unsigned int)bnv[i]->depth();
      for(size_t i=0;i!=bnv.size();i++) {
         tempavgdepth+=tempdepth[i];
         mi.tmaxd=std::max(mi.tmaxd,tempdepth[i]);
         mi.tmind=std::min(mi.tmind,tempdepth[i]);
      }
      mi.tavgd+=((double)tempavgdepth)/((double)bnv.size());
   }
}
//--------------------------------------------------
//adapt the proposal widths for perturb proposals,
//bd or rot proposals and b or d proposals.
void brt::adapt()
{
   //double pert_rate,b_rate,rot_rate,m_rate,chgv_rate;
   double pert_rate,m_rate,chgv_rate;

   pert_rate=((double)mi.pertaccept)/((double)mi.pertproposal);
   chgv_rate=((double)mi.chgvaccept)/((double)mi.chgvproposal);
//   pert_rate=((double)(mi.pertaccept+mi.baccept+mi.daccept+mi.rotaccept))/((double)(mi.pertproposal+mi.dproposal+mi.bproposal+mi.rotproposal));
//   b_rate=((double)mi.baccept)/((double)mi.bproposal);
//   d_rate=((double)mi.daccept)/((double)mi.dproposal);
//   bd_rate=((double)(mi.baccept+mi.daccept))/((double)(mi.dproposal+mi.bproposal));
//   rot_rate=((double)mi.rotaccept)/((double)mi.rotproposal);
   m_rate=((double)(mi.baccept+mi.daccept+mi.rotaccept))/((double)(mi.dproposal+mi.bproposal+mi.rotproposal));

   //update pbd
   // a mixture between calibrating to m_rate (25%) and not moving too quickly away from
   // the existing probability of birth/death (75%):
//   mi.pbd=0.25*mi.pbd*m_rate/0.24+0.75*mi.pbd;
   // avoid too small or large by truncating to 0.1,0.9 range:
//   mi.pbd=std::max(std::min(0.9,mi.pbd),0.1);

   //update pb
//old   mi.pb=mi.pb*bd_rate/0.24;
//old   mi.pb=mi.pb*(b_rate+d_rate)/2.0/bd_rate;
   // a mixture between calibrating to the (bd_rate and m_rate) and existing probability of birth
   // in other words, don't move too quickly away from existing probability of birth
   // and when we do move, generally we favor targetting bd_rate (90%) but also target m_rate to
   // a small degree (10%):
//   mi.pb=0.25*(0.9*mi.pb*(b_rate+d_rate)/2.0/bd_rate + 0.1*mi.pb*m_rate/0.24)+0.75*mi.pb;
   // avoid too small or large by truncating to 0.1,0.9 range:
//   mi.pb=std::max(std::min(0.9,mi.pb),0.1);

   //update pertalpha
   mi.pertalpha=mi.pertalpha*pert_rate/0.44;
//   if(mi.pertalpha>2.0) mi.pertalpha=2.0;
//   if(mi.pertalpha>(1.0-1.0/ncp1)) mi.pertalpha=(1.0-1.0/ncp1);
   if(mi.pertalpha>2.0) mi.pertalpha=2.0;
   if(mi.pertalpha<(1.0/ncp1)) mi.pertalpha=(1.0/ncp1);

   mi.pertaccept=0; mi.baccept=0; mi.rotaccept=0; mi.daccept=0;
   mi.pertproposal=1; mi.bproposal=1; mi.rotproposal=1; mi.dproposal=1;
   //if(mi.dostats) {
   COUT << "pert_rate=" << pert_rate << " pertalpha=" << mi.pertalpha << " chgv_rate=" << chgv_rate;
   // COUT << "   b_rate=" << b_rate << endl;
   // COUT << "   d_rate=" << d_rate << endl;
   // COUT << "   bd_rate=" << bd_rate << endl;
   // COUT << " rot_rate=" << rot_rate << endl;
   COUT << "   m_rate=" << m_rate;
   //   COUT << "mi.pbd=" << mi.pbd << "  mi.pb=" << mi.pb<< "  mi.pertalpha=" << mi.pertalpha << endl;
   //   COUT << endl;
   //}
}
//--------------------------------------------------
//draw all the bottom node theta's for the brt model
void brt::drawtheta(rn& gen)
{
   tree::npv bnv;
//   std::vector<sinfo> siv;
   std::vector<sinfo*>& siv = newsinfovec();

   allsuff(bnv,siv);
   for(size_t i=0;i<bnv.size();i++) {
      bnv[i]->settheta(drawnodetheta(*(siv[i]),gen));
      delete siv[i]; //set it, then forget it!
   }
   delete &siv;  //and then delete the vector of pointers.
}
//--------------------------------------------------
//draw theta for a single bottom node for the brt model
double brt::drawnodetheta(sinfo& si, rn& gen)
{
//   return 1.0;
   gen.uniform(); //dummy use of gen
   return si.n;
}
//--------------------------------------------------
//pr for brt
void brt::pr()
{
   COUT << "***** brt object:\n";
   if(xi) {
      size_t p = xi->size();
      COUT  << "**xi cutpoints set:\n";
      COUT << "\tnum x vars: " << p << endl;
      COUT << "\tfirst x cuts, first and last " << (*xi)[0][0] << ", ... ," << 
              (*xi)[0][(*xi)[0].size()-1] << endl;
      COUT << "\tlast x cuts, first and last " << (*xi)[p-1][0] << ", ... ," << 
              (*xi)[p-1][(*xi)[p-1].size()-1] << endl;
   } else {
      COUT << "**xi cutpoints not set\n";
   }
   if(di) {
      COUT << "**data set, n,p: " << di->n << ", " << di->p << endl;
   } else {
      COUT << "**data not set\n";
   }
   COUT << "**the tree:\n";
   t.pr();
}
//--------------------------------------------------
//lm: log of integrated likelihood, depends on prior and suff stats
double brt::lm(sinfo& si)
{
   COUT << "in brt::lm " << si.n << "\n"; // dummy use of si to get rid of unused argument warning
   return 0.0;  //just drawing from prior for now.
}
//--------------------------------------------------
//getsuff used for birth.
void brt::local_getsuff(diterator& diter, tree::tree_p nx, size_t v, size_t c, sinfo& sil, sinfo& sir)    
{
   double *xx;//current x
   sil.n=0; sir.n=0;

   for(;diter<diter.until();diter++)
   {
      xx = diter.getxp();
      if(nx==t.bn(diter.getxp(),*xi)) { //does the bottom node = xx's bottom node
         if(xx[v] < (*xi)[v][c]) {
               //sil.n +=1;
               add_observation_to_suff(diter,sil);
          } else {
               //sir.n +=1;
               add_observation_to_suff(diter,sir);
          }
      }
   }
}
//--------------------------------------------------
//getsuff used for death
void brt::local_getsuff(diterator& diter, tree::tree_p l, tree::tree_p r, sinfo& sil, sinfo& sir)
{
   sil.n=0; sir.n=0;

   for(;diter<diter.until();diter++)
   {
      tree::tree_cp bn = t.bn(diter.getxp(),*xi);
      if(bn==l) {
         //sil.n +=1;
         add_observation_to_suff(diter,sil);
      }
      if(bn==r) {
         //sir.n +=1;
         add_observation_to_suff(diter,sir);
      }
   }
}
//--------------------------------------------------
//Add in an observation, this has to be changed for every model.
//Note that this may well depend on information in brt with our leading example
//being double *sigma in cinfo for the case of e~N(0,sigma_i^2).
// Note that we are using the training data and the brt object knows the training data
//     so all we need to specify is the row of the data (argument size_t i).
void brt::add_observation_to_suff(diterator& diter, sinfo& si)
{
   COUT << "in brt::add_observation_to_suff, diter.gety() is " << diter.gety() << "\n";  // dummy use to get rid of warning
   si.n+=1; //in add_observation_to_suff
}
//--------------------------------------------------
//getsuff wrapper used for birth.  Calls serial or parallel code depending on how
//the code is compiled.
void brt::getsuff(tree::tree_p nx, size_t v, size_t c, sinfo& sil, sinfo& sir)
{
   #ifdef _OPENMP
#     pragma omp parallel num_threads(tc)
      local_ompgetsuff(nx,v,c,*di,sil,sir); //faster if pass dinfo by value.
   #else
         diterator diter(di);
         local_getsuff(diter,nx,v,c,sil,sir);
   #endif
}
//--------------------------------------------------
//allsuff (1)
void brt::allsuff(tree::npv& bnv,std::vector<sinfo*>& siv)
{
   #ifdef _OPENMP
   typedef tree::npv::size_type bvsz;
   #endif

   //get bots once and pass them around
   bnv.clear();
   t.getbots(bnv);

   #ifdef _OPENMP
      siv.clear(); //need to setup space threads will add into
      siv.resize(bnv.size());
      for(bvsz i=0;i!=bnv.size();i++) siv[i]=newsinfo();
#     pragma omp parallel num_threads(tc)
      local_ompallsuff(*di,bnv,siv); //faster if pass di and bnv by value.
   #else
         diterator diter(di);
         local_allsuff(diter,bnv,siv); //will resize siv
   #endif
}
//--------------------------------------------------
//local_subsuff
void brt::local_subsuff(diterator& diter, tree::tree_p nx, tree::npv& path, tree::npv& bnv, std::vector<sinfo*>& siv)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observation
   size_t ni;         //the  index into vector of the current bottom node
   size_t index;      //the index into the path vector.
   double *x;
   tree::tree_p root=path[path.size()-1];

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   siv.clear();
   siv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) { bnmap[bnv[i]]=i; siv[i]=newsinfo(); }

   for(;diter<diter.until();diter++) {
      index=path.size()-1;
      x=diter.getxp();
      if(root->xonpath(path,index,x,*xi)) { //x is on the subtree, 
         tbn = nx->bn(x,*xi);              //so get the right bn below interior node n.
         ni = bnmap[tbn];
         //siv[ni].n +=1;
         add_observation_to_suff(diter, *(siv[ni]));
      }
      //else this x doesn't map to the subtree so it's not added into suff stats.
   }
}
//--------------------------------------------------
//local_ompsubsuff
void brt::local_ompsubsuff(dinfo di, tree::tree_p nx, tree::npv& path, tree::npv bnv,std::vector<sinfo*>& siv)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = di.n;
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);

   std::vector<sinfo*>& tsiv = newsinfovec(); //will be sized in local_subsuff
   diterator diter(&di,beg,end);
   local_subsuff(diter,nx,path,bnv,tsiv);

#  pragma omp critical
   {
      for(size_t i=0;i<siv.size();i++) *(siv[i]) += *(tsiv[i]);
   }

   for(size_t i=0;i<tsiv.size();i++) delete tsiv[i];
   delete &tsiv;
#endif
}

//--------------------------------------------------
//get suff stats for bots that are only below node n.
//NOTE!  subsuff is the only method for computing suff stats that does not
//       assume the root of the tree you're interested is brt.t.  Instead,
//       it takes the root of the tree to be the last entry in the path
//       vector.  In other words, for MCMC proposals that physically
//       construct a new proposed tree, t', suff stats must be computed
//       on t' using subsuff.  Using getsuff or allsuff is WRONG and will
//       result in undefined behaviour since getsuff/allsuff *assume* the 
//       the root of the tree is brt.t.
void brt::subsuff(tree::tree_p nx, tree::npv& bnv, std::vector<sinfo*>& siv)
{
   #ifdef _OPENMP
   typedef tree::npv::size_type bvsz;
   #endif

   tree::npv path;

   bnv.clear();
   nx->getpathtoroot(path);  //path from n back to root
   nx->getbots(bnv);  //all bots ONLY BELOW node n!!

   #ifdef _OPENMP
      siv.clear(); //need to setup space threads will add into
      siv.resize(bnv.size());
      for(bvsz i=0;i!=bnv.size();i++) siv[i]=newsinfo();
#     pragma omp parallel num_threads(tc)
      local_ompsubsuff(*di,nx,path,bnv,siv); //faster if pass di and bnv by value.
   #else
      diterator diter(di);
      local_subsuff(diter,nx,path,bnv,siv);
   #endif
}

//--------------------------------------------------
//allsuff (2)
void brt::local_ompallsuff(dinfo di, tree::npv bnv,std::vector<sinfo*>& siv)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = di.n;
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);

   std::vector<sinfo*>& tsiv = newsinfovec(); //will be sized in local_allsuff

   diterator diter(&di,beg,end);
   local_allsuff(diter,bnv,tsiv);

#  pragma omp critical
   {
      for(size_t i=0;i<siv.size();i++) *(siv[i]) += *(tsiv[i]);
   }
   
   for(size_t i=0;i<tsiv.size();i++) delete tsiv[i];
   delete &tsiv;
#endif
}
//--------------------------------------------------
//allsuff (3)
void brt::local_allsuff(diterator& diter, tree::npv& bnv,std::vector<sinfo*>& siv)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   siv.clear();
   siv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) { bnmap[bnv[i]]=i; siv[i]=newsinfo(); }

   for(;diter<diter.until();diter++) {
      tbn = t.bn(diter.getxp(),*xi);
      ni = bnmap[tbn];
      //siv[ni].n +=1; 
      add_observation_to_suff(diter, *(siv[ni]));
   }
}
//--------------------------------------------------
//get suff stats for nodes related to change of variable proposal
//this is simply the allsuff for all nodes under the perturb node, not the entire tree.
void brt::getchgvsuff(tree::tree_p pertnode, tree::npv& bnv, size_t oldc, size_t oldv, bool didswap, 
                  std::vector<sinfo*>& sivold, std::vector<sinfo*>& sivnew)
{
   subsuff(pertnode,bnv,sivnew);
   if(didswap) pertnode->swaplr();  //undo the swap so we can calculate the suff stats for the original variable, cutpoint.
   pertnode->setv(oldv);
   pertnode->setc(oldc);
   subsuff(pertnode,bnv,sivold);
}

//--------------------------------------------------
//get suff stats for nodes related to perturb proposal
//this is simply the allsuff for all nodes under the perturb node, not the entire tree.
void brt::getpertsuff(tree::tree_p pertnode, tree::npv& bnv, size_t oldc, 
                  std::vector<sinfo*>& sivold, std::vector<sinfo*>& sivnew)
{
   subsuff(pertnode,bnv,sivnew);
   pertnode->setc(oldc);
   subsuff(pertnode,bnv,sivold);
}
//--------------------------------------------------
//getsuff wrapper used for death.  Calls serial or parallel code depending on how
//the code is compiled.
void brt::getsuff(tree::tree_p l, tree::tree_p r, sinfo& sil, sinfo& sir)
{
   #ifdef _OPENMP
#     pragma omp parallel num_threads(tc)
      local_ompgetsuff(l,r,*di,sil,sir); //faster if pass dinfo by value.
   #else
         diterator diter(di);
         local_getsuff(diter,l,r,sil,sir);
   #endif
}

//--------------------------------------------------
//--------------------------------------------------
//#ifdef _OPENMP
//--------------------------------------------------
//openmp version of getsuff for birth
void brt::local_ompgetsuff(tree::tree_p nx, size_t v, size_t c, dinfo di, sinfo& sil, sinfo& sir)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = di.n;
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);

   sinfo& tsil = *newsinfo();
   sinfo& tsir = *newsinfo();

   diterator diter(&di,beg,end);
   local_getsuff(diter,nx,v,c,tsil,tsir);

#  pragma omp critical
   {
      sil+=tsil; sir+=tsir;
   }

   delete &tsil;
   delete &tsir;
#endif
}
//--------------------------------------------------
//opemmp version of getsuff for death
void brt::local_ompgetsuff(tree::tree_p l, tree::tree_p r, dinfo di, sinfo& sil, sinfo& sir)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = di.n;
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);

//   sinfo tsil, tsir;
   sinfo& tsil = *newsinfo();
   sinfo& tsir = *newsinfo();

   diterator diter(&di,beg,end);
   local_getsuff(diter,l,r,tsil,tsir);

#  pragma omp critical
   {
      sil+=tsil; sir+=tsir;
   }

   delete &tsil;
   delete &tsir;
#endif
}
//#endif

//--------------------------------------------------
//--------------------------------------------------
//set the vector of predicted values
void brt::setf() {
   #ifdef _OPENMP
#     pragma omp parallel num_threads(tc)
      local_ompsetf(*di); //faster if pass dinfo by value.
   #else
         diterator diter(di);
         local_setf(diter);
   #endif
}
void brt::local_ompsetf(dinfo di)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = di.n;
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);

   diterator diter(&di,beg,end);
   local_setf(diter);
#endif
}
void brt::local_setf(diterator& diter)
{
   tree::tree_p bn;

   for(;diter<diter.until();diter++) {
      bn = t.bn(diter.getxp(),*xi);
      yhat[*diter] = bn->gettheta();
   }
}
//--------------------------------------------------
//set the vector of residual values
void brt::setr() {
   #ifdef _OPENMP
#     pragma omp parallel num_threads(tc)
      local_ompsetr(*di); //faster if pass dinfo by value.
   #else
         diterator diter(di);
         local_setr(diter);
   #endif
}
void brt::local_ompsetr(dinfo di)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = di.n;
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);

   diterator diter(&di,beg,end);
   local_setr(diter);
#endif
}
void brt::local_setr(diterator& diter)
{
   tree::tree_p bn;

   for(;diter<diter.until();diter++) {
      bn = t.bn(diter.getxp(),*xi);
      resid[*diter] = 0.0 - bn->gettheta();
//      resid[*diter] = di->y[*diter] - bn->gettheta();
   }
}
//--------------------------------------------------
//predict the response at the (npred x p) input matrix *x
//Note: the result appears in *dipred.y.
void brt::predict(dinfo* dipred) {
   #ifdef _OPENMP
#     pragma omp parallel num_threads(tc)
      local_omppredict(*dipred); //faster if pass dinfo by value.
   #else
         diterator diter(dipred);
         local_predict(diter);
   #endif
}
void brt::local_omppredict(dinfo dipred)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = dipred.n;
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);

   diterator diter(&dipred,beg,end);
   local_predict(diter);
#endif
}
void brt::local_predict(diterator& diter)
{
   tree::tree_p bn;

   for(;diter<diter.until();diter++) {
      bn = t.bn(diter.getxp(),*xi);
      diter.sety(bn->gettheta());
   }
}
//--------------------------------------------------
//save/load tree to/from vector format
//Note: for single tree models the parallelization just directs
//      back to the serial path (ie no parallel execution occurs).
//      For multi-tree models, the parallelization occurs in the
//      definition of that models class.
//void brt::savetree(int* id, int* v, int* c, double* theta)
void brt::savetree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
   #ifdef _OPENMP
#    pragma omp parallel num_threads(tc)
     local_ompsavetree(iter,m,nn,id,v,c,theta);
   #else
     int beg=0;
     int end=(int)m;
     local_savetree(iter,beg,end,nn,id,v,c,theta);
   #endif
}
//void brt::local_ompsavetree(int* id, int* v, int* c, double* theta)
void brt::local_ompsavetree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = (int)m; //1 tree in brt version of save/load tree(s)
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);
   if(end>my_rank)
      local_savetree(iter,beg,end,nn,id,v,c,theta);
#endif
}
void brt::local_savetree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, 
     std::vector<std::vector<int> >& v, std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
   //beg,end are not used in the single-tree models.
   COUT << "in brt::local_savetree, beg, end: " << beg << ", " << end << "\n";  //dummy to get rid of unused warning
   nn[iter]=t.treesize();
   id[iter].resize(nn[iter]);
   v[iter].resize(nn[iter]);
   c[iter].resize(nn[iter]);
   theta[iter].resize(nn[iter]);
   t.treetovec(&id[iter][0],&v[iter][0],&c[iter][0],&theta[iter][0]);
}
void brt::loadtree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
   #ifdef _OPENMP
#    pragma omp parallel num_threads(tc)
     local_omploadtree(iter,m,nn,id,v,c,theta);
   #else
     int beg=0;
     int end=(int)m;
     local_loadtree(iter,beg,end,nn,id,v,c,theta);
   #endif
}
//void brt::local_omploadtree(size_t nn, int* id, int* v, int* c, double* theta)
void brt::local_omploadtree(size_t iter, size_t m, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
#ifdef _OPENMP
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int n = (int)m; //1 tree in brt version of save/load tree(s)
   int beg=0;
   int end=0;
   calcbegend(n,my_rank,thread_count,&beg,&end);
   if(end>my_rank)
      local_loadtree(iter,beg,end,nn,id,v,c,theta);
#endif
}
void brt::local_loadtree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
   //beg,end are not used in the single-tree models.
   COUT << "in brt::local_savetree, beg, end: " << beg << ", " << end << "\n";  //dummy to get rid of unused warning
   t.vectotree(nn[iter],&id[iter][0],&v[iter][0],&c[iter][0],&theta[iter][0]);
}

//--------------------------------------------------
//--------------------------------------------------
//bd: birth/death
void brt::bd(rn& gen)
{
//   COUT << "--------------->>into bd" << endl;
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(t,*xi,mi.pb,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death
      mi.bproposal++;
      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(t,*xi,tp,mi.pb,goodbots,PBx,nx,v,c,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      sinfo& sil = *newsinfo();
      sinfo& sir = *newsinfo();
      sinfo& sit = *newsinfo();

      getsuff(nx,v,c,sil,sir);
      // sit = sil + sir; NO! The + operator cannot be overloaded, so instead we do this:
      sit += sil;
      sit += sir;

      //--------------------------------------------------
      //compute alpha
      bool hardreject=true;
      double lalpha=0.0;
      double lml, lmr, lmt;  // lm is the log marginal left,right,total
      if((sil.n>=mi.minperbot) && (sir.n>=mi.minperbot)) { 
         lml=lm(sil); lmr=lm(sir); lmt=lm(sit);
         hardreject=false;
         lalpha = log(pr) + (lml+lmr-lmt);
         lalpha = std::min(0.0,lalpha);
      }
      //--------------------------------------------------
      //try metrop
      double thetal,thetar; //parameters for new bottom nodes, left and right
      double uu = gen.uniform();
      if( !hardreject && (log(uu) < lalpha) ) {
         thetal = 0.0;//drawnodetheta(sil,gen);
         thetar = 0.0;//drawnodetheta(sir,gen);
         t.birthp(nx,v,c,thetal,thetar);
         mi.baccept++;
      }
      delete &sil;
      delete &sir;
      delete &sit;
   } else {
      mi.dproposal++;
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree::tree_p nx; //nog node to death at
      dprop(t,*xi,tp,mi.pb,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      //sinfo sil,sir,sit;
      sinfo& sil = *newsinfo();
      sinfo& sir = *newsinfo();
      sinfo& sit = *newsinfo();
      getsuff(nx->getl(),nx->getr(),sil,sir);
      // sit = sil + sir; NO! The + operator cannot be overloaded, so instead we do this:
      sit += sil;
      sit += sir;

      //--------------------------------------------------
      //compute alpha
      double lml, lmr, lmt;  // lm is the log marginal left,right,total
      lml=lm(sil); lmr=lm(sir); lmt=lm(sit);
      double lalpha = log(pr) + (lmt - lml - lmr);
      lalpha = std::min(0.0,lalpha);

      //--------------------------------------------------
      //try metrop
      double theta;
      if(log(gen.uniform()) < lalpha) {
         theta = 0.0;//drawnodetheta(sit,gen);
         t.deathp(nx,theta);
         mi.daccept++;
      }
      delete &sil;
      delete &sir;
      delete &sit;
   }
}
//--------------------------------------------------
//peturb proposal for internal node cut points.
void brt::pertcv(rn& gen)
{
//   COUT << "--------------->>into pertcv" << endl;
   tree::tree_p pertnode;
   if(t.treesize()==1) // nothing to perturb if the tree is a single terminal node
      return;

   // Get interior nodes and propose new split value
   tree::npv intnodes;
   t.getintnodes(intnodes);
   for(size_t pertdx=0;pertdx<intnodes.size();pertdx++)
   if(di->p > 1 && gen.uniform()<mi.pchgv) {
      mi.chgvproposal++;
      pertnode = intnodes[pertdx];

      //get L,U for the old variable and save it as well as oldc
      int Lo,Uo;
      getLU(pertnode,*xi,&Lo,&Uo);
      size_t oldc = pertnode->getc();

      //update correlation matrix
      std::vector<std::vector<double> > chgv;
      chgv= *mi.corv; //initialize it
      bool didswap=false;
      size_t oldv=pertnode->getv();
      updatecormat(pertnode,*xi,chgv);
      normchgvrow(oldv,chgv);

      //choose new variable randomly
      size_t newv=getchgv(oldv,chgv,gen);
      pertnode->setv(newv);
      if(chgv[oldv][newv]<0.0) {
         pertnode->swaplr();
         didswap=true;
      }

      //get L,U for the new variable and save it and set newc
      int Ln,Un;
      getLU(pertnode,*xi,&Ln,&Un);
      size_t newc = Ln + (size_t)(floor(gen.uniform()*(Un-Ln+1.0)));
      pertnode->setc(newc);

      //now we also need to update the row of chgv for newv->oldv to calc MH correctly
      updatecormat(pertnode,*xi,chgv);
      normchgvrow(newv,chgv);
      
      //sanity check:
      if(chgv[newv][oldv]==0.0)
         COUT << "Proposal newv cannot return to oldv!  This is not possible!" << endl;
      double alpha0=chgv[newv][oldv]/chgv[oldv][newv];  //proposal ratio for newv->oldv and oldv->newv

      //get sufficient statistics and calculate lm
      std::vector<sinfo*>& sivold = newsinfovec();
      std::vector<sinfo*>& sivnew = newsinfovec();
      tree::npv bnv;
      getchgvsuff(pertnode,bnv,oldc,oldv,didswap,sivold,sivnew);

      typedef tree::npv::size_type bvsz;
      double lmold,lmnew;
      bool hardreject=false;
      lmold=0.0;
      for(bvsz j=0;j!=sivold.size();j++) {
         if(sivold[j]->n < mi.minperbot)
            COUT << "Error: old tree has some bottom nodes with <minperbot observations!" << endl;
         lmold += lm(*(sivold[j]));
      }

      lmnew=0.0;
      for(bvsz j=0;j!=sivnew.size();j++) {
         if(sivnew[j]->n < mi.minperbot)
            hardreject=true;
         lmnew += lm(*(sivnew[j]));
      }
      double alpha1 = ((double)(Uo-Lo+1.0))/((double)(Un-Ln+1.0)); //from prior for cutpoints
      double alpha2=alpha0*alpha1*exp(lmnew-lmold);
      double alpha = std::min(1.0,alpha2);
      if(hardreject) alpha=0.0;  //change of variable led to an bottom node with <minperbot observations in it, we reject this.

      if(gen.uniform()<alpha) {
         mi.chgvaccept++;
         if(didswap) pertnode->swaplr();  //because the call to getchgvsuff unswaped if they were swapped
         pertnode->setv(newv); //because the call to getchgvsuff changes it back to oldv to calc the old lil
         pertnode->setc(newc); //because the call to getchgvsuff changes it back to oldc to calc the old lil
      }
      //else nothing, pertnode->c and pertnode->v is already reset to the old values and if a swap was done in the 
      //proposal it was already undone by getchgvsuff.
      for(bvsz j=0;j<sivold.size();j++) delete sivold[j];
      for(bvsz j=0;j<sivnew.size();j++) delete sivnew[j];
      delete &sivold;
      delete &sivnew;
   }
   else {
      mi.pertproposal++;
      pertnode = intnodes[pertdx];

      // Get allowable range for perturbing cv at pertnode
      int L,U;
      bool hardreject=false;
      getLU(pertnode,*xi,&L,&U);
      size_t oldc = pertnode->getc();
      int ai,bi,oldai,oldbi;
      ai=(int)(floor(oldc-mi.pertalpha*(U-L+1)/2.0));
      bi=(int)(floor(oldc+mi.pertalpha*(U-L+1)/2.0));
      ai=std::max(ai,L);
      bi=std::min(bi,U);
      size_t propc = ai + (size_t)(floor(gen.uniform()*(bi-ai+1.0)));
      pertnode->setc(propc);
      oldai=(int)(floor(propc-mi.pertalpha*(U-L+1)/2.0));
      oldbi=(int)(floor(propc+mi.pertalpha*(U-L+1)/2.0));
      oldai=std::max(oldai,L);
      oldbi=std::min(oldbi,U);

      //std::vector<sinfo> sivold, sivnew;
      std::vector<sinfo*>& sivold = newsinfovec();
      std::vector<sinfo*>& sivnew = newsinfovec();

      tree::npv bnv;
      getpertsuff(pertnode,bnv,oldc,sivold,sivnew);

      typedef tree::npv::size_type bvsz;
      double lmold,lmnew;
      lmold=0.0;
      for(bvsz j=0;j!=sivold.size();j++) {
         if(sivold[j]->n < mi.minperbot)
            COUT << "Error: old tree has some bottom nodes with <minperbot observations!" << endl;
         lmold += lm(*(sivold[j]));
      }

      lmnew=0.0;
      for(bvsz j=0;j!=sivnew.size();j++) {
         if(sivnew[j]->n < mi.minperbot)
            hardreject=true;
         lmnew += lm(*(sivnew[j]));
      }
      double alpha1 = ((double)(bi-ai+1.0))/((double)(oldbi-oldai+1.0)); //anything from the prior?
      double alpha2=alpha1*exp(lmnew-lmold);
      double alpha = std::min(1.0,alpha2);
      if(hardreject) alpha=0.0;  //perturb led to an bottom node with <minperbot observations in it, we reject this.

      if(gen.uniform()<alpha) {
         mi.pertaccept++;
         pertnode->setc(propc); //because the call to getpertsuff changes it back to oldc to calc the old lil.
      }
      //else nothing, pertnode->c is already reset to the old value.
      for(bvsz j=0;j<sivold.size();j++) delete sivold[j];
      for(bvsz j=0;j<sivnew.size();j++) delete sivnew[j];
      delete &sivold;
      delete &sivnew;
   }
}

//--------------------------------------------------
//do a rotation proposal at a randomly selected internal node.
bool brt::rot(tree::tree_p tnew, tree& x, rn& gen)
{
//   COUT << "--------------->>into rot" << endl;
   tree::tree_p rotp,temp;
   tree::tree_cp xp;
   tree::npv subtold, subtnew, nbold, nbnew;
   double Qold_to_new, Qnew_to_old;
   unsigned int rdx=0;
   bool twowaystoinvert=false;
   double prinew=1.0,priold=1.0;
   size_t rotid;
   bool hardreject=false;
   std::vector<size_t> goodvars; //variables an internal node can split on

   mi.rotproposal++;

   // Get rot nodes
   tree::npv rnodes;
   tnew->getrotnodes(rnodes);
   if(rnodes.size()==0)  return false;  //no rot nodes so that's a reject.

   rdx = (unsigned int)floor(gen.uniform()*rnodes.size()); //which rotatable node will we rotate at?
   rotp = rnodes[rdx];
   rotid=rotp->nid();
   xp=x.getptr(rotid);

/* Can check the funcitonality of getpathtoroot:  
   tree::npv path;
   rotp->getpathtoroot(path);
   COUT << "rot id=" << rotid << endl;
   tnew->pr();
   for(size_t i=0;i<path.size();i++)
      COUT << "i=" << i << ", node id=" << path[i]->nid() << endl;
*/

   int nwaysm1=0,nwaysm2=0,nwayss1=0,nwayss2=0;
   double pm1=1.0,pm2=1.0,ps1=1.0,ps2=1.0;
   if(rotp->isleft()) {
      if(rotp->v==rotp->p->v) //special case, faster to handle it direclty
      {
         rotright(rotp);
         rotp=tnew->getptr(rotid);
         delete rotp->r;
         temp=rotp->l;
         rotp->p->l=temp;
         temp->p=rotp->p;
         rotp->r=0;
         rotp->l=0;
         rotp->p=0;
         delete rotp;
         rotp=tnew->getptr(rotid);
         //pm1=pm2=ps1=ps2=1.0 in this case
      }
      else
      {
         rotright(rotp);
         rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
         reduceleft(rotp,rotp->p->v,rotp->p->c);
         rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
         reduceright(rotp->p->r,rotp->p->v,rotp->p->c);
         rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
         splitleft(rotp->r,rotp->p->v,rotp->p->c);
         splitright(rotp->p->r->r,rotp->p->v,rotp->p->c);

         mergecount(rotp->r,rotp->p->r->r,rotp->p->v,rotp->p->c,&nwayss1);
         ps1=1.0/nwayss1;

         mergecount(rotp->l,rotp->p->r->l,rotp->p->v,rotp->p->c,&nwayss2);
         ps2=1.0/nwayss2;

         tree::tree_p tmerged=new tree;
         tmerged->p=rotp->p;

         mergecount(rotp->p->r->l,rotp->p->r->r,rotp->p->r->v,rotp->p->r->c,&nwaysm1);
         pm1=1.0/nwaysm1;
         merge(rotp->p->r->l,rotp->p->r->r,tmerged,rotp->p->r->v,rotp->p->r->c,gen);
         rotp->p->r->p=0;
         delete rotp->p->r;
         rotp->p->r=tmerged;

         tmerged=new tree;
         tmerged->p=rotp->p;

         mergecount(rotp->l,rotp->r,rotp->v,rotp->c,&nwaysm2);
         pm2=1.0/nwaysm2;
         size_t v,c;
         v=rotp->v;
         c=rotp->c;
         merge(rotp->l,rotp->r,tmerged,rotp->v,rotp->c,gen);
         rotp->p->l=tmerged;
         rotp->p=0;
         delete rotp;
         rotp=tnew->getptr(rotid);

      //end of merge code if rotp isleft.
      //there are some "extra" isleaf's here because we don't explicitly reset v,c if node becomes leaf so we need to check.
         if( !isleaf(rotp) && !isleaf(rotp->p->r) && (rotp->v!=v && rotp->c!=c) && (rotp->p->r->v != v && rotp->p->r->c != c))
            hardreject=true;
         if( isleaf(rotp) && isleaf(rotp->p->r))
            hardreject=true;
         if(rotp->p->r->v==rotp->v && rotp->p->r->c==rotp->c && !isleaf(rotp->p->r) && !isleaf(rotp))
            twowaystoinvert=true;

      }
   }
   else { //isright
      if(rotp->v==rotp->p->v) //special case, faster to handle it directly
      {
         rotleft(rotp);
         rotp=tnew->getptr(rotid);
         delete rotp->l;
         temp=rotp->r;
         rotp->p->r=temp;
         temp->p=rotp->p;
         rotp->r=0;
         rotp->l=0;
         rotp->p=0;
         delete rotp;
         rotp=tnew->getptr(rotid);
         //pm1=pm2=ps1=ps2=1.0 in this case
      }
      else
      {
         rotleft(rotp);
         rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
         reduceleft(rotp->p->l,rotp->p->v,rotp->p->c);
         rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
         reduceright(rotp,rotp->p->v,rotp->p->c);
         rotp=tnew->getptr(rotid); //just in case the above changed the pointer.
         splitleft(rotp->p->l->l,rotp->p->v,rotp->p->c);
         splitright(rotp->l,rotp->p->v,rotp->p->c);

         mergecount(rotp->p->l->l,rotp->l,rotp->p->v,rotp->p->c,&nwayss1);
         ps1=1.0/nwayss1;

         mergecount(rotp->p->l->r,rotp->r,rotp->p->v,rotp->p->c,&nwayss2);
         ps2=1.0/nwayss2;

         tree::tree_p tmerged=new tree;
         tmerged->p=rotp->p;

         mergecount(rotp->p->l->l,rotp->p->l->r,rotp->p->l->v,rotp->p->l->c,&nwaysm1);
         pm1=1.0/nwaysm1;
         merge(rotp->p->l->l,rotp->p->l->r,tmerged,rotp->p->l->v,rotp->p->l->c,gen);
         rotp->p->l->p=0;
         delete rotp->p->l;
         rotp->p->l=tmerged;

         tmerged=new tree;
         tmerged->p=rotp->p;

         mergecount(rotp->l,rotp->r,rotp->v,rotp->c,&nwaysm2);
         pm2=1.0/nwaysm2;
         size_t v,c;
         v=rotp->v;
         c=rotp->c;
         merge(rotp->l,rotp->r,tmerged,rotp->v,rotp->c,gen);
         rotp->p->r=tmerged;
         rotp->p=0;
         delete rotp;
         rotp=tnew->getptr(rotid);

      //end of merge code if rotp isright
      //there are some "extra" isleaf's here because we don't explicitly reset v,c if node becomes leaf so we need to check.
         if( !isleaf(rotp) && !isleaf(rotp->p->l) && (rotp->v!=v && rotp->c!=c) && (rotp->p->l->v != v && rotp->p->l->c != c))
            hardreject=true;
         if( isleaf(rotp) && isleaf(rotp->p->l))
            hardreject=true;
         if(rotp->p->l->v==rotp->v && rotp->p->l->c==rotp->c && !isleaf(rotp->p->l) && !isleaf(rotp))
            twowaystoinvert=true;
      }
   }

   // Calculate prior probabilities, we just need to use the subtree where the rotation occured of tnew and x.
   subtold.clear();
   subtnew.clear();
   xp->p->getnodes(subtold);
   rotp->p->getnodes(subtnew);

   for(size_t i=0;i<subtold.size();i++) {
      if(subtold[i]->l) { //interior node
         priold*=tp.alpha/pow(1.0 + subtold[i]->depth(),tp.beta);
         goodvars.clear();
         getinternalvars(subtold[i],*xi,goodvars);
         priold*=1.0/((double)goodvars.size()); //prob split on v 
         priold*=1.0/((double)getnumcuts(subtold[i],*xi,subtold[i]->v)); //prob split on v at c is 1/numcutpoints
      }
      else //terminal node
         priold*=(1.0-tp.alpha/pow(1.0 + subtold[i]->depth(),tp.beta)); 
   }
   for(size_t i=0;i<subtnew.size();i++) {
      if(subtnew[i]->l) { //interior node
         prinew*=tp.alpha/pow(1.0 + subtnew[i]->depth(),tp.beta);
         goodvars.clear();
         getinternalvars(subtnew[i],*xi,goodvars);
         prinew*=1.0/((double)goodvars.size()); //prob split on v
         prinew*=1.0/((double)getnumcuts(subtnew[i],*xi,subtnew[i]->v)); //prob split on v at c is 1/numcutpoints
         if(getnumcuts(subtnew[i],*xi,subtnew[i]->v)<1)
         {
            x.pr(true);
            tnew->pr(true);
         }
      }
      else //terminal node
         prinew*=(1.0-tp.alpha/pow(1.0 + subtnew[i]->depth(),tp.beta)); 
   }

   Qold_to_new=1.0/((double)rnodes.size()); //proposal probability of rotating from x to tnew
   
   rnodes.clear();
   tnew->getrotnodes(rnodes);  //this is very inefficient, could make it much nicer later on.
//   if(rnodes.size()==0) hardreject=true; //if we're back down to just a root node we can't transition back, so this is a hard reject.

   if(!twowaystoinvert)
      Qnew_to_old=1.0/((double)rnodes.size()); //proposal probability of rotating from tnew back to x
   else
      Qnew_to_old=2.0/((double)rnodes.size());

   // Calculate log integrated likelihoods for the subtree where the rotation occured of tnew and x.
   double lmold=0.0,lmnew=0.0;
//   std::vector<sinfo> sold,snew;
   std::vector<sinfo*>& sold = newsinfovec();
   std::vector<sinfo*>& snew = newsinfovec();
   nbold.clear();
   nbnew.clear();
   sold.clear();
   snew.clear();
   x.getbots(nbold);
   tnew->getbots(nbnew);

   //get sufficient statistics for subtree involved in rotation
   //which is just everything below rotp->p.
   //Use subsuff here, which will get the suff stats for both the
   //orignal tree and the proposed tree without needed to explicitly
   //know the root node of either tree as this is recovered when
   //finding the path to rotp->p within the subsuff method.
   rotp=x.getptr(rotid);
   subsuff(rotp->p,nbold,sold);
   rotp=tnew->getptr(rotid);
   subsuff(rotp->p,nbnew,snew);

   for(size_t i=0;i<nbold.size();i++)
         lmold += lm(*(sold[i]));

   for(size_t i=0;i<nbnew.size();i++) {
      if( (snew[i]->n) >= mi.minperbot )
         lmnew += lm(*(snew[i]));
      else 
         hardreject=true;
   }

   for(size_t i=0;i<sold.size();i++) delete sold[i];
   for(size_t i=0;i<snew.size();i++) delete snew[i];
   delete &sold;
   delete &snew;

   double alpha1;
   alpha1=prinew*Qnew_to_old/priold/Qold_to_new/pm1/pm2*ps1*ps2;
   double alpha2 = alpha1*exp(lmnew-lmold);
   double alpha = std::min(1.0,alpha2);

   if(hardreject)
      alpha=0.0;

   if(gen.uniform()<alpha) {
      mi.rotaccept++;
      x = *tnew;
      return true;
   }
   else {
      return false;
   }

   return false;  // we never actually get here.
}

#endif

