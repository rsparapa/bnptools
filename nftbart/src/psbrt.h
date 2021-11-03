/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         and Hugh A. Chipman
 *  
 * This file is part of nftbart.
 * psbrt.h
 *
 * nftbart is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * nftbart is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author contact information
 * Matthew T. Pratola: mpratola@gmail.com
 * Robert E. McCulloch: robert.e.mculloch@gmail.com
 * Hugh A. Chipman: hughchipman@gmail.com
 *
 */

#ifndef GUARD_psbrt_h
#define GUARD_psbrt_h

/*
#include "tree.h"
#include "treefuns.h"
#include "dinfo.h"
#include "sbrt.h"
//#include "brtfuns.h"
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
*/

class psbrt : public sbrt 
{
public:
   //--------------------
   //classes
   // tprior and mcmcinfo are same as in sbrt

   //--------------------
   //constructors/destructors
   psbrt(): sbrt(),m(10),sb(m),notjsigmavs(m),divec(m) {}
   psbrt(size_t im): sbrt(),m(im),sb(m),notjsigmavs(m),divec(m) {}
   psbrt(size_t im, double itheta): sbrt(pow(itheta,1/im)),m(im),sb(m),notjsigmavs(m),divec(m) {}
   virtual ~psbrt() {
      if(!notjsigmavs.empty()) {
         for(size_t j=0;j<m;j++) notjsigmavs[j].clear();
         notjsigmavs.clear();
         for(size_t j=0;j<m;j++) delete divec[j];
      }
   }

   //--------------------
   //methods
   void draw(rn& gen);
   void adapt();
   void setci(double nu, double lambda) { ci.nu=nu; ci.lambda=lambda; for(size_t j=0;j<m;j++) sb[j].setci(nu,lambda); }
   void settc(int tc) { this->tc = tc; for(size_t j=0;j<m;j++) sb[j].settc(tc); }
   void setxi(xinfo *xi) { this->xi=xi; for(size_t j=0;j<m;j++) sb[j].setxi(xi); }
   void setdata(dinfo *di);
   void settp(double alpha, double beta) { tp.alpha=alpha;tp.beta=beta; for(size_t j=0;j<m;j++) sb[j].settp(alpha,beta); }
   tree::tree_p gettree(size_t i) { return &sb[i].t; } 
   void setmi(double pbd, double pb, size_t minperbot, bool dopert, double pertalpha, double pchgv, std::vector<std::vector<double> >* chgv)
             { mi.pbd=pbd; mi.pb=pb; mi.minperbot=minperbot; mi.dopert=dopert;
               mi.pertalpha=pertalpha; mi.pchgv=pchgv; mi.corv=chgv; 
               for(size_t j=0;j<m;j++) sb[j].setmi(pbd,pb,minperbot,dopert,pertalpha,pchgv,chgv); }
   void setstats(bool dostats) { mi.dostats=dostats; for(size_t j=0;j<m;j++) sb[j].setstats(dostats); if(dostats) mi.varcount=new int[xi->size()];
   this->resetstats(); }
   //void setstats(bool dostats) { mi.dostats=dostats; for(size_t j=0;j<m;j++) sb[j].setstats(dostats); if(dostats) mi.varcount=new unsigned int[xi->size()]; }
   void pr();
   // drawnodetheta, lm, add_observation_to_suff and newsinfo/newsinfovec unused here.

   //--------------------
   //data
   //--------------------------------------------------
   //stuff that maybe should be protected
protected:
   //--------------------
   //model information
   size_t m;  //number of trees in product representation
   std::vector<sbrt> sb;  // the vector of individual sigma trees for product representation
   //--------------------
   //data
   std::vector<std::vector<double> > notjsigmavs;
   std::vector<dinfo*> divec;
   //--------------------
   //mcmc info
   //--------------------
   //methods
   virtual void local_setf(diterator& diter);  //set the vector of predicted values
   virtual void local_setr(diterator& diter);  //set the vector of residuals
   virtual void local_predict(diterator& diter); // predict y at the (npred x p) settings *di.x
   virtual void local_savetree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);
   virtual void local_loadtree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta);

};


//--------------------------------------------------
//a single iteration of the MCMC for brt model
void psbrt::draw(rn& gen)
{
  for(size_t j=0;j<m;j++) {
    //update this row of notjsigmavs
//    for(size_t i=0;i<di->n;i++) {
//      notjsigmavs[j][i]=r(i)*sb[j].f(i);
//    }
    *divec[j]= *getr();
    *divec[j]*= *sb[j].getf();

    // do the draw for jth component
    sb[j].draw(gen);

    // Update the in-sample predicted vector
    setf();

    // Update the in-sample residual vector
    setr();
  }
  // overall statistics from the subtrees.  Need to divide by m*N to get
  // useful numbers after the MCMC is done.
  if(mi.dostats) {
    resetstats();
    for(size_t j=0;j<m;j++)
      sb[j].addstats(mi.varcount,&mi.tavgd,&mi.tmaxd,&mi.tmind);
  }
}
//--------------------------------------------------
//adapt the proposal widths for perturb proposals,
//bd or rot proposals and b or d proposals.
void psbrt::adapt()
{
  for(size_t j=0;j<m;j++) {
    //COUT << "\nAdapt sbrt[" << j << "]:";
    sb[j].adapt();
  }
}
//--------------------------------------------------
//setdata for psbrt
void psbrt::setdata(dinfo *di) {
  this->di=di;

  // initialize notjsigmavs.
  for(size_t j=0;j<m;j++)
      notjsigmavs[j].resize(this->di->n,1.0);
  for(size_t j=0;j<m;j++)
    for(size_t i=0;i<di->n;i++)
      notjsigmavs[j][i]=pow(std::abs(this->di->y[i]/.8),1.0/m);  //E(|Y|) = .8sigma for normal

  for(size_t j=0;j<m;j++)
    divec[j]=new dinfo(this->di->p,this->di->n,this->di->x,&notjsigmavs[j][0],this->di->tc);

  // each sb[j]'s data is the appropriate row in notjsigmavs
  for(size_t j=0;j<m;j++)
    sb[j].setdata(divec[j]);

  resid.resize(di->n);
  yhat.resize(di->n);
  setf();
  setr();
}
//--------------------------------------------------
//set vector of predicted values for psbrt model
void psbrt::local_setf(diterator& diter)
{
   for(;diter<diter.until();diter++) {
      yhat[*diter]=1.0;
      for(size_t j=0;j<m;j++)
        yhat[*diter]*=sb[j].f(*diter);
   }
}
//--------------------------------------------------
//set vector of residuals for psbrt model
void psbrt::local_setr(diterator& diter)
{
   for(;diter<diter.until();diter++) {
      resid[*diter]=di->y[*diter]/f(*diter);
   }
}
//--------------------------------------------------
//predict the response at the (npred x p) input matrix *x
//Note: the result appears in *dipred.y.
void psbrt::local_predict(diterator& diter)
{
  tree::tree_p bn;
  double temp;

  for(;diter<diter.until();diter++) {
    temp=1.0;
    for(size_t j=0;j<m;j++) {
      bn = sb[j].t.bn(diter.getxp(),*xi);
      temp*=bn->gettheta();
    }
    diter.sety(temp);
  }
}
void psbrt::local_savetree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, 
     std::vector<std::vector<int> >& v, std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
  size_t indx=iter*m;
  for(size_t i=(indx+(size_t)beg);i<(indx+(size_t)end);i++) {
    nn[i]=sb[i-indx].t.treesize();
    id[i].resize(nn[i]);
    v[i].resize(nn[i]);
    c[i].resize(nn[i]);
    theta[i].resize(nn[i]);
    sb[i-indx].t.treetovec(&id[i][0],&v[i][0],&c[i][0],&theta[i][0]);
  }
}
void psbrt::local_loadtree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
  size_t indx=iter*m;
  for(size_t i=(indx+(size_t)beg);i<(indx+(size_t)end);i++)
    sb[i-indx].t.vectotree(nn[i],&id[i][0],&v[i][0],&c[i][0],&theta[i][0]);
}
//--------------------------------------------------
//pr for brt
void psbrt::pr()
{
   COUT << "***** psbrt object:\n";
   COUT << "Number of trees in product representation:" << endl;
   COUT << "        m:   m=" << m << endl;
   COUT << "Conditioning info on each individual tree:" << endl;
   COUT << "      dof:  nu=" << ci.nu << endl;
   COUT << "    scale:  lambda=" << ci.lambda << endl;
   brt::pr();
   COUT << "**************Trees in product representation*************:" << endl;
   for(size_t j=0;j<m;j++) sb[j].t.pr();
}

#endif
