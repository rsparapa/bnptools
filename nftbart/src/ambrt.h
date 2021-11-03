/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         and Hugh A. Chipman
 *  
 * This file is part of nftbart.
 * ambrt.h
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

#ifndef GUARD_ambrt_h
#define GUARD_ambrt_h

/*
#include "tree.h"
#include "treefuns.h"
#include "dinfo.h"
#include "mbrt.h"
//#include "brtfuns.h"
#include <iostream>
#include <map>
#include <vector>
*/

class ambrt : public mbrt 
{
public:
   //--------------------
   //classes
   // tprior and mcmcinfo are same as in brt
   // cinfo same as in mbrt
   //--------------------
   //constructors/destructors
   ambrt(): mbrt(),m(200),mb(m),notjmus(m),divec(m) {}
   ambrt(size_t im): mbrt(),m(im),mb(m),notjmus(m),divec(m) {}
   virtual ~ambrt() {
      if(!notjmus.empty()) {
         for(size_t j=0;j<m;j++) notjmus[j].clear();
         notjmus.clear();
         for(size_t j=0;j<m;j++) delete divec[j];
      }
   }

   //--------------------
   //methods
   void draw(rn& gen);
   void adapt();
   void setci(double tau, double* sigma) { ci.tau=tau; ci.sigma=sigma; for(size_t j=0;j<m;j++) mb[j].setci(tau,sigma); }
   void settc(int tc) { this->tc = tc; for(size_t j=0;j<m;j++) mb[j].settc(tc); }
   void setxi(xinfo *xi) { this->xi=xi; for(size_t j=0;j<m;j++) mb[j].setxi(xi); }
   void setdata(dinfo *di);
   void settp(double alpha, double beta) { tp.alpha=alpha;tp.beta=beta; for(size_t j=0;j<m;j++) mb[j].settp(alpha,beta); }
   tree::tree_p gettree(size_t i) { return &mb[i].t; }
   void setmi(double pbd, double pb, size_t minperbot, bool dopert, double pertalpha, double pchgv, std::vector<std::vector<double> >* chgv)
             { mi.pbd=pbd; mi.pb=pb; mi.minperbot=minperbot; mi.dopert=dopert;
               mi.pertalpha=pertalpha; mi.pchgv=pchgv; mi.corv=chgv; 
               for(size_t j=0;j<m;j++) mb[j].setmi(pbd,pb,minperbot,dopert,pertalpha,pchgv,chgv); }
   void setstats(bool dostats) { mi.dostats=dostats; for(size_t j=0;j<m;j++) mb[j].setstats(dostats); if(dostats) mi.varcount=new int[xi->size()]; 
     this->resetstats(); }
   //void setstats(bool dostats) { mi.dostats=dostats; for(size_t j=0;j<m;j++) mb[j].setstats(dostats); if(dostats) mi.varcount=new unsigned int[xi->size()]; }
   void pr();
   // drawnodetheta, lm, add_observation_to_suff and newsinfo/newsinfovec unused here.

   //--------------------
   //data
   //--------------------------------------------------
   //stuff that maybe should be protected
   // std::vector<mbrt> mb;  // the vector of individual mu trees for sum representation
protected:
   //--------------------
   //model information
   size_t m;  //number of trees in sum representation
   std::vector<mbrt> mb;  // the vector of individual mu trees for sum representation
   //--------------------
   //data
   std::vector<std::vector<double> > notjmus;
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
void ambrt::draw(rn& gen)
{
  for(size_t j=0;j<m;j++) {
    //update this row of notjmus
    // for(size_t i=0;i<di->n;i++) {
    //   notjmus[j][i]=di->y[i]-f(i)+mb[j].f(i);  //res_j=y-sum_{k!=j}tree_j
    // }
   *divec[j]= *di;
   *divec[j]-= *getf();
   *divec[j]+= *mb[j].getf();

    // do the draw for jth component
    mb[j].draw(gen);

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
      mb[j].addstats(mi.varcount,&mi.tavgd,&mi.tmaxd,&mi.tmind);
  }
}
//--------------------------------------------------
//adapt the proposal widths for perturb proposals,
//bd or rot proposals and b or d proposals.
void ambrt::adapt()
{
  for(size_t j=0;j<m;j++) {
    //COUT << "\nAdapt ambrt[" << j << "]:";
    mb[j].adapt();
  }
}
//--------------------------------------------------
//setdata for ambrt
void ambrt::setdata(dinfo *di) {
  this->di=di;

  // initialize notjsigmavs.
  for(size_t j=0;j<m;j++)
      notjmus[j].resize(this->di->n,0.0);
  for(size_t j=0;j<m;j++)
    for(size_t i=0;i<di->n;i++)
      notjmus[j][i]=this->di->y[i]/((double)m);

  for(size_t j=0;j<m;j++)
    divec[j]=new dinfo(this->di->p,this->di->n,this->di->x,&notjmus[j][0],this->di->tc);

  // each mb[j]'s data is the appropriate row in notjmus
  for(size_t j=0;j<m;j++)
    mb[j].setdata(divec[j]);

  resid.resize(di->n);
  yhat.resize(di->n);
  setf();
  setr();
}
//--------------------------------------------------
//set vector of predicted values for psbrt model
void ambrt::local_setf(diterator& diter)
{
   for(;diter<diter.until();diter++) {
      yhat[*diter]=0.0;
      for(size_t j=0;j<m;j++)
        yhat[*diter]+=mb[j].f(*diter);
   }
}
//--------------------------------------------------
//set vector of residuals for psbrt model
void ambrt::local_setr(diterator& diter)
{
   for(;diter<diter.until();diter++) {
      resid[*diter]=di->y[*diter]-f(*diter);
   }
}
//--------------------------------------------------
//predict the response at the (npred x p) input matrix *x
//Note: the result appears in *dipred.y.
void ambrt::local_predict(diterator& diter)
{
  tree::tree_p bn;
  double temp;

  for(;diter<diter.until();diter++) {
    temp=0.0;
    for(size_t j=0;j<m;j++) {
      bn = mb[j].t.bn(diter.getxp(),*xi);
      temp+=bn->gettheta();
    }
    diter.sety(temp);
  }
}
void ambrt::local_savetree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, 
     std::vector<std::vector<int> >& v, std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
  size_t indx=iter*m;
  for(size_t i=(indx+(size_t)beg);i<(indx+(size_t)end);i++) {
    nn[i]=mb[i-indx].t.treesize();
    id[i].resize(nn[i]);
    v[i].resize(nn[i]);
    c[i].resize(nn[i]);
    theta[i].resize(nn[i]);
    mb[i-indx].t.treetovec(&id[i][0],&v[i][0],&c[i][0],&theta[i][0]);
  }
}
void ambrt::local_loadtree(size_t iter, int beg, int end, std::vector<int>& nn, std::vector<std::vector<int> >& id, std::vector<std::vector<int> >& v,
                  std::vector<std::vector<int> >& c, std::vector<std::vector<double> >& theta)
{
  size_t indx=iter*m;
  for(size_t i=(indx+(size_t)beg);i<(indx+(size_t)end);i++)
    mb[i-indx].t.vectotree(nn[i],&id[i][0],&v[i][0],&c[i][0],&theta[i][0]);
}

//--------------------------------------------------
//pr for brt
void ambrt::pr()
{
   COUT << "***** ambrt object:\n";
   COUT << "Number of trees in product representation:" << endl;
   COUT << "        m:   m=" << m << endl;
   COUT << "Conditioning info on each individual tree:" << endl;
   COUT << "   mean:   tau=" << ci.tau << endl;
   if(!ci.sigma)
     COUT << "         sigma=[]" << endl;
   else
     COUT << "         sigma=[" << ci.sigma[0] << ",...," << ci.sigma[di->n-1] << "]" << endl;
   brt::pr();
   COUT << "**************Trees in sum representation*************:" << endl;
   for(size_t j=0;j<m;j++) mb[j].t.pr();
}

#endif
