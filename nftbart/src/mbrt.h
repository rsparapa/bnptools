/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         Hugh A. Chipman and Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * mbrt.h
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
 * Rodney A. Sparapani: rsparapa@mcw.edu
 *
 */

#ifndef GUARD_mbrt_h
#define GUARD_mbrt_h

/*
#include "tree.h"
#include "treefuns.h"
#include "dinfo.h"
#include "brt.h"
//#include "brtfuns.h"
#include <iostream>
#include <map>
#include <vector>
*/

class msinfo : public sinfo { //sufficient statistics (will depend on end node model)
public:
   msinfo():sinfo(),sumw(0.0),sumwy(0.0) {}
   msinfo(const msinfo& is):sinfo(is),sumw(is.sumw),sumwy(is.sumwy) {}
   virtual ~msinfo() {}  //need this so memory is properly freed in derived classes.
   double sumw;
   double sumwy;
   // compound addition operator needed when adding suff stats
   virtual sinfo& operator+=(const sinfo& rhs) {
      sinfo::operator+=(rhs);
      const msinfo& mrhs=static_cast<const msinfo&>(rhs);
      sumw+=mrhs.sumw;
      sumwy+=mrhs.sumwy;
      return *this;
   }
   // assignment operator for suff stats
   virtual sinfo& operator=(const sinfo& rhs)
   {
      if(&rhs != this) {
         sinfo::operator=(rhs);
         const msinfo& mrhs=static_cast<const msinfo&>(rhs);
         this->sumw = mrhs.sumw;
         this->sumwy = mrhs.sumwy;
      }
      return *this;
   }
   // addition opertor is defined in terms of compound addition
   const msinfo operator+(const msinfo& other) const {
      msinfo result = *this; //copy of myself.
      result += other;
      return result;
   }
};

class mbrt : public brt 
{
public:
   //--------------------
   //classes
   // tprior and mcmcinfo are same as in brt
   class cinfo { //parameters for end node model prior
   public:
      cinfo():tau(1.0),sigma(0) {}
      double tau;
      double* sigma;
   };
   //--------------------
   //constructors/destructors
   mbrt():brt() {}
   //--------------------
   //methods
   void draw(rn& gen);
   void setci(double tau, double* sigma) { ci.tau=tau; ci.sigma=sigma; }
   virtual double drawnodetheta(sinfo& si, rn& gen);
   virtual double lm(sinfo& si);
   virtual void add_observation_to_suff(diterator& diter, sinfo& si);
   virtual sinfo* newsinfo() { return new msinfo; }
   virtual std::vector<sinfo*>& newsinfovec() { std::vector<sinfo*>* si= new std::vector<sinfo*>; return *si; }
   virtual std::vector<sinfo*>& newsinfovec(size_t dim) { std::vector<sinfo*>* si = new std::vector<sinfo*>; si->resize(dim); for(size_t i=0;i<dim;i++) si->push_back(new msinfo); return *si; }
   void pr();

   //--------------------
   //data
   //--------------------------------------------------
   //stuff that maybe should be protected
protected:
   //--------------------
   //model information
   cinfo ci; //conditioning info (e.g. other parameters and prior and end node models)
   //--------------------
   //data
   //--------------------
   //mcmc info
   //--------------------
   //methods
};



//--------------------------------------------------
//a single iteration of the MCMC for brt model
void mbrt::draw(rn& gen)
{
   //All the usual steps
   brt::draw(gen);

   // Update the in-sample predicted vector
   setf();

   // Update the in-sample residual vector
   setr();

}
//--------------------------------------------------
//draw theta for a single bottom node for the brt model
double mbrt::drawnodetheta(sinfo& si, rn& gen)
{
   msinfo& msi=static_cast<msinfo&>(si);
   double muhat = msi.sumwy/msi.sumw;
   double a = 1.0/(ci.tau*ci.tau);
   return (msi.sumw*muhat)/(a+msi.sumw) + gen.normal()/sqrt(a+msi.sumw);
}
//--------------------------------------------------
//lm: log of integrated likelihood, depends on prior and suff stats
double mbrt::lm(sinfo& si)
{
   msinfo& msi=static_cast<msinfo&>(si);
   double t2 =ci.tau*ci.tau;
   double k = msi.sumw*t2+1;
   return -.5*log(k)+.5*msi.sumwy*msi.sumwy*t2/k;
}
//--------------------------------------------------
//Add in an observation, this has to be changed for every model.
//Note that this may well depend on information in brt with our leading example
//being double *sigma in cinfo for the case of e~N(0,sigma_i^2).
// Note that we are using the training data and the brt object knows the training data
//     so all we need to specify is the row of the data (argument size_t i).
void mbrt::add_observation_to_suff(diterator& diter, sinfo& si)
{
   msinfo& msi=static_cast<msinfo&>(si);
   double w;
   w=1.0/(ci.sigma[*diter]*ci.sigma[*diter]);
   msi.n+=1;
   msi.sumw+=w;
   msi.sumwy+=w*diter.gety();
}
//--------------------------------------------------
//pr for brt
void mbrt::pr()
{
   COUT << "***** mbrt object:\n";
   COUT << "Conditioning info:" << endl;
   COUT << "   mean:   tau=" << ci.tau << endl;
   if(!ci.sigma)
     COUT << "         sigma=[]" << endl;
   else
     COUT << "         sigma=[" << ci.sigma[0] << ",...," << ci.sigma[di->n-1] << "]" << endl;
   brt::pr();
}

#endif
