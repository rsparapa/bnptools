/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         and Hugh A. Chipman
 *  
 * This file is part of nftbart.
 * sbrt.h
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

#ifndef GUARD_sbrt_h
#define GUARD_sbrt_h

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

class ssinfo : public sinfo { //sufficient statistics (will depend on end node model)
public:
   ssinfo():sinfo(),sumy2(0.0) {}
   ssinfo(const ssinfo& is):sinfo(is),sumy2(is.sumy2) {}
   virtual ~ssinfo() {}  //need this so memory is properly freed in derived classes.
   double sumy2;
   // compound addition operator needed when adding suff stats
   virtual sinfo& operator+=(const sinfo& rhs) {
      sinfo::operator+=(rhs);
      const ssinfo& srhs=static_cast<const ssinfo&>(rhs);
      sumy2+=srhs.sumy2;
      return *this;
   }
   // assignment operator for suff stats
   virtual sinfo& operator=(const sinfo& rhs)
   {
      if(&rhs != this) {
         sinfo::operator=(rhs);
         const ssinfo& srhs=static_cast<const ssinfo&>(rhs);
         this->sumy2 = srhs.sumy2;
      }
      return *this;
   }
   // addition opertor is defined in terms of compound addition
   const ssinfo operator+(const ssinfo& other) const {
      ssinfo result = *this; //copy of myself.
      result += other;
      return result;
   }
};

class sbrt : public brt 
{
public:
   //--------------------
   //classes
   // tprior and mcmcinfo are same as in brt
   class cinfo { //parameters for end node model prior
   public:
      cinfo():nu(1.0),lambda(1.0) {}
      double nu; //dof of variance prior
      double lambda; //scale of variance prior
   };
   //--------------------
   //constructors/destructors
   sbrt():brt() { t.settheta(1.0); }
   sbrt(double itheta):brt() { t.settheta(itheta); }
   //--------------------
   //methods
   void draw(rn& gen);
   void setci(double nu, double lambda) { ci.nu=nu; ci.lambda=lambda; }
   virtual double drawnodetheta(sinfo& si, rn& gen);
   virtual double lm(sinfo& si);
   virtual void add_observation_to_suff(diterator& diter, sinfo& si);
   virtual sinfo* newsinfo() { return new ssinfo; }
   virtual std::vector<sinfo*>& newsinfovec() { std::vector<sinfo*>* si= new std::vector<sinfo*>; return *si; }
   virtual std::vector<sinfo*>& newsinfovec(size_t dim) { std::vector<sinfo*>* si = new std::vector<sinfo*>; si->resize(dim); for(size_t i=0;i<dim;i++) si->push_back(new ssinfo); return *si; }
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
   double logam(double x);

};


//--------------------------------------------------
//a single iteration of the MCMC for brt model
void sbrt::draw(rn& gen)
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
double sbrt::drawnodetheta(sinfo& si, rn& gen)
{
   ssinfo& ssi=static_cast<ssinfo&>(si);
   int nupost=ssi.n+(int)ci.nu;
   double nulampost=ci.nu*ci.lambda+ssi.sumy2;
   gen.set_df(nupost);

   return sqrt((nulampost)/gen.chi_square());
}
//--------------------------------------------------
//lm: log of integrated likelihood, depends on prior and suff stats
double sbrt::lm(sinfo& si)
{
   ssinfo& ssi=static_cast<ssinfo&>(si);
   double val;
   double nstar;
   double nudiv2;

   nudiv2=ci.nu/2.0;
   nstar=(ci.nu+ssi.n)/2.0;
   val = nudiv2*log(ci.nu*ci.lambda);
   val+= -(ci.nu+ssi.n)/2.0*log(ci.nu*ci.lambda+ssi.sumy2);
   val+= logam(nstar)-logam(nudiv2);

   return val;
}
//--------------------------------------------------
//Add in an observation, this has to be changed for every model.
//Note that this may well depend on information in brt with our leading example
//being double *sigma in cinfo for the case of e~N(0,sigma_i^2).
// Note that we are using the training data and the brt object knows the training data
//     so all we need to specify is the row of the data (argument size_t i).
void sbrt::add_observation_to_suff(diterator& diter, sinfo& si)
{
   ssinfo& ssi=static_cast<ssinfo&>(si);
   ssi.n+=1;
   ssi.sumy2+=diter.gety()*diter.gety();
}
//--------------------------------------------------
//pr for brt
void sbrt::pr()
{
   COUT << "***** sbrt object:\n";
   COUT << "Conditioning info:" << endl;
   COUT << "      dof:  nu=" << ci.nu << endl;
   COUT << "    scale:  lambda=" << ci.lambda << endl;
   brt::pr();
}
//--------------------------------------------------
// compute the logarithm of the Gamma function
// people.sc.fsu.edu/~jburkardt/cpp_src/toms291/toms291.html
double sbrt::logam (double x)
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
