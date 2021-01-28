// rrn.h: Random number generator class for using BART in R package.
// Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
//
// This file is part of hbart.
//
// hbart is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// hbart is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author contact information
// Matthew T. Pratola: mpratola@gmail.com
// Robert E. McCulloch: robert.e.mculloch@gmail.com
// Hugh A. Chipman: hughchipman@gmail.com


#ifndef RRN_H
#define RRN_H

#ifdef Rcpp_hpp
using R::rchisq;
#else
//extern "C" {
#include <R.h>
#include <Rmath.h>
//};
#endif

//#include "rn.h"

class rrn: public rn
{
public:
//constructor
   rrn():df(1) {}
//virtual
   virtual ~rrn() {}
   virtual double normal() {return norm_rand();}
   virtual double uniform() { return unif_rand();}
   virtual double chi_square() {return rchisq((double)df);}
   virtual double exp() {return exp_rand();}
   virtual void set_df(int df) {this->df=df;}
//get,set
   int get_df() {return df;}
private:
   int df;
};

#endif
