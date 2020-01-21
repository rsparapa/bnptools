/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017-2018 Robert McCulloch, Rodney Sparapani
 *                          and Charles Spanbauer
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include <ctime>
#include "info.h"

class dart {
public:
   //------------------------------
   //friends
   //------------------------------
   //constructor/destructor
   dart();
   ~dart();
   //------------------------------
   //operators
   dart& operator=(const dart&);
   //------------------------------
   //setters and getters
   //------------------------------
   void setdart(size_t _p, double _theta, size_t _thetaDraw, 
		double _a, double _b, double _rho){
     this->p=_p; this->thetaDraw=_thetaDraw;
     this->a=_a;this->b=_b;this->rho=_rho;
     if(_theta==0.) {
       this->const_theta=false;
       this->theta=1.;
     }
     else {
       this->const_theta=true;
       this->theta=_theta;
     }
     for(size_t j=0;j<p;j++) {
       pv.push_back(1/(double)p);
       lpv.push_back(-std::log((double)p));
       nv.push_back(0);
     }
   } 
   std::vector<double>& getlpv() {return lpv;}
   std::vector<double>& getpv() {return pv;}
   std::vector<size_t>& getnv() {return nv;}
   double gettheta() {return theta;}
   void setnv(std::vector<size_t> _nv) {this->nv=_nv;}
   //--------------------------------------------------
   //draw variable splitting probabilities from Dirichlet (Linero, 2018)
   void draw_s(rn& gen);
   //--------------------------------------------------
   //draw Dirichlet sparsity parameter from posterior using grid
   void draw_theta0(rn& gen);   
protected:
   bool const_theta;
   double a,b,rho,theta,omega;
   size_t p,thetaDraw;
   std::vector<size_t> nv;
   std::vector<double> pv,lpv;
};
