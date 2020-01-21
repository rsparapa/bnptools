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

#ifndef GUARD_bart_h
#define GUARD_bart_h

#include <ctime>

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "dart.h"

class bart {
public:
   //------------------------------
   //friends
   friend bool bd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double sigma,
		  dart& drt, bool aug, rn& gen);
   //------------------------------
   //constructor/destructor
   bart();
   bart(size_t m);
   bart(const bart&);
   ~bart();
   //------------------------------
   //operators
   bart& operator=(const bart&);
   //------------------------------
   //get,set
   size_t getm() {return m;}
   void setm(size_t m);
   void setdata(size_t p, size_t n, double *x, double *y, 
		   double theta, size_t thetaDraw,
		   double a, double b, double rho, size_t numcut);
   void setdata(size_t p, size_t n, double *x, double *y, 
		double theta, size_t thetaDraw,
		double a, double b, double rho, int* nc);
   void setpi(pinfo& pi) {this->pi = pi;}
   void setprior(double alpha, double beta, double tau)
      {pi.alpha=alpha; pi.mybeta = beta; pi.tau=tau;}
   void startdart() {this->useDart=true;}
   void settau(double tau) {pi.tau=tau;}
   tree& gettree(size_t i ) { return t[i];}
   xinfo& getxinfo() {return xi;}
   void setxinfo(xinfo& _xi);
   dart getdart() {return drt;}
   //------------------------------
   //public methods
   void birth(size_t i, size_t nid,size_t v, size_t c, double ml, double mr)
         {t[i].birth(nid,v,c,ml,mr);}
   void death(size_t i,size_t nid, double mu)
         {t[i].death(nid,mu);}
   void pr();
   void tonull() {for(size_t i=0;i!=t.size();i++) t[i].tonull();}
   void predict(size_t p, size_t n, double *x, double *fp);
   void draw(double sigma, rn& gen);
//   void draw_s(rn& gen);
   double f(size_t i) {return allfit[i];}
protected:
   dart drt;
   bool useDart;
   size_t m;  //number of trees
   std::vector<tree> t; //the trees
   pinfo pi; //prior and mcmc info
   //data
   size_t p,n; //x has dim p, n obserations
   double *x,*y;  //x is column stack, pxn
   xinfo xi; //cutpoint info
   //working
   double *allfit; //if the data is set, should be f(x)
   double *r;
   double *ftemp;
   dinfo di;
   std::vector<double> pv;
   std::vector<size_t> nv;
};

#endif
