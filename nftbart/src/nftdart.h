/*
 * Copyright (C) 2023 Mehri BagheriMohmadiPour 
 *  
 * This file is part of nftbart.
 * nftdart.h
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
 *
 */

#ifndef GUARD_nftdart_h
#define GUARD_nftdart_h

class nftdart {
public:
  nftdart() { }
  void setdart(void) { }
  void startdart(void) { }
  double gettheta(void) { return theta; }
   std::vector<size_t>& getnv() {return nv;}
   std::vector<double>& getpv() {return pv;}
   void setpv(double *varprob) {
     for(size_t j=0;j<p;j++) pv[j]=varprob[j];
   }
  //draw function
protected:
  int p;
bool dart,dartOn,aug,const_theta;
double a,b,rho,theta,omega;  
double *grp;
   std::vector<size_t> nv;
   std::vector<double> pv, lpv;
};
  
#endif
