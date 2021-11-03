/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         and Hugh A. Chipman
 *  
 * This file is part of nftbart.
 * dinfo.h
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

#ifndef GUARD_dinfo_h
#define GUARD_dinfo_h

//#ifdef _OPENMP
//#   include <omp.h>
//#endif

#include "rn.h"

class brt; //let compiler know this class will be defined so we can define brtMethodWrapper:
class brtMethodWrapper {
public:
  typedef double (brt::* brtMethodp)(size_t i);
  brtMethodWrapper(const brtMethodp ptomethod, brt& object) { p=ptomethod; o=&object; }
  //eg: instantiate as brtMethodWrapper mywrapper(&brt::f,bm);
  ~brtMethodWrapper() {}
  double callMethod(size_t i) { return ((o)->*(p))(i); }

  brtMethodp p;
  brt* o;
};

class dinfo {
public:
   dinfo() {p=0;n=0;x=0;y=0;tc=1;}
   dinfo(const dinfo& d) : p(d.p),n(d.n),x(d.x),y(d.y),tc(d.tc) {}
   dinfo(size_t ip, size_t in, double *ix, double *iy) : p(ip), n(in), x(ix), y(iy), tc(1) {}
   dinfo(size_t ip, size_t in, double *ix, double *iy,int itc) : p(ip), n(in), x(ix), y(iy), tc(itc) {}
   size_t p;  //number of vars
   size_t n;  //number of observations
   double *x; // jth var of ith obs is *(x + p*i+j)
   double *y; // ith y is *(y+i) or y[i]
   int tc; //thread count

   // compound addition operator
   dinfo& operator+=(const dinfo& rhs) {
      if(this->x==rhs.x) { //sanity check: they have to be dinfo's for the same data!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]+=rhs.y[i];
      }
      return *this;
   }
  // compound addition operator with vector rhs
   dinfo& operator+=(const std::vector<double>& rhs) {
      if(this->n==rhs.size()) { //sanity check: they have to be the same size!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]+=rhs[i];
      }
      return *this;
   }
  // compound addition operator with scalar rhs
   dinfo& operator+=(const double rhs) {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]+=rhs;
    return *this;
   }
   // compound addition opertor with function pointer rhs
   dinfo& operator+=(brtMethodWrapper& bmw)
   {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]+=bmw.callMethod(i);
      return *this;
   }
   // compound subtraction operator
   dinfo& operator-=(const dinfo& rhs) {
      if(this->x==rhs.x) { //sanity check: they have to be dinfo's for the same data!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]-=rhs.y[i];
      }
      return *this;
   }
  // compound subtraction operator with vector rhs
   dinfo& operator-=(const std::vector<double>& rhs) {
      if(this->n==rhs.size()) { //sanity check: they have to be the same size!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]-=rhs[i];
      }
      return *this;
   }
   // compound subtraction operator with scalar rhs
   dinfo& operator-=(const double rhs) {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]-=rhs;
    return *this;
   }
   // compound subtraction opertor with function pointer rhs
   dinfo& operator-=(brtMethodWrapper& bmw)
   {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]-=bmw.callMethod(i);
      return *this;
   }
   // compound multiplication operator
   dinfo& operator*=(const dinfo& rhs) {
      if(this->x==rhs.x) { //sanity check: they have to be dinfo's for the same data!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]*=rhs.y[i];
      }
      return *this;
   }
  // compound multiplication operator with vector rhs
   dinfo& operator*=(const std::vector<double>& rhs) {
      if(this->n==rhs.size()) { //sanity check: they have to be the same size!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]*=rhs[i];
      }
      return *this;
   }
  // compound multiplication operator with scalar rhs
   dinfo& operator*=(const double rhs) {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]*=rhs;
    return *this;
   }
   // compound multiplication opertor with function pointer rhs
   dinfo& operator*=(brtMethodWrapper& bmw)
   {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]*=bmw.callMethod(i);
      return *this;
   }
   // compound division operator
   dinfo& operator/=(const dinfo& rhs) {
      if(this->x==rhs.x) { //sanity check: they have to be dinfo's for the same data!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]/=rhs.y[i];
      }
      return *this;
   }
  // compound division operator with vector rhs
   dinfo& operator/=(const std::vector<double>& rhs) {
      if(this->n==rhs.size()) { //sanity check: they have to be the same size!
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]/=rhs[i];
      }
      return *this;
   }
  // compound division operator with scalar rhs
   dinfo& operator/=(const double rhs) {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]/=rhs;
    return *this;
   }
   // compound division opertor with function pointer rhs
   dinfo& operator/=(brtMethodWrapper& bmw)
   {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]/=bmw.callMethod(i);
      return *this;
   }

   // assignment operator
   dinfo& operator=(const dinfo& rhs) {
      if(&rhs != this) {
        this->n = rhs.n;
        this->p = rhs.p;
        this->x = rhs.x;
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]=rhs.y[i];     
      }
      return *this; 
   }
   // assignment operator with vector rhs
   dinfo& operator=(const std::vector<double>& rhs) {
      if(this->x && this->n == rhs.size()) {
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]=rhs[i];     
      }
      return *this; 
   }
   // assignment operator with scalar rhs
   dinfo& operator=(const double rhs) {
      if(this->x) {
        #ifdef _OPENMP
        #pragma omp parallel for num_threads(tc)
        #endif
        for(size_t i=0;i<n;i++)
          y[i]=rhs;     
      }
      return *this; 
   }
   // assignment operator with function pointer rhs
   dinfo& operator=(brtMethodWrapper& bmw)
   {
      #ifdef _OPENMP
      #pragma omp parallel for num_threads(tc)
      #endif
      for(size_t i=0;i<n;i++)
        y[i]=bmw.callMethod(i);
      return *this;
   }

};

// a dinfo iterator, of sorts.
class diterator : public std::iterator<std::input_iterator_tag, size_t>
{
	size_t i, end;
	dinfo di;
public:
  diterator(dinfo* d) : i(0),end((*d).n),di(*d) {}  //copy of d, helps OPENMP speed apparently.
  diterator(const diterator& dit) : i(dit.i),end(dit.end),di(dit.di) {}
  diterator(dinfo* d, size_t first, size_t last) : i(first),end(last),di(*d) {}
  diterator& operator++() {++i;return *this;}
  diterator operator++(int) {diterator tmp(*this); operator++(); return tmp;}
  double* getxp() { return di.x+i*di.p; }
  double getx() { return *(getxp()); }
  double* getyp() { return di.y+i; }
  double gety() { return *(getyp()); }
  void sety(double val) { di.y[i]=val; }
  size_t until() { return end; }
  bool operator==(const diterator& rhs) { return i==rhs.i; }
  bool operator==(size_t last) { return i==last; }
  bool operator!=(const diterator& rhs) { return i!=rhs.i; }
  bool operator!=(size_t last) { return i!=last; }
  bool operator<=(const diterator& rhs) { return i<=rhs.i; }
  bool operator<=(size_t last) { return i<=last; }
  bool operator<(const diterator& rhs) { return i<rhs.i; }
  bool operator<(size_t last) { return i<last; }
  bool operator>(const diterator& rhs) { return i>rhs.i; }
  bool operator>(size_t last) { return i>last; }
  size_t operator*() { return i; }
};


#endif
