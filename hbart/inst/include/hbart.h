
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <map>
#include <vector>
#include <string>
#include <limits>

#ifndef NotInR
#include <Rcpp.h>
#define COUT Rcpp::Rcout 
#else
#define COUT std::cout 
using std::cout;
#endif

using std::endl;

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hbart/rn.h"
#include "hbart/tree.h"
#include "hbart/treefuns.h"
#include "hbart/dinfo.h"
#include "hbart/brt.h"
//#include "hbart/brtfuns.h"
#include "hbart/mbrt.h"
#include "hbart/ambrt.h"
#include "hbart/sbrt.h"
#include "hbart/psbrt.h"
/*
#define VEC(x, i) x[i]

double DPMLIOmutau_F(double y, double mu, double tau);

double DPMLIOmutau_G0tau(double a0, double b0, rn &gen);

double DPMLIOmutau_G0mu(double tau, double m0, double k0, rn &gen);

void DPMLIOmutau_P0(size_t row_c, SEXP _Y, SEXP _phi, 
		    double m0, double k0, double a0, double b0, rn &gen);

double DPMtau_F(double y, double tau);

double DPMtau_G0(double a0, double b0, rn &gen);

void DPMtau_P0(size_t row_c, SEXP _Y, SEXP _phi, 
		    double a0, double b0, rn &gen);

#include "hbart/DPMLIOneal8.h"
#include "hbart/DPMtauneal8.h"
*/
