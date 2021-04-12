
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
