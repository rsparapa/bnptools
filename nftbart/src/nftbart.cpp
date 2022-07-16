
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

//#include "rtnorm.h"
#include "rn.h"
#include "tree.h"
#include "treefuns.h"
#include "dinfo.h"
#include "brt.h"
//#include "brtfuns.h"
#include "mbrt.h"
#include "ambrt.h"
#include "sbrt.h"
#include "psbrt.h"
#include "DPMLIOmutau.h"
#include "DPMLIOneal8.h"
#include "cnft.h"
#include "cprnft.h"
#include "cpsambrt_predict.h"
#include "cnft2.h"
#include "cprnft2.h"
#include "cpsambrt_predict2.h"
