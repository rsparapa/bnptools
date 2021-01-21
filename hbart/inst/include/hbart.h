
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
#include "hbart/rrn.h"
#include "hbart/tree.h"
#include "hbart/treefuns.h"

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, double pipb, tree::npv& goodbots);
//--------------------------------------------------
/*
//bprop: function to generate birth proposal
void bprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, size_t& v, size_t& c, double& pr, rn& gen);
//--------------------------------------------------
// death proposal
void dprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, double& pr, rn& gen);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, brt::tprior& tp);
*/
//--------------------------------------------------
//calculate beginning and end points of data vector to be accessed in parallel computations
void calcbegend(int n, int my_rank, int thread_count, int* beg, int* end);

//--------------------------------------------------
// Functions to support change-of-variable proposal
//--------------------------------------------------
// update the correlation matrix for chgv move taking into account that not all
// variables may be eligible at pertnode.
void updatecormat(tree::tree_p pertnode, xinfo& xi, std::vector<std::vector<double> >& chgv);
//--------------------------------------------------
// renormalize the correlation matrix so that the probability of row sums to 1.
void normchgvrow(size_t row, std::vector<std::vector<double> >& chgv);
//--------------------------------------------------
// randomly choose a new variable to transition to from oldv
size_t getchgv(size_t oldv, std::vector<std::vector<double> >& chgv, rn& gen);

//--------------------------------------------------
// Functions to support rotate proposal
//--------------------------------------------------
//setup the initial right rotation
void rotright(tree::tree_p n);
//--------------------------------------------------
//setup the initial left rotation
void rotleft(tree::tree_p n);
//--------------------------------------------------
//eliminate immediate dead ends from the rotate
void reduceleft(tree::tree_p n, size_t v, size_t c);
//--------------------------------------------------
//eliminate immediate dead ends from the rotate
void reduceright(tree::tree_p n, size_t v, size_t c);
//--------------------------------------------------
//split tree along variable v at cutpoint c retaining only 
//part of the tree that is ``left'' of this v,c rule
void splitleft(tree::tree_p t, size_t v, size_t c);
//--------------------------------------------------
//split tree along variable v at cutpoint c retaining only 
//part of the tree that is ``right'' of this v,c rule
void splitright(tree::tree_p t, size_t v, size_t c);
//--------------------------------------------------
//does an actual merge (randomly chosen) 
bool merge(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c, rn& gen);
//--------------------------------------------------
// only to get nways, not to actually do the merge.
bool mergecount(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways);
//--------------------------------------------------
// End of functions to support rotate proposal
//--------------------------------------------------

#include "hbart/dinfo.h"
#include "hbart/brt.h"
#include "hbart/brtfuns.h"
#include "hbart/mbrt.h"
#include "hbart/ambrt.h"
#include "hbart/sbrt.h"
#include "hbart/psbrt.h"
