#include <cmath>
#include <iostream>
#include "tree.h"
#include "info.h"
#include "rng.h" 

using std::cout;
using std::endl;

//pi and log(2*pi)
#define PI 3.1415926535897931
#define LTPI 1.83787706640934536

//--------------------------------------------------
//normal density
double pn(
   double x,    //variate
   double m,    //mean
   double v     //variance
);
//--------------------------------------------------
//normal cdf --ADDED!!!
double phi(double x);
//--------------------------------------------------
//draw from a discrete distribution
int rdisc(
   double *p,   //vector of probabilities
   RNG& gen     //random number generator
);
//--------------------------------------------------
//evaluate tree tr on grid xi, write to os
void grm(tree& tr, xinfo& xi, std::ostream& os);
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi, int* vartype);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, pinfo& pi, tree::npv& goodbots,int* vartype);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi, std::vector<size_t>& goodvars, int* vartype);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else a/(1+d)^b
double pgrow(tree::tree_p n, xinfo& xi, pinfo& pi,int* vartype);
//--------------------------------------------------
//get sufficients stats for all bottom nodes
void allsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<sinfo>& sv);
//--------------------------------------------------
//get sufficient stats for children (v,c) of node nx in tree x
void getsuff(tree& x, tree::tree_cp nx, size_t v, size_t c, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
//--------------------------------------------------
//get sufficient stats for pair of bottom children nl(left) and nr(right) in tree x
void getsuff(tree& x, tree::tree_cp nl, tree::tree_cp nr, xinfo& xi, dinfo& di, sinfo& sl, sinfo& sr);
//--------------------------------------------------
//log of the integreted likelihood
double lil(size_t n, double sy, double sy2, double sigma, double tau);
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, std::vector<double>& fv);
//--------------------------------------------------
//fit
void fit(tree& t, xinfo& xi, dinfo& di, double* fv);
//--------------------------------------------------
//partition
void partition(tree& t, xinfo& xi, dinfo& di, std::vector<size_t>& pv);
//--------------------------------------------------
// draw all the bottom node mu's
#ifdef MPIBART
void MPImasterdrmu(tree& t, xinfo& xi, pinfo& pi, RNG& gen, size_t numslaves);
#else
void drmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, RNG& gen);
#endif
//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi);
//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc,int* vartype);
//get min/max for p predictors needed to make cutpoints.
void makeminmax(size_t p, size_t n, double *x, std::vector<double> &minx, std::vector<double> &maxx);
//make xinfo = cutpoints given minx/maxx vectors
void makexinfominmax(size_t p, xinfo& xi, size_t nc, std::vector<double> &minx, std::vector<double> &maxx,int* vartype);
#ifdef MPIBART
//MPI Calls
void MPImasterallsuff(tree& x, tree::npv& bnv, std::vector<sinfo>& sv, size_t numslaves);
void MPIslaveallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv);
void MPIslavedrmu(tree& t, xinfo& xi, dinfo& di);
void MPImastergetsuff(tree::tree_cp nl, tree::tree_cp nr, sinfo &sl, sinfo &sr, size_t numslaves);
void MPImastergetsuffvc(tree::tree_cp nx, size_t v, size_t c, xinfo& xi, sinfo& sl, sinfo& sr, size_t numslaves);
void MPImastersendbirth(tree::tree_p nx, size_t v, size_t c, double mul, double mur, size_t numslaves);
void MPImastersenddeath(tree::tree_p nx, double mu, size_t numslaves);
void MPImastersendnobirthdeath(size_t numslaves);
void MPIslaveupdatebirthdeath(tree& x);
void MPIslavegetsuff(tree& x, xinfo& xi, dinfo& di);
void makepred(dinfo dip, xinfo &xi, std::vector<tree> &t, double *ppredmean);
void makeypred(dinfo dip, xinfo &xi, std::vector<tree> &t, double sigma, double *ppredmean);
void makepostpred(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp);
void makepostpred2(dinfo dip, xinfo &xi, std::vector< std::vector<tree> > &t, double *postp, double *postp2);
double logcalp(std::vector<double> &theta, dinfo dip, xinfo &xi, std::vector<tree> &t, double sigmae, double sigma, size_t pth, size_t myrank);
#endif
