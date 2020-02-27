#ifndef GUARD_QFUNS_h
#define GUARD_QFUNS_h

#include "BART3/tree.h"
#include "BART3/treefuns.h"
#include "BART3/info.h"

//compute n and \sum y_i for left and right give bot and v,c
void qgetsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);

//compute n and \sum y_i for left and right bots
void qgetsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr);

//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void qallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<size_t>& nv, std::vector<double>& syv);

// draw all the bottom node mu's
void qdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen);


//compute n and \sum y_i for left and right give bot and v,c
void getsuff(tree& x, tree::tree_p nx, size_t v, size_t c, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr)
{
   double *xx;//current x
   nl=0; syl=0.0;
   nr=0; syr=0.0;

   for(size_t i=0;i<di.n;i++) if(di.q == 1) {
      xx = di.x + i*di.p;
      if(nx==x.bn(xx,xi)) { //does the bottom node = xx's bottom node
         if(xx[v] < xi[v][c]) {
               nl++;
               syl += di.y[i];
          } else {
               nr++;
               syr += di.y[i];
          }
      }
   }
}


//compute n and \sum y_i for left and right bots
void qgetsuff(tree& x, tree::tree_p l, tree::tree_p r, xinfo& xi, dinfo& di, size_t& nl, double& syl, size_t& nr, double& syr)
{
   double *xx;//current x
   nl=0; syl=0.0;
   nr=0; syr=0.0;

   for(size_t i=0;i<di.n;i++) if(di.q == 1) {
      xx = di.x + i*di.p;
      tree::tree_cp bn = x.bn(xx,xi);
      if(bn==l) {
         nl++;
         syl += di.y[i];
      }
      if(bn==r) {
         nr++;
         syr += di.y[i];
      }
   }
}


//get sufficients stats for all bottom nodes, this way just loop through all the data once.
void qallsuff(tree& x, xinfo& xi, dinfo& di, tree::npv& bnv, std::vector<size_t>& nv, std::vector<double>& syv)
{
   tree::tree_cp tbn; //the pointer to the bottom node for the current observations
   size_t ni;         //the  index into vector of the current bottom node
   double *xx;        //current x

   bnv.clear();
   x.getbots(bnv);

   typedef tree::npv::size_type bvsz;
   bvsz nb = bnv.size();
   nv.resize(nb);
   syv.resize(nb);

   std::map<tree::tree_cp,size_t> bnmap;
   for(bvsz i=0;i!=bnv.size();i++) {bnmap[bnv[i]]=i;nv[i]=0;syv[i]=0.0;}

   for(size_t i=0;i<di.n;i++) if(di.q == 1) {
      xx = di.x + i*di.p;
      tbn = x.bn(xx,xi);
      ni = bnmap[tbn];

      ++(nv[ni]);
      syv[ni] += di.y[i];
   }
}


// draw all the bottom node mu's
void qdrmu(tree& t, xinfo& xi, dinfo& di, pinfo& pi, double sigma, rn& gen)
{
   tree::npv bnv;
   std::vector<size_t> nv;
   std::vector<double> syv;
   qallsuff(t,xi,di,bnv,nv,syv);

   for(tree::npv::size_type i=0;i!=bnv.size();i++) 
      bnv[i]->settheta(drawnodemu(nv[i],syv[i],pi.tau,sigma,gen));
}

#endif
