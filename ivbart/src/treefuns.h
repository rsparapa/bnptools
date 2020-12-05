#ifndef GUARD_treefuns_h
#define GUARD_treefuns_h

#include <iostream>
#include "tree.h"

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi);
//--------------------------------------------------
//evaluate tree tr on grid xi, write to os
void grm(tree& tr, xinfo& xi, std::ostream& os);
//--------------------------------------------------
//fit tree at matrix of x, matrix is stacked columns x[i,j] is *(x+p*i+j)
void fit(tree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv);
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);

//#include "treefuns.h"

//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
   cout << "xinfo: \n";
   for(size_t v=0;v!=xi.size();v++) {
      cout << "v: " << v << std::endl;
      for(size_t j=0;j!=xi[v].size();j++) cout << "j,xi[v][j]: " << j << ", " << xi[v][j] << std::endl;
   }
   cout << "\n\n";
}
//--------------------------------------------------
//evalute tree tr on grid given by xi and write to os
void grm(tree& tr, xinfo& xi, std::ostream& os)
{
   size_t p = xi.size();
   if(p!=2) {
      cout << "error in grm, p !=2\n";
      return;
   }
   size_t n1 = xi[0].size();
   size_t n2 = xi[1].size();
   tree::tree_p bp; //pointer to bottom node
   double *x = new double[2];
   for(size_t i=0;i!=n1;i++) {
      for(size_t j=0;j!=n2;j++) {
         x[0] = xi[0][i];
         x[1] = xi[1][j];
         bp = tr.bn(x,xi);
         os << x[0] << " " << x[1] << " " << bp->gettheta() << " " << bp->nid() << std::endl;
      }
   }
   delete[] x;
}
//--------------------------------------------------
//fit tree at matrix of x, matrix is stacked columns x[i,j] is *(x+p*i+j)
void fit(tree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv)
{
   tree::tree_p bn;
   for(size_t i=0;i<n;i++) {
      bn = t.bn(x+i*p,xi);
      fv[i] = bn->gettheta();
   }
}
//--------------------------------------------------
//does this bottom node n have any variables it can split on.
bool cansplit(tree::tree_p n, xinfo& xi)
{
   int L,U;
   bool v_found = false; //have you found a variable you can split on
   size_t v=0;
   while(!v_found && (v < xi.size())) { //invar: splitvar not found, vars left
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) v_found=true;
      v++;
   }
   return v_found;
}
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
   goodvars.clear();
   int L,U;
   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      n->rg(v,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}

#endif
