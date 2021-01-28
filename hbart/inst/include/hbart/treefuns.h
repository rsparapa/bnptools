//     treefuns.h: BART tree class helper functions header.
//     Copyright (C) 2012-2016 Matthew T. Pratola, Robert E. McCulloch and Hugh A. Chipman
//
//     This file is part of hbart.
//
//     hbart is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 2 of the License, or
//     (at your option) any later version.
//
//     hbart is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//     Author contact information
//     Matthew T. Pratola: mpratola@gmail.com
//     Robert E. McCulloch: robert.e.mculloch@gmail.com
//     Hugh A. Chipman: hughchipman@gmail.com


#ifndef GUARD_treefuns_h
#define GUARD_treefuns_h

/*
#include <iostream>
#include "tree.h"

#include "treefuns.h"

#ifndef NotInR
#include <Rcpp.h>
#define COUT Rcpp::Rcout 
#else
#define COUT std::cout  << "std::cout "
#endif
*/

//--------------------------------------------------
//make xinfo which has nc cutpoints uniform on [0,1] for each x variable.
// this is in brt: void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
void makeUnifXinfo(size_t p,size_t nc,xinfo& xi);
//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi);
//--------------------------------------------------
//evaluate tree tr on grid xgrid, write to os
void grm(tree& tr, xinfo& xgrid, std::ostream& os);
//--------------------------------------------------
//fit tree at matrix of x, matrix is stacked columns x[i,j] is *(x+p*i+j)
void fit(tree& t, xinfo& xi, size_t p, size_t n, double *x,  double* fv);
//--------------------------------------------------
//does a (bottom) node have variables you can split on?
bool cansplit(tree::tree_p n, xinfo& xi);
//--------------------------------------------------
//find variables n can split on, put their indices in goodvars
void getgoodvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);
//--------------------------------------------------
// Get the L,U values for a node in the tree *given* the tree
// structure both above and below that node.
void getLU(tree::tree_p node, xinfo& xi, int* L, int* U);
//--------------------------------------------------
// Get the L,U values for a node in the tree *given* the tree
// structure both above and below that node.
void getvarLU(tree::tree_p node, size_t var, xinfo& xi, int* L, int* U);

//--------------------------------------------------
// These ones support the rotate code, but could be generally useful too.
//--------------------------------------------------
// Does the tree split on variable v at cutpoint c?
bool hasvcsplit(tree::tree_p t, size_t v, size_t c);
//--------------------------------------------------
// Does the node split on variable v?
bool splitsonv(tree::tree_p t, size_t v);
//--------------------------------------------------
// Do both nodes split on variable v?
bool splitsonv(tree::tree_p nl, tree::tree_p nr, size_t v);
//--------------------------------------------------
// Is this a leaf node?
bool isleaf(tree::tree_p t);
//--------------------------------------------------
// Are these two nodes equal?
bool arenodesequal(tree::tree_p nl, tree::tree_p nr);
//--------------------------------------------------
// Are these two nodes leaf nodes?
bool arenodesleafs(tree::tree_p nl, tree::tree_p nr);
//--------------------------------------------------
// number of available cutpoints at node n for variable var
int getnumcuts(tree::tree_p n, xinfo& xi, size_t var);
//--------------------------------------------------
// Find number of variables internal tree node n can split on
void getinternalvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars);
//--------------------------------------------------
// End of rotate helper functions
//--------------------------------------------------


//--------------------------------------------------
//make xinfo which has nc cutpoints uniform on [0,1] for each x variable.
//We put this in tree to make it easier to test tree methods without having any x variables
//which is helpful in tree/doc.cpp.  But usually in the model classes we are using the makexinfo
//method(s) in brt (or whatever model class).
void makeUnifXinfo(size_t p,size_t nc,xinfo& xi)
{
   double xinc = 1.0/(nc+1.0);
   xi.resize(p);
   for(size_t i=0;i<p;i++) xi[i].resize(nc);
   double cut;
   for(size_t j=0;j<nc;j++) {
      cut = 0.0 + (j+1)*xinc;
      for(size_t i=0;i<p;i++) xi[i][j]=cut;
   }
}
//--------------------------------------------------
//write cutpoint information to screen
void prxi(xinfo& xi)
{
   COUT << "xinfo: \n";
   for(size_t v=0;v!=xi.size();v++) {
      COUT << "v: " << v << std::endl;
      for(size_t j=0;j!=xi[v].size();j++) COUT << "j,xi[v][j]: " << j << ", " << xi[v][j] << std::endl;
   }
   COUT << "\n\n";
}
//--------------------------------------------------
//evalute tree tr on grid given by xgrid and write to os
void grm(tree& tr, xinfo& xgrid, std::ostream& os)
{
   size_t p = xgrid.size();
   if(p!=2) {
      COUT << "error in grm, p !=2\n";
      return;
   }
   size_t n1 = xgrid[0].size();
   size_t n2 = xgrid[1].size();
   tree::tree_p bp; //pointer to bottom node
   double *x = new double[2];
   for(size_t i=0;i!=n1;i++) {
      for(size_t j=0;j!=n2;j++) {
         x[0] = xgrid[0][i];
         x[1] = xgrid[1][j];
         bp = tr.bn(x,xgrid);
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
//--------------------------------------------------
// Get the L,U values for a node in the tree *given* the tree
// structure both above and below that node.
void getLU(tree::tree_p node, xinfo& xi, int* L, int* U)
{
   getvarLU(node,node->getv(),xi,L,U);
}
//--------------------------------------------------
// Get the L,U values for a node in the tree *given* the tree
// structure both above and below that node for a given variable var.
void getvarLU(tree::tree_p node, size_t var, xinfo& xi, int* L, int* U)
{
   tree::tree_p l,r;

   *L=0; *U = xi[var].size()-1;
   l=node->getl();
   r=node->getr();

   bool usel,user;
   usel=l->nuse(var);
   user=r->nuse(var);
   if(usel && user)
   {
      l->rl(var,L);
      r->ru(var,U);
   }
   else if(usel)
   {
      node->rg(var,L,U);
      l->rl(var,L);
   }
   else //handles both user case and !usel,!user case since ru() will do nothing for the latter.
   {
      node->rg(var,L,U);
      r->ru(var,U);
   }
}

//--------------------------------------------------
// These ones support the rotate code, but could be generally useful too.
//--------------------------------------------------
// Does the tree split on variable v at cutpoint c?
bool hasvcsplit(tree::tree_p t, size_t v, size_t c)
{
   tree::npv tnodes;

   tnodes.clear();
   t->getnodesonvc(tnodes,v,c);
   if(tnodes.size())
      return true;
   return false;
}
//--------------------------------------------------
// Does the node split on variable v?
bool splitsonv(tree::tree_p t, size_t v)
{
   if(!t->getl()) //its a leaf
      return false;
   if(t->getv()==v) return true;
   return false;
}
//--------------------------------------------------
// Do both nodes split on variable v?
bool splitsonv(tree::tree_p nl, tree::tree_p nr, size_t v)
{
   if(!nl->getl() || !nr->getl()) //one is a leaf
      return false;
   if(nl->getv()==v && nr->getv()==v) return true;
   return false;
}
//--------------------------------------------------
// Is this a leaf node?
bool isleaf(tree::tree_p t)
{
   if(!t->getl()) return true;
   return false;
}
//--------------------------------------------------
// Are these two nodes equal?
bool arenodesequal(tree::tree_p nl, tree::tree_p nr)
{
   if(!nl->getl() || !nr->getl())  //one is a leaf
      return false;
   if(nl->getv()==nr->getv() && nl->getc()==nr->getc()) //both split on v at c
      return true;
   return false;
}
//--------------------------------------------------
// Are these two nodes leaf nodes?
bool arenodesleafs(tree::tree_p nl, tree::tree_p nr)
{
   if(!nl->getl() && !nr->getl()) //both are leaves
      return true;
   return false;
}
//--------------------------------------------------
// number of available cutpoints at node n for variable var
int getnumcuts(tree::tree_p n, xinfo& xi, size_t var)
{
   int L,U;

   getvarLU(n,var,xi,&L,&U);

   return std::max(0,U-L+1);
}
//--------------------------------------------------
// Find number of variables internal tree node n can split on
void getinternalvars(tree::tree_p n, xinfo& xi,  std::vector<size_t>& goodvars)
{
   int L,U;

   for(size_t v=0;v!=xi.size();v++) {//try each variable
      L=0; U = xi[v].size()-1;
      getvarLU(n,v,xi,&L,&U);
      if(U>=L) goodvars.push_back(v);
   }
}
//--------------------------------------------------
// End of rotate helper functions
//--------------------------------------------------

#endif
