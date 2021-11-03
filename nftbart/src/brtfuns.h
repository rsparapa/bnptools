/*
 * Copyright (C) 2012-2021 Matthew T. Pratola, Robert E. McCulloch,
 *                         and Hugh A. Chipman
 *  
 * This file is part of nftbart.
 * brtfuns.h
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

#ifndef GUARD_brtfuns_h
#define GUARD_brtfuns_h

/*
#include <iostream>
#include "tree.h"
#include "treefuns.h"
#include "brt.h"

using std::cout;
using std::endl;

#ifndef NotInR
#include <Rcpp.h>
#define COUT Rcpp::Rcout 
#else
#define COUT std::cout  << "std::cout "
#endif

#include "brtfuns.h"
*/

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc);
//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, double pipb, tree::npv& goodbots);
//--------------------------------------------------
//bprop: function to generate birth proposal
void bprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, size_t& v, size_t& c, double& pr, rn& gen);
//--------------------------------------------------
// death proposal
void dprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, double& pr, rn& gen);
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, brt::tprior& tp);
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

//--------------------------------------------------
//make xinfo = cutpoints
void makexinfo(size_t p, size_t n, double *x, xinfo& xi, size_t nc)
{
   double xinc;

   //compute min and max for each x
   std::vector<double> minx(p,INFINITY);
   std::vector<double> maxx(p,-INFINITY);
   double xx;
   for(size_t i=0;i<p;i++) {
      for(size_t j=0;j<n;j++) {
         xx = *(x+p*j+i);
         if(xx < minx[i]) minx[i]=xx;
         if(xx > maxx[i]) maxx[i]=xx;
      }
   }
   //make grid of nc cutpoints between min and max for each x.
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
      xinc = (maxx[i]-minx[i])/(nc+1.0);
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = minx[i] + (j+1)*xinc;
   }
}

//--------------------------------------------------
//compute prob of a birth, goodbots will contain all the good bottom nodes
double getpb(tree& t, xinfo& xi, double pipb, tree::npv& goodbots)
{
   double pb;  //prob of birth to be returned
   tree::npv bnv; //all the bottom nodes
   t.getbots(bnv);
   for(size_t i=0;i!=bnv.size();i++)
      if(cansplit(bnv[i],xi)) goodbots.push_back(bnv[i]);
   if(goodbots.size()==0) { //are there any bottom nodes you can split on?
      pb=0.0;
   } else {
      if(t.treesize()==1) pb=1.0; //is there just one node?
      else pb=pipb;
   }
   return pb;
}
//--------------------------------------------------
//bprop: function to generate birth proposal
void bprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, size_t& v, size_t& c, double& pr, rn& gen)
{

      //draw bottom node, choose node index ni from list in goodbots
      size_t ni = floor(gen.uniform()*goodbots.size());
      nx = goodbots[ni]; //the bottom node we might birth at

      //draw v,  the variable
      std::vector<size_t> goodvars; //variables nx can split on
      getgoodvars(nx,xi,goodvars);
      size_t vi = floor(gen.uniform()*goodvars.size()); //index of chosen split variable
      v = goodvars[vi];

      //draw c, the cutpoint
      int L,U;
      L=0; U = xi[v].size()-1;
      nx->rg(v,&L,&U);
      c = L + floor(gen.uniform()*(U-L+1)); //U-L+1 is number of available split points

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double Pbotx = 1.0/goodbots.size(); //proposal prob of choosing nx
      size_t dnx = nx->depth();
      double PGnx = tp.alpha/pow(1.0 + dnx,tp.beta); //prior prob of growing at nx

      double PGly, PGry; //prior probs of growing at new children (l and r) of proposal
      if(goodvars.size()>1) { //know there are variables we could split l and r on
         PGly = tp.alpha/pow(1.0 + dnx+1.0,tp.beta); //depth of new nodes would be one more
         PGry = PGly;
      } else { //only had v to work with, if it is exhausted at either child need PG=0
         if((int)(c-1)<L) { //v exhausted in new left child l, new upper limit would be c-1
            PGly = 0.0;
         } else {
            PGly = tp.alpha/pow(1.0 + dnx+1.0,tp.beta);
         }
         if(U < (int)(c+1)) { //v exhausted in new right child r, new lower limit would be c+1
            PGry = 0.0;
         } else {
            PGry = tp.alpha/pow(1.0 + dnx+1.0,tp.beta);
         }
      }

      double PDy; //prob of proposing death at y
      if(goodbots.size()>1) { //can birth at y because splittable nodes left
         PDy = 1.0 - pb;
      } else { //nx was the only node you could split on
         if((PGry==0) && (PGly==0)) { //cannot birth at y
            PDy=1.0;
         } else { //y can birth at either l or r
            PDy = 1.0 - pb;
         }
      }

      double Pnogy; //death prob of choosing the nog node at y
      size_t nnogs = x.nnogs();
      tree::tree_p nxp = nx->getp();
      if(nxp==0) { //no parent, nx is the top and only node
         Pnogy=1.0;
      } else {
         if(nxp->ntype() == 'n') { //if parent is a nog, number of nogs same at x and y
            Pnogy = 1.0/nnogs;
         } else { //if parent is not a nog, y has one more nog.
           Pnogy = 1.0/(nnogs+1.0);
         }
      }

      pr = (PGnx*(1.0-PGly)*(1.0-PGry)*PDy*Pnogy)/((1.0-PGnx)*Pbotx*PBx);
}
//--------------------------------------------------
// death proposal
void dprop(tree& x, xinfo& xi, brt::tprior& tp, double pb, tree::npv& goodbots, double& PBx, tree::tree_p& nx, double& pr, rn& gen)
{
      //draw nog node, any nog node is a possibility
      tree::npv nognds; //nog nodes
      x.getnogs(nognds);
      size_t ni = floor(gen.uniform()*nognds.size());
      nx = nognds[ni]; //the nog node we might kill children at

      //--------------------------------------------------
      //compute things needed for metropolis ratio

      double PGny; //prob the nog node grows
      size_t dny = nx->depth();
      PGny = tp.alpha/pow(1.0+dny,tp.beta);

      //better way to code these two?
      double PGlx = pgrow(nx->getl(),xi,tp);
      double PGrx = pgrow(nx->getr(),xi,tp);

      double PBy;  //prob of birth move at y
      if(nx->ntype()=='t') { //is the nog node nx the top node
         PBy = 1.0;
      } else {
         PBy = pb;
      }

      double Pboty;  //prob of choosing the nog as bot to split on when y
      int ngood = goodbots.size();
      if(cansplit(nx->getl(),xi)) --ngood; //if can split at left child, lose this one
      if(cansplit(nx->getr(),xi)) --ngood; //if can split at right child, lose this one
      ++ngood;  //know you can split at nx
      Pboty=1.0/ngood;

      double PDx = 1.0-PBx; //prob of a death step at x
      double Pnogx = 1.0/nognds.size();

      pr =  ((1.0-PGny)*PBy*Pboty)/(PGny*(1.0-PGlx)*(1.0-PGrx)*PDx*Pnogx);
}
//--------------------------------------------------
//get prob a node grows, 0 if no good vars, else alpha/(1+d)^beta
double pgrow(tree::tree_p n, xinfo& xi, brt::tprior& tp)
{
   if(cansplit(n,xi)) {
      return tp.alpha/pow(1.0+n->depth(),tp.beta);
   } else {
      return 0.0;
   }
}
//--------------------------------------------------
//calculate beginning and end points of data vector to be accessed in parallel computations
void calcbegend(int n, int my_rank, int thread_count, int* beg, int* end)
{
   if(n>=thread_count) {
      int h = n/thread_count;
      *beg = my_rank*h;
      *end = *beg+h;
      if(my_rank==(thread_count-1)) *end=n;
   }
   else // n < thread_count
   {
      *beg=my_rank;
      *end=my_rank+1;
      if(my_rank>=n) {
         *beg=0;
         *end=0;
      }
   }
}
//--------------------------------------------------
// Functions to support change-of-variable proposal
//--------------------------------------------------
// update the correlation matrix for chgv move taking into account that not all
// variables may be eligible at pertnode.
void updatecormat(tree::tree_p pertnode, xinfo& xi, std::vector<std::vector<double> >& chgv)
{
   int Ln,Un; //L,U for the ``new'' variable
   size_t oldv=pertnode->getv();
   size_t p=chgv.size();

   for(size_t i=0;i<p;i++) {
      if(i!=oldv && std::abs(chgv[oldv][i])>0.0) {
         if(chgv[oldv][i]<0.0)  //swap left,right branches
            pertnode->swaplr();
         getvarLU(pertnode,i,xi,&Ln,&Un);
         if(chgv[oldv][i]<0.0)  //undo the swap
            pertnode->swaplr();
         if(Un<Ln) //we can't transition to variable i here according to the tree structure
            chgv[oldv][i]=0.0;
      }
   }
}
//--------------------------------------------------
// renormalize the correlation matrix so that the probability of row sums to 1.
void normchgvrow(size_t row, std::vector<std::vector<double> >& chgv)
{
   double tmp=0.0;
   size_t p=chgv.size();

   for(size_t i=0;i<p;i++)
      tmp+=std::abs(chgv[row][i]);
   for(size_t i=0;i<p;i++)
      chgv[row][i]/=tmp;
}
//--------------------------------------------------
// randomly choose a new variable to transition to from oldv
size_t getchgv(size_t oldv, std::vector<std::vector<double> >& chgv, rn& gen)
{
   double cp=gen.uniform();
   size_t p=chgv.size();
   size_t newv=oldv;
   std::vector<double> cumprob;

   cumprob=chgv[oldv];
   cumprob[1]=std::abs(cumprob[1]);
   for(size_t i=1;i<p;i++)
      cumprob[i]=std::abs(cumprob[i])+cumprob[i-1];

   for(size_t i=0;i<p;i++) {
      if(cumprob[i] >= cp) { //propose transition to this variable
         newv=i;
         i=p;  //break out of loop early
      }        
   }
   return newv;
}

//--------------------------------------------------
// Functions to support rotate proposal
//--------------------------------------------------
// Rotate a given rotatable node in the tree.
// node n must be the left child of n->p and also not a root or terminal leaf node.
void rotright(tree::tree_p n)
{
   tree::tree_cp ctstar;
   tree::tree_p tstar;
   tree::tree_p newt = new tree;
   size_t vprime;
   size_t cprime;

   tstar=n->p->r;
   ctstar=n->p->r;
   tree::tree_p newnr = new tree(*ctstar);
   tstar->p=0;

   vprime=n->v;
   cprime=n->c;   
   n->v=n->p->v;
   n->c=n->p->c;
   n->p->v=vprime;
   n->p->c=cprime;
   newt->v=n->v;
   newt->c=n->c;

   newt->p=n->p;
   newt->l=n->r;
   newt->l->p=newt;
   newt->r=tstar;
   tstar->p=newt;
   n->p->r=newt;
   
   n->r=newnr;
   newnr->p=n;
}
//--------------------------------------------------
// Rotate a given rotatable node in the tree.
// node n must be the right child of n->p and also not a root or terminal leaf node.
void rotleft(tree::tree_p n)
{
   tree::tree_cp ctstar;
   tree::tree_p tstar;
   tree::tree_p newt = new tree;
   size_t vprime;
   size_t cprime;

   tstar=n->p->l;
   ctstar=n->p->l;
   tree::tree_p newnl = new tree(*ctstar);
   tstar->p=0;

   vprime=n->v;
   cprime=n->c;
   n->v=n->p->v;
   n->c=n->p->c;
   n->p->v=vprime;
   n->p->c=cprime;
   newt->v=n->v;
   newt->c=n->c;

   newt->p=n->p;
   newt->r=n->l;
   newt->r->p=newt;
   newt->l=tstar;
   tstar->p=newt;
   n->p->l=newt;

   n->l=newnl;
   newnl->p=n;
}
//--------------------------------------------------
// reduce the left sub-tree of the node that was rotated to the top
void reduceleft(tree::tree_p n, size_t v, size_t c)
{
   tree::tree_p temp;

   if(n->r->l && n->r->v==v) //right is not terminal and splits on v
      if(n->r->c >= c) { //then only keep left branch
         delete n->r->r;
         temp=n->r;
         n->r=temp->l;
         temp->l->p=n;
         temp->r=0;
         temp->l=0;
         temp->p=0;
         delete temp;
      }
   if(n->l->l && n->l->v==v) //left is not terminal and splits on v
      if(n->l->c >= c) { // then only keep left branch
         delete n->l->r;
         temp=n->l;
         n->l=temp->l;
         temp->l->p=n;
         temp->r=0;
         temp->l=0;
         temp->p=0;
         delete temp;
      }
}
//--------------------------------------------------
// reduce the right sub-tree of the node that was rotated to the top
void reduceright(tree::tree_p n, size_t v, size_t c)
{
   tree::tree_p temp;

   if(n->r->v==v && n->r->l) //right is not terminal and splits on v
      if(n->r->c <= c) { //then only keep right branch
         delete n->r->l;
         temp=n->r;
         n->r=temp->r;
         temp->r->p=n;
         temp->r=0;
         temp->l=0;
         temp->p=0;
         delete temp;
      }
   if(n->l->v==v && n->l->l) //left is not terminal and splits on v
      if(n->l->c <= c) { // then only keep right branch
         delete n->l->l;
         temp=n->l;
         n->l=temp->r;
         temp->r->p=n;
         temp->r=0;
         temp->l=0;
         temp->p=0;
         delete temp;
      }
}
//--------------------------------------------------
//split tree along variable v at cutpoint c retaining only 
//part of the tree that is ``left'' of this v,c rule
void splitleft(tree::tree_p t, size_t v, size_t c)
{
   tree::tree_p temp;

   if(t->l) //not terminal node
   {
      if(t->v==v && t->c >= c)
      {
         temp=t->l;
         if(t->isleft())
         {
            t->p->l=temp;
            temp->p=t->p;
         }
         else //isright
         {
            t->p->r=temp;
            temp->p=t->p;
         }
         delete t->r;
         t->p=0;
         t->r=0;
         t->l=0;
         delete t;
         t=temp;
         splitleft(t,v,c);
      }
      else
      {
         splitleft(t->l,v,c);
         splitleft(t->r,v,c);
      }
   }
}
//--------------------------------------------------
//split tree along variable v at cutpoint c retaining only 
//part of the tree that is ``right'' of this v,c rule
void splitright(tree::tree_p t, size_t v, size_t c)
{
   tree::tree_p temp;

   if(t->l) //not terminal node
   {
      if(t->v==v && t->c <= c)
      {
         temp=t->r;
         if(t->isleft())
         {
            t->p->l=temp;
            temp->p=t->p;
         }
         else //isright
         {
            t->p->r=temp;
            temp->p=t->p;
         }
         delete t->l;
         t->p=0;
         t->l=0;
         t->r=0;
         delete t;
         t=temp;
         splitright(t,v,c);
      }
      else
      {
         splitright(t->l,v,c);
         splitright(t->r,v,c);
      }
   }
}
//--------------------------------------------------
//does an actual merge (randomly chosen) 
bool merge(tree::tree_p tl, tree::tree_p tr, tree::tree_p t, size_t v, size_t c, rn& gen)
{
   bool m1,m2;
   tree::tree_cp temptl,temptr;
   int tnwl=0,tnwr=0;
   double u;

   u=gen.uniform();

   if(arenodesleafs(tl,tr)) {  //merging type 3
      if(u<0.5) {
         t->v=tl->v;
         t->c=tl->c;
         t->theta=tl->theta;  //doesn't matter actually it will be overwritten.
         t->l=0;
         t->r=0;
      }
      else 
      {
         temptl=tl;
         temptr=tr;
         t->v=v;
         t->c=c;
         t->l=new tree(*temptl);
         t->r=new tree(*temptr);
         t->l->p=t;
         t->r->p=t;
      }
      return true;
   }
   else if(arenodesequal(tl,tr) && !splitsonv(tl,tr,v)) {  //merging type 4
      m1=mergecount(tl->l,tr->l,v,c,&tnwl);
      m2=mergecount(tl->r,tr->r,v,c,&tnwr);
      if(u < (1.0/(tnwl+tnwr+1.0)) )
      {
         temptl=tl;
         temptr=tr;
         t->v=v;
         t->c=c;
         t->l=new tree(*temptl);
         t->r=new tree(*temptr);
         t->l->p=t;
         t->r->p=t;
      }
      else
      {
         t->v=tl->v;
         t->c=tl->c;
         t->l=new tree;
         t->r=new tree;
         t->l->p=t;
         t->r->p=t;
         tnwl=0;
         tnwr=0;
         m1=merge(tl->l,tr->l,t->l,v,c,gen);
         m2=merge(tl->r,tr->r,t->r,v,c,gen);
      }
      return (m1 & m2);
   }
   else if(splitsonv(tl,tr,v)) {  //merging type 7
      m1=mergecount(tl->r,tr,v,c,&tnwr);
      m2=mergecount(tl,tr->l,v,c,&tnwl);
      if(u < (1.0/(tnwr+tnwl+1.0)) )
      {
         temptl=tl;
         temptr=tr;
         t->v=v;
         t->c=c;
         t->l=new tree(*temptl);
         t->r=new tree(*temptr);
         t->l->p=t;
         t->r->p=t;
      }
      else if(u < ((1.0+tnwr)/(1.0+tnwr+tnwl)) )
      {
         temptl=tl->l;
         t->v=tl->v;
         t->c=tl->c;
         t->l=new tree(*temptl);
         t->l->p=t;
         t->r=new tree;
         t->r->p=t;
         m2=merge(tl->r,tr,t->r,v,c,gen);
      }
      else
      {
         temptr=tr->r;
         t->v=tr->v;
         t->c=tr->c;
         t->r=new tree(*temptr);
         t->r->p=t;
         t->l=new tree;
         t->l->p=t;
         m1=merge(tl,tr->l,t->l,v,c,gen);
      }
      if(!m1) COUT << "doh7a" << endl;
      if(!m2) COUT << "doh7b" << endl; 
      return (m1 & m2);
   }
   else if(splitsonv(tl,v) && isleaf(tr)) //merging type 1
   {
      m1=mergecount(tl->r,tr,v,c,&tnwr);
      if(u < (1.0/(tnwr + 1.0)) )
      {
         temptl=tl;
         temptr=tr;
         t->v=v;
         t->c=c;
         t->l=new tree(*temptl);
         t->r=new tree(*temptr);
         t->l->p=t;
         t->r->p=t;
      }
      else
      {
         temptl=tl->l;
         t->v=tl->v;
         t->c=tl->c;
         t->l=new tree(*temptl);
         t->l->p=t;
         t->r=new tree;
         t->r->p=t;
         m1=merge(tl->r,tr,t->r,v,c,gen);
      }
      if(!m1) COUT << "doh1(m1)" << endl;
      return m1;
   }
   else if(splitsonv(tr,v) && isleaf(tl)) //merging type 2
   {
      m2=mergecount(tl,tr->l,v,c,&tnwl);
      if(u < (1.0/(tnwl+1.0)) )
      {
         temptl=tl;
         temptr=tr;
         t->v=v;
         t->c=c;
         t->l=new tree(*temptl);
         t->r=new tree(*temptr);
         t->l->p=t;
         t->r->p=t;        
      }
      else
      {
         temptr=tr->r;
         t->v=tr->v;
         t->c=tr->c;
         t->r=new tree(*temptr);
         t->r->p=t;
         t->l=new tree;
         t->l->p=t;
         m2=merge(tl,tr->l,t->l,v,c,gen);
      }
      if(!m2) COUT << "doh2(m2)" << endl;
      return m2;
   }
   else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tr,v)) { //merge type 6(i)
      m1=mergecount(tl,tr->l,v,c,&tnwr);
      if(u < (1.0/(tnwr + 1.0)) )
      {
         temptl=tl;
         temptr=tr;
         t->v=v;
         t->c=c;
         t->l=new tree(*temptl);
         t->r=new tree(*temptr);
         t->l->p=t;
         t->r->p=t;        
      }
      else
      {
         temptr=tr->r;
         t->v=tr->v;
         t->c=tr->c;
         t->r=new tree(*temptr);
         t->r->p=t;
         t->l=new tree;
         t->l->p=t;
         m1=merge(tl,tr->l,t->l,v,c,gen);
      }
      if(!m1) COUT << "doh6i(m1)" << endl;
      return m1;
   }
   else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tl,v)) { //merge type 6(ii)
      m2=mergecount(tl->r,tr,v,c,&tnwl);
      if(u < (1.0/(tnwl + 1.0)) )
      {
         temptl=tl;
         temptr=tr;
         t->v=v;
         t->c=c;
         t->l=new tree(*temptl);
         t->r=new tree(*temptr);
         t->l->p=t;
         t->r->p=t;
      }
      else
      {
         temptl=tl->l;
         t->v=tl->v;
         t->c=tl->c;
         t->l=new tree(*temptl);
         t->l->p=t;
         t->r=new tree;
         t->r->p=t;
         m2=merge(tl->r,tr,t->r,v,c,gen);
      }
      if(!m2) COUT << "doh6ii(m2)" << endl;
      return m2;
   }
   else if(!splitsonv(tl,v) && isleaf(tr)) { //merge type 5(i)
      temptl=tl;
      temptr=tr;
      t->v=v;
      t->c=c;
      t->l=new tree(*temptl);
      t->r=new tree(*temptr);
      t->l->p=t;
      t->r->p=t;
      return true;
   }
   else if(!splitsonv(tr,v) && isleaf(tl)) { //merge type 5(ii)
      temptl=tl;
      temptr=tr;
      t->v=v;
      t->c=c;
      t->l=new tree(*temptl);
      t->r=new tree(*temptr);
      t->l->p=t;
      t->r->p=t;
      return true;
   }
   else // default type aka type 8
   {
      temptl=tl;
      temptr=tr;
      t->v=v;
      t->c=c;
      t->l=new tree(*temptl);
      t->r=new tree(*temptr);
      t->l->p=t;
      t->r->p=t;
      return true;
   }

   return false;
}
//--------------------------------------------------
// only to get nways, not to actually do the merge.
bool mergecount(tree::tree_p tl, tree::tree_p tr, size_t v, size_t c, int* nways)
{
   bool m1,m2;
   int tnwl=0,tnwr=0;

   if(arenodesleafs(tl,tr)) {  //merging type 3
      *nways += 2;
      return true;
   }
   else if(arenodesequal(tl,tr) && !splitsonv(tl,tr,v)) {  //merging type 4
      *nways += 1;
      m1=mergecount(tl->l,tr->l,v,c,&tnwl);
      m2=mergecount(tl->r,tr->r,v,c,&tnwr);
      *nways += (tnwr*tnwl);
      return (m1 & m2);
   }
   else if(splitsonv(tl,tr,v)) {  //merging type 7
      *nways += 1;  //for this one
      m1=mergecount(tl->r,tr,v,c,&tnwr);
      m2=mergecount(tl,tr->l,v,c,&tnwl);
      *nways+= (tnwr+tnwl);
      if(!m1) COUT << "doh7a" << endl;
      if(!m2) COUT << "doh7b" << endl; 
      return (m1 & m2);
   }
   else if(splitsonv(tl,v) && isleaf(tr)) //merging type 1
   {
      *nways += 1; //for this one
      m1=mergecount(tl->r,tr,v,c,&tnwr);
      *nways += tnwr;
      if(!m1) COUT << "doh1(m1)" << endl;
      return m1;
   }
   else if(splitsonv(tr,v) && isleaf(tl)) //merging type 2
   {
      *nways += 1; //for this one
      m2=mergecount(tl,tr->l,v,c,&tnwl);
      *nways += tnwl;
      if(!m2) COUT << "doh2(m2)" << endl;
      return m2;
   }
   else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tr,v)) { //merge type 6(i)
      *nways += 1; //for this one
      m1=mergecount(tl,tr->l,v,c,&tnwr);
      *nways += tnwr;
      if(!m1) COUT << "doh6i(m1)" << endl;
      return m1;
   }
   else if(!isleaf(tl) && !isleaf(tr) && splitsonv(tl,v)) { //merge type 6(ii)
      *nways +=1 ; //for this one
      m2=mergecount(tl->r,tr,v,c,&tnwl);
      *nways += tnwl;
      if(!m2) COUT << "doh6ii(m2)" << endl;
      return m2;
   }
   else if(!splitsonv(tl,v) && isleaf(tr)) { //merge type 5(i)
      *nways += 1; //for this one
      return true;
   }
   else if(!splitsonv(tr,v) && isleaf(tl)) { //merge type 5(ii)
      *nways += 1; //for this one
      return true;
   }
   else // default type aka type 8
   {
      *nways += 1; //for this one
      return true;
   }

   return false;
}

#endif
