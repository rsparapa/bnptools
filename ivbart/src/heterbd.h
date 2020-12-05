#ifndef GUARD_heterbd_h
#define GUARD_heterbd_h

#include "rrn.h"
#include "info.h"
#include "tree.h"
#include "treefuns.h"
#include "bartfuns.h"
#include "heterbartfuns.h"

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen);

#include <iostream>
//#include "heterbd.h"
#include "heterbartfuns.h"

//using std::cout;
using std::endl;

bool heterbd(tree& x, xinfo& xi, dinfo& di, pinfo& pi, double *sigma, rn& gen)
{
   tree::npv goodbots;  //nodes we could birth at (split on)
   double PBx = getpb(x,xi,pi,goodbots); //prob of a birth at x

   if(gen.uniform() < PBx) { //do birth or death

      //--------------------------------------------------
      //draw proposal
      tree::tree_p nx; //bottom node
      size_t v,c; //variable and cutpoint
      double pr; //part of metropolis ratio from proposal and prior
      bprop(x,xi,pi,goodbots,PBx,nx,v,c,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      size_t nr,nl; //counts in proposed bots
      double bl,br; //sums of weights
      double Ml, Mr; //weighted sum of y in proposed bots
      hetergetsuff(x,nx,v,c,xi,di,nl,bl,Ml,nr,br,Mr,sigma);

      //--------------------------------------------------
      //compute alpha
      double alpha=0.0, lalpha=0.0;
      double lhl, lhr, lht;
      if((nl>=5) && (nr>=5)) { //cludge?
         lhl = heterlh(bl,Ml,pi.tau);
         lhr = heterlh(br,Mr,pi.tau);
         lht = heterlh(bl+br,Ml+Mr,pi.tau);
   
         alpha=1.0;
         lalpha = log(pr) + (lhl+lhr-lht); 
         lalpha = std::min(0.0,lalpha);
      }

      //--------------------------------------------------
      //try metrop
      double mul,mur; //means for new bottom nodes, left and right
      double uu = gen.uniform();
      bool dostep = (alpha > 0) && (log(uu) < lalpha);
      if(dostep) {
         mul = heterdrawnodemu(bl,Ml,pi.tau,gen);
         mur = heterdrawnodemu(br,Mr,pi.tau,gen);
         x.birthp(nx,v,c,mul,mur);
         return true;
      } else {
         return false;
      }
   } else {
      //--------------------------------------------------
      //draw proposal
      double pr;  //part of metropolis ratio from proposal and prior
      tree::tree_p nx; //nog node to death at
      dprop(x,xi,pi,goodbots,PBx,nx,pr,gen);

      //--------------------------------------------------
      //compute sufficient statistics
      double br,bl; //sums of weights
      double Ml, Mr; //weighted sums of y
      hetergetsuff(x, nx->getl(), nx->getr(), xi, di, bl, Ml, br, Mr, sigma);

      //--------------------------------------------------
      //compute alpha
      double lhl, lhr, lht;
      lhl = heterlh(bl,Ml,pi.tau);
      lhr = heterlh(br,Mr,pi.tau);
      lht = heterlh(bl+br,Ml+Mr,pi.tau);

      double lalpha = log(pr) + (lht - lhl - lhr);
      lalpha = std::min(0.0,lalpha);

      //--------------------------------------------------
      //try metrop
      //double a,b,s2,yb;
      double mu;
      if(log(gen.uniform()) < lalpha) {
         mu = heterdrawnodemu(bl+br,Ml+Mr,pi.tau,gen);
         x.deathp(nx,mu);
         return true;
      } else {
         return false;
      }
   }
}

#endif
