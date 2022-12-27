/*
 * Copyright (C) 2012-2023 Matthew T. Pratola, Robert E. McCulloch,
 *                         Hugh A. Chipman and Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * cpsambrt_predict2.h
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
 * Rodney A. Sparapani: rsparapa@mcw.edu
 *
 */

// Draw predictive realizations at the prediciton points, xp.
// Requires the original data locations, x, the number of
// mean and sd trees, m and mh, the number of draws
// saved from the MCMC, nd, the number of cutpoints used
// in the MCMC, numcut, the number of threads to run the
// predicitons in parallel, tc and the R output of the fitted
// MCMC, fit.
RcppExport SEXP cpsambrt_predict2(
   SEXP _ixf, 
   SEXP _ixs, 
   SEXP _im,  //number of trees in mean model
   SEXP _imh, //number of trees in variance model
   SEXP _ind, //number of draws saved from the posterior
   SEXP _ixif,
   SEXP _ixis, 
   SEXP _itc, //number of parallel compute threads
   SEXP _ifit //saved fitted model object returned from cpsambrt
)
{
   //--------------------------------------------------
   //process args
   //xp prediction points
   Rcpp::NumericMatrix xfm(_ixf), xsm(_ixs);
   size_t np = xfm.ncol();
   size_t pf = xfm.nrow(), ps = xsm.nrow();
   double *xf = nullptr, *xs = nullptr;
   if(np)  {
     xf = &xfm[0];
     xs = &xsm[0];
   }

   //x training points
   //Rcpp::NumericMatrix xm(_ix);

   //make xinfo
   xinfo xif, xis;
   Rcpp::List ixif(_ixif), ixis(_ixis);
   xif.resize(pf);
   for(size_t i=0;i<pf;i++) xif[i]=Rcpp::as< std::vector<double> >(ixif[i]);
   xis.resize(ps);
   for(size_t i=0;i<ps;i++) xis[i]=Rcpp::as< std::vector<double> >(ixis[i]);

   //number of trees
   size_t m = Rcpp::as<int>(_im);
   size_t mh = Rcpp::as<int>(_imh);

   //nd
   size_t nd = Rcpp::as<int>(_ind);

   //thread count
   int tc = Rcpp::as<int>(_itc);

   Rcpp::NumericMatrix tedraw(nd,np);
   double *fp = new double[np], *sp = new double[np];
   dinfo dip, dis;
   dip.x = xf; dip.y=fp; dip.p = pf; dip.n=np; dip.tc=tc;
   dis.x = xs; dis.y=sp; dis.p = ps; dis.n=np; dis.tc=tc;

   // set up ambrt object
   ambrt ambm(m);
   ambm.settc(tc);  //set the number of threads when using OpenMP, etc.
   ambm.setxi(&xif); //set the cutpoints for this model object

   //cpsambrt fitted model object (contains posterior MCMC samples)
   Rcpp::List fit(_ifit);
   Rcpp::XPtr< std::vector<int> > e_ots(Rcpp::as<SEXP>(fit["ots"]));
   Rcpp::XPtr< std::vector<int> > e_oid(Rcpp::as<SEXP>(fit["oid"]));
   Rcpp::XPtr< std::vector<int> > e_ovar(Rcpp::as<SEXP>(fit["ovar"]));
   Rcpp::XPtr< std::vector<int> > e_oc(Rcpp::as<SEXP>(fit["oc"]));
   Rcpp::XPtr< std::vector<double> > e_otheta(Rcpp::as<SEXP>(fit["otheta"]));

   // Temporary vectors used for loading one model realization at a time.
   std::vector<int> onn(m,1);
   std::vector<std::vector<int> > oid(m, std::vector<int>(1));
   std::vector<std::vector<int> > ov(m, std::vector<int>(1));
   std::vector<std::vector<int> > oc(m, std::vector<int>(1));
   std::vector<std::vector<double> > otheta(m, std::vector<double>(1));

   // Draw realizations of the posterior predictive.
   size_t curdx=0;
   size_t cumdx=0;
   // Mean trees first
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<m;j++) {
         onn[j]=e_ots->at(i*m+j);
         oid[j].resize(onn[j]);
         ov[j].resize(onn[j]);
         oc[j].resize(onn[j]);
         otheta[j].resize(onn[j]);
         for(int k=0;k<onn[j];k++) {
            oid[j][k]=e_oid->at(cumdx+curdx+k);
            ov[j][k]=e_ovar->at(cumdx+curdx+k);
            oc[j][k]=e_oc->at(cumdx+curdx+k);
            otheta[j][k]=e_otheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)onn[j];
      }
      cumdx+=curdx;

      ambm.loadtree(0,m,onn,oid,ov,oc,otheta);
      // draw realization
      ambm.predict(&dip);
      for(size_t j=0;j<np;j++) tedraw(i,j) = fp[j];
   }

   // Save the draws and return to R.
   Rcpp::List ret;
   ret["f.test."]=tedraw;

if(mh>0) {
   Rcpp::NumericMatrix tedrawh(nd,np);

   //setup psbrt object
   psbrt psbm(mh);
   psbm.settc(tc);  //set the number of threads when using OpenMP, etc.
   psbm.setxi(&xis); //set the cutpoints for this model object

   Rcpp::XPtr< std::vector<int> > e_sts(Rcpp::as<SEXP>(fit["sts"]));
   Rcpp::XPtr< std::vector<int> > e_sid(Rcpp::as<SEXP>(fit["sid"]));
   Rcpp::XPtr< std::vector<int> > e_svar(Rcpp::as<SEXP>(fit["svar"]));
   Rcpp::XPtr< std::vector<int> > e_sc(Rcpp::as<SEXP>(fit["sc"]));
   Rcpp::XPtr< std::vector<double> > e_stheta(Rcpp::as<SEXP>(fit["stheta"]));

   std::vector<int> snn(mh,1);
   std::vector<std::vector<int> > sid(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sv(mh, std::vector<int>(1));
   std::vector<std::vector<int> > sc(mh, std::vector<int>(1));
   std::vector<std::vector<double> > stheta(mh, std::vector<double>(1));

   // Variance trees second
   cumdx=0;
   curdx=0;
   for(size_t i=0;i<nd;i++) {
      curdx=0;
      for(size_t j=0;j<mh;j++) {
         snn[j]=e_sts->at(i*mh+j);
         sid[j].resize(snn[j]);
         sv[j].resize(snn[j]);
         sc[j].resize(snn[j]);
         stheta[j].resize(snn[j]);
         for(int k=0;k<snn[j];k++) {
            sid[j][k]=e_sid->at(cumdx+curdx+k);
            sv[j][k]=e_svar->at(cumdx+curdx+k);
            sc[j][k]=e_sc->at(cumdx+curdx+k);
            stheta[j][k]=e_stheta->at(cumdx+curdx+k);
         }
         curdx+=(size_t)snn[j];
      }
      cumdx+=curdx;

      psbm.loadtree(0,mh,snn,sid,sv,sc,stheta);
      // draw realization
      psbm.predict(&dis);
      for(size_t j=0;j<np;j++) tedrawh(i,j) = sp[j];
      ret["s.test."]=tedrawh;
   }
}
   if(fp) delete [] fp;
   if(sp) delete [] sp;

   return ret;
}



