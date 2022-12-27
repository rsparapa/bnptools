/*
 * Copyright (C) 2012-2023 Matthew T. Pratola, Robert E. McCulloch,
 *                         Hugh A. Chipman and Rodney A. Sparapani
 *  
 * This file is part of nftbart.
 * getpred.h
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

typedef std::vector<tree> vtree;

void getpred(int beg, int end, size_t p, size_t m, size_t np, xinfo& xi, std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat)
{
   double *fptemp = new double[np];

   for(int i=beg;i<=end;i++) {
      for(size_t j=0;j<m;j++) {
         fit(tmat[i][j],xi,p,np,px,fptemp);
         for(size_t k=0;k<np;k++) yhat(i,k) += fptemp[k];
      }
   }

   delete [] fptemp;
}
#ifdef _OPENMP
void local_getpred(size_t nd, size_t p, size_t m, size_t np, xinfo& xi, 
std::vector<vtree>& tmat, double *px, Rcpp::NumericMatrix& yhat)
{
   int my_rank = omp_get_thread_num();
   int thread_count = omp_get_num_threads();
   int h = nd/thread_count; int beg = my_rank*h; int end = beg+h-1;
   
   getpred(beg,end,p,m,np,xi,tmat,px,yhat);
}
#endif
