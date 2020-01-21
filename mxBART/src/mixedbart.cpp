/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include "mixedbart.h"
#include "RcppEigen.h"

// Draw inverse-Wishart distribution with df and Psi
// Generally, df and mat_size are small so this method is sufficient for generating Inverse-Wishart RVs
void draw_InverseWishart(Eigen::MatrixXd& iwish, double df, Eigen::MatrixXd& Psi, arn gen){
  size_t mat_size = Psi.rows();
  Eigen::MatrixXd Psi_inv = Psi.inverse();
  Eigen::MatrixXd unitwish(mat_size,mat_size);
  Eigen::VectorXd gauss(mat_size);
  Eigen::LLT<Eigen::MatrixXd> cholesky(Psi_inv);
  unitwish.setZero();
  for(size_t d=0;d<df;d++){
    for(size_t i=0;i<mat_size;i++) {
      gauss(i)=gen.normal();
    }
  unitwish+=gauss*gauss.transpose();
  }
  iwish=(cholesky.matrixL()*unitwish*cholesky.matrixU()).inverse();
}

// Draw from the posterior of the variance-covariance matrix of the random effects
void draw_varianceComponents(Eigen::MatrixXd& varcov, Eigen::MatrixXd& sum_rereT, v1d scales, double df, size_t L, size_t prior_no, arn gen){
  const size_t mat_size = varcov.rows();
  Eigen::MatrixXd varcov_inv = varcov.inverse();
  Eigen::MatrixXd psi(mat_size,mat_size);
  psi.setZero();
  // Hierarchical inverse-Wishart (or half-T if mat_size==1)
  if(prior_no==1){
    psi=sum_rereT;
    for(size_t i=0;i<mat_size;i++) psi(i,i)=psi(i,i)+2*df*gen.gamma(.5*(df+(double)mat_size),df*varcov_inv(i,i)+pow(scales[i],-2.));
    draw_InverseWishart(varcov,L+df+mat_size-1,psi,gen);
  }
  // inverse-Wishart (or inverse-gamma if mat_size==1)
  else{
    psi=sum_rereT;
    for(size_t i=0;i<mat_size;i++) psi(i,i)=psi(i,i)+2*df*scales[i];
    draw_InverseWishart(varcov,L+df+mat_size,psi,gen);
  }
}

// Get individual random effect parameter values varcov and mean
void get_reParams(Eigen::MatrixXd& RE_varcov, Eigen::VectorXd& RE_mean, Eigen::MatrixXd vc, Eigen::MatrixXd zzti, Eigen::MatrixXd zdiffi, double sigma){
  RE_varcov=(zzti/pow(sigma,2.)+vc.inverse()).inverse();
  RE_mean=RE_varcov/pow(sigma,2.)*zdiffi;
}

// Draw the random effects
void draw_randomEffects(Eigen::VectorXd& randomEffect, Eigen::VectorXd RE_mean, Eigen::MatrixXd RE_varcov, arn gen){
  size_t mat_size=randomEffect.size();
  Eigen::LLT<Eigen::MatrixXd> cholesky(RE_varcov);
  Eigen::VectorXd gauss(mat_size);
  for(size_t i=0;i<mat_size;i++) gauss(i)=gen.normal();
  randomEffect=cholesky.matrixL()*gauss+RE_mean;
}

