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

#include "rn.h"
#include "RcppEigen.h"

typedef std::vector<double> v1d;
typedef std::vector<v1d> v2d;
typedef std::vector<v2d> v3d;

void draw_InverseWishart(Eigen::MatrixXd& iwish, double df, Eigen::MatrixXd& Psi, arn gen);

void draw_varianceComponents(Eigen::MatrixXd& varcov, Eigen::MatrixXd& sum_rereT, v1d scales, double df, size_t L, size_t prior_no, arn gen);

void get_reParams(Eigen::MatrixXd& RE_varcov, Eigen::VectorXd& RE_mean, Eigen::MatrixXd vc, Eigen::MatrixXd zzti, Eigen::MatrixXd zdiffi, double sigma);

void draw_randomEffects(Eigen::VectorXd& randomEffect, Eigen::VectorXd RE_mean, Eigen::MatrixXd RE_varcov, arn gen);
