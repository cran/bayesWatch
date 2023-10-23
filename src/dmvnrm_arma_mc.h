#ifndef dmvnrm_arma_mc_H
#define dmvnrm_arma_mc_H

#include <RcppArmadillo.h>
#include <numeric>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace Rcpp;
using namespace arma;

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat);
double log_dmvnrm_arma_regular(arma::mat const &data_x,  
                               arma::rowvec const &mean,  
                               arma::mat const &prec);

#endif
