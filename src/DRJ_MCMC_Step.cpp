#include <omp.h>
#include <math.h>
#include <RcppArmadillo.h>
#include <boost/math/special_functions/erf.hpp>
//#include "omp_set_num_cores.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Completion algorithm for precision matrix estimate.  See Murph et al 2023 for explaination.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::mat complete_lambda(const arma::mat& orig_chol_mat, const arma::mat& current_G, int p, int cores){
  // Complete a precision matrices' cholesky decomposition to match graph structure G.
//  omp_set_num_cores( cores );
  arma::mat chol_mat = orig_chol_mat;
  int rplus1;
  for(int r = 0; r < p; r++){
    if(r == 0){
      for(int i = (r+1); i < p; i++){
        if(current_G(r, i) == 0){
          chol_mat(r,i) = 0;
        }
      }
    } else {
      rplus1 = r + 1;
//     #pragma omp parallel
//     {
//     #pragma omp for
        for(int s = rplus1; s < p; s++){
          double chol_rs  = 0;
          if(current_G(r, s) == 0){
            for(int i = 0; i < r; i++){
              chol_rs    += chol_mat(i,r) * chol_mat(i,s);
            }
            chol_rs      *= -(1/chol_mat(r,r));
            chol_mat(r,s) = chol_rs;
          }
        }
//     }
    }
  }
  
  return chol_mat;
}

//' Single step in Lenkowski's Double Reversible Metropolis-Hastings algorithm.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List DRJ_MCMC_singlestep(const arma::mat& current_lambda, const arma::mat& lambda_0, const arma::mat& current_G, 
                   const int p, const int cores, const int edge_updated_i, const int edge_updated_j, const arma::mat& scale_matrix, 
                   const int n_regime, const arma::vec mean_vector_regime, const arma::mat nS2,
                   const int b, const double spread_parameter_sd2, const arma::vec& mean_hyperparameter,
                   const double lambda_hyperparameter, const double g_prior) {
  // I decided to move the lambda_0 calculation to the R level (using Rcpp) because the FORTRAN used
  // in that method conflicts with armadillo.
  
  arma::mat lambda_new_tilde;
  arma::mat new_G                       = current_G;
  new_G(edge_updated_i, edge_updated_j) = 1 - current_G(edge_updated_i, edge_updated_j);
  arma::mat lambda_0_cholesky           = arma::chol(lambda_0, "upper");
  double log_MH_ratio;
  
  List return_items;
  
  if(current_G.at(edge_updated_i, edge_updated_j)){
    // | -------------------- REMOVE -------------------------- |
    // In this case, our new G will REMOVE an edge.
    arma::mat lambda_cholesky = arma::chol(current_lambda, "upper");
    double nu                 = 0;
    for(int r = 0; r < edge_updated_i; r++){
      nu += lambda_cholesky.at(r, edge_updated_i) * lambda_cholesky.at(r, edge_updated_j);
    }
    nu                       *= -(1/lambda_cholesky.at(edge_updated_i, edge_updated_i));
    
    double nu_tilde           = lambda_0_cholesky.at(edge_updated_i, edge_updated_j);
    double gamma              = lambda_cholesky.at(edge_updated_i, edge_updated_j) - nu;
    double gamma_tilde        = nu_tilde + arma::randn<double>() * sqrt(spread_parameter_sd2);
    
    // Move Lambda_cholesky and Lambda_0 cholesky opposite graph structures:
    arma::mat lambda_cholesky_tilde_new                       = lambda_cholesky;
    lambda_cholesky_tilde_new(edge_updated_i, edge_updated_j) = gamma;
    arma::mat lambda_tilde_new_cholesky                       = complete_lambda(lambda_cholesky_tilde_new, 
                                                                                 new_G, p,  cores);
    
    arma::mat lambda_0_cholesky_new                       = lambda_0_cholesky;
    lambda_0_cholesky_new(edge_updated_i, edge_updated_j) = gamma_tilde;
    arma::mat lambda_0_new_cholesky                       = complete_lambda(lambda_0_cholesky_new, 
                                                                             current_G, p, cores);
    lambda_new_tilde           = lambda_tilde_new_cholesky.t() * lambda_tilde_new_cholesky;
    arma::mat lambda_0_new     = lambda_0_new_cholesky.t() * lambda_0_new_cholesky;
    
    // Calculate the double MH ratio
    arma::mat diff_lambdas       = lambda_new_tilde - current_lambda;
    arma::mat scale_mat_shifted  = scale_matrix + nS2;
    arma::mat diff_lambdas_0     = lambda_0_cholesky - lambda_0_new;
    double log_phi_ratio         = log(lambda_cholesky.at(edge_updated_i, edge_updated_i)) - 
                                    log(lambda_0_new_cholesky.at(edge_updated_i, edge_updated_i));
    // I think that this is the only thing that changes code-wise in the MH ratio.
    double log_normal_transition = -(-pow((gamma - nu), 2.0) + pow((gamma_tilde - nu_tilde), 2))/(2*spread_parameter_sd2);
    //
    double log_trace_terms       = -0.5*arma::trace(diff_lambdas.t()*scale_mat_shifted) + 0.5*arma::trace(diff_lambdas_0.t()*scale_matrix);
    
    log_MH_ratio                 = log_trace_terms + log_normal_transition + log_phi_ratio + log(1.0-g_prior) - log(g_prior);
    
  } else {
    // | -------------------- ADD -------------------------- |
    // In this case, our new G will ADD an edge.
    arma::mat lambda_cholesky = arma::chol(current_lambda, "upper");
    double nu_tilde           = 0;
    for(int r = 0; r < edge_updated_i; r++){
      nu_tilde += lambda_0_cholesky.at(r, edge_updated_i) * lambda_0_cholesky.at(r, edge_updated_j);
    }
    nu_tilde                 *= -(1/lambda_0_cholesky.at(edge_updated_i, edge_updated_i));
    
    double nu                 = lambda_cholesky.at(edge_updated_i, edge_updated_j);
    double gamma              = nu + arma::randn<double>() * sqrt(spread_parameter_sd2);
    double gamma_tilde        = lambda_0_cholesky.at(edge_updated_i, edge_updated_j) - nu_tilde;
    
    // Move Lambda_cholesky and Lambda_0 cholesky opposite graph structures:
    arma::mat lambda_cholesky_tilde_new                       = lambda_cholesky;
    lambda_cholesky_tilde_new(edge_updated_i, edge_updated_j) = gamma;
    arma::mat lambda_tilde_new_cholesky                       = complete_lambda(lambda_cholesky_tilde_new, 
                                                                                 new_G, p, cores);
    
    arma::mat lambda_0_cholesky_new                       = lambda_0_cholesky;
    lambda_0_cholesky_new(edge_updated_i, edge_updated_j) = gamma_tilde;
    arma::mat lambda_0_new_cholesky                       = complete_lambda(lambda_0_cholesky_new, 
                                                                             current_G, p, cores);
    lambda_new_tilde           = lambda_tilde_new_cholesky.t() * lambda_tilde_new_cholesky;
    arma::mat lambda_0_new     = lambda_0_new_cholesky.t() * lambda_0_new_cholesky;
    
    // Calculate the double MH ratio
    arma::mat diff_lambdas       = lambda_new_tilde - current_lambda;
    arma::mat scale_mat_shifted  = scale_matrix + nS2;
    arma::mat diff_lambdas_0     = lambda_0_cholesky - lambda_0_new;
    double log_phi_ratio         = log(lambda_cholesky.at(edge_updated_i, edge_updated_i)) - 
                                     log(lambda_0_new_cholesky.at(edge_updated_i, edge_updated_i));
    // I think that this is the only thing that changes code-wise in the MH ratio.
    double log_normal_transition = -(pow((gamma - nu), 2) - pow((gamma_tilde - nu_tilde), 2))/(2*spread_parameter_sd2);
    //
    double log_trace_terms       = -0.5*arma::trace(diff_lambdas.t()*scale_mat_shifted) + 0.5*arma::trace(diff_lambdas_0.t()*scale_matrix);
    
    // Note: I changed the prior distribution on the graph structure G.  The 'g.prior' fed into this algorithm now accounts for if it is an add/subtract prob.
    //       I should not longer do it at this level.
    log_MH_ratio                 = log_trace_terms + log_normal_transition + log_phi_ratio + log(1.0-g_prior) - log(g_prior); //+ log(g_prior) - log(1.0-g_prior);
  }
  
  double new_lambda_hyperparameter  = lambda_hyperparameter + n_regime;
  arma::vec new_mean_hyperparameter = (1/(n_regime + lambda_hyperparameter)) * 
                                        (n_regime*mean_vector_regime + mean_hyperparameter*lambda_hyperparameter);
  
  arma::mat new_covariance_rescaled = arma::inv(new_lambda_hyperparameter * lambda_new_tilde);
  arma::mat covariance_cholesky     = arma::chol(new_covariance_rescaled, "lower");
  arma::vec std_norm_vals           = arma::randn<arma::vec>(p);
  arma::vec new_mu                  = new_mean_hyperparameter + covariance_cholesky * std_norm_vals;
  
  return_items["log_MH_ratio"]      = log_MH_ratio;
  return_items["new_lambda"]        = lambda_new_tilde;
  return_items["new_mu"]            = new_mu;

return return_items;
}


//' Unnormalized kernel log value for a NW posterior distribution, given the data.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
double log_dNormalWishart_posterior_unnormalized( const arma::mat& data_matrix, const arma::vec& m_hyperparameter,
                                               const arma::mat& scale_matrix, const double lambda_hyperparameter,
                                               const double nu_wishartDF,
                                               const arma::vec& observed_mu, const arma::mat& observed_precision) {
  int n = data_matrix.n_rows;
  int p = data_matrix.n_cols;
  arma::vec col_means(p);
  arma::mat nS2;
  arma::mat centered_data;
  arma::vec mean_minus_m;
  double log_likelihood, log_det_term, sign;

  centered_data.zeros(n,p);
  for(int i = 0; i < p; i++){
    col_means(i)                  = arma::mean(data_matrix.col(i));
    centered_data.col(i)          = data_matrix.col(i) - col_means(i);
  }
  nS2                             = (centered_data.t()) * centered_data;
  mean_minus_m                    = col_means - m_hyperparameter;
  
  arma::vec posterior_mu          = (lambda_hyperparameter * m_hyperparameter + n * col_means) / (n + lambda_hyperparameter);
  arma::mat posterior_inv_scale   = arma::inv(scale_matrix) + nS2 + 
                                     ((n*lambda_hyperparameter) / (n + lambda_hyperparameter)) * (mean_minus_m * (mean_minus_m.t()));
  double posterior_wishartDF      = nu_wishartDF + n;
  arma::vec posterior_mu_centered = observed_mu - posterior_mu;
  arma::log_det(log_det_term, sign, observed_precision);
  
  log_likelihood  = 0.5* (posterior_wishartDF-2) * (log_det_term);
  log_likelihood += -0.5 * arma::trace(posterior_inv_scale * observed_precision);
  log_likelihood += 0.5 * (log_det_term);
  log_likelihood += -0.5 * arma::as_scalar((posterior_mu_centered.t()) * observed_precision * posterior_mu_centered);
  
  return(log_likelihood);
}

//' Unnormalized kernel log value for a NW distribution.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
double log_dNormalWishart_unnormalized( const arma::vec& m_hyperparameter, const arma::mat& posterior_inv_scale, 
                                        const double lambda_hyperparameter, const double nu_wishartDF,
                                        const arma::vec& observed_mu, const arma::mat& observed_precision) {
  double log_likelihood, log_det_term, sign;
  
  // arma::mat posterior_inv_scale   = arma::inv(scale_matrix);
  double posterior_wishartDF      = nu_wishartDF;
  arma::vec posterior_mu_centered = observed_mu - m_hyperparameter;
  int p                           = observed_precision.n_cols;
  arma::log_det(log_det_term, sign, observed_precision);
 
  log_likelihood  = 0.5 * (posterior_wishartDF-2) * (log_det_term);
  log_likelihood += -0.5 * p * log(2*M_PI);
  log_likelihood += -0.5 * arma::trace(posterior_inv_scale * observed_precision);
  log_likelihood += p * 0.5 * log(lambda_hyperparameter) + 0.5 * (log_det_term);
  log_likelihood += -0.5 * lambda_hyperparameter * arma::as_scalar((posterior_mu_centered.t()) * observed_precision * posterior_mu_centered);
  
  return(log_likelihood);
}


//' Sample a new mu according to a NW distribution.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::vec rmu_0( const arma::mat& sigma_0, const arma::mat& sum_precision_matrices,
                 const arma::vec& sum_precision_times_mu, const arma::vec& m_hyperparameter) {
  
  int p                          = sigma_0.n_cols;
  arma::mat posterior_precision  = sigma_0 + sum_precision_matrices;
  arma::mat posterior_covariance = arma::inv(posterior_precision);
  arma::mat chol_covariance      = arma::chol(posterior_covariance, "lower");
    
  arma::vec posterior_mu         = posterior_covariance * (sum_precision_times_mu + m_hyperparameter);
  arma::vec std_norm_vals        = arma::randn<arma::vec>(p);
  arma::vec new_mu               = posterior_mu + chol_covariance * std_norm_vals;
  return(new_mu);
}



