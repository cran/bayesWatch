#include <math.h>
#include <cmath>
#include <RcppArmadillo.h>
#include <numeric>
#include <limits>
#include "dmvnrm_arma_mc.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Calculates the probability of a component for the Gibbs sweep update.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
double calc_logprob_Gibbs_comp(const arma::mat& current_precision, const arma::vec& current_mu, const arma::vec& regime_comp_log_probs,
                             const arma::mat& current_data, int proposed_component) {
  double log_prob = 0;
  if(isnan(regime_comp_log_probs.at(proposed_component)) ){
    log_prob = -std::numeric_limits<double>::infinity();
  } else {
    log_prob = (double)arma::as_scalar(regime_comp_log_probs.at(proposed_component)) + 
                  log_dmvnrm_arma_regular(current_data,current_mu.t(),current_precision);
  }
  return log_prob;
}

//' Calculates relative probability of two possible components according to Dirichlet process.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
List gibbs_swap_btwn_two(const arma::mat& first_precision, const arma::mat& second_precision, 
                              const arma::vec& first_mu, const arma::vec& second_mu, 
                              const arma::vec& component_log_probs, const arma::vec& indices_of_split_component, 
                              const arma::mat& data_points_of_state, arma::vec& assignments_launch,
                              int first_component, int second_component, int num_gibbs_sweeps) {

  arma::vec data_value;
  int data_value_index;
  double log_prob_first, log_prob_second, prob_first, prob_second, log_max_log_prob, unif_draw;
  List return_items;
  
  double total_log_prob = 0;
  int num_elements = indices_of_split_component.n_elem;
  int length_launch = assignments_launch.n_elem;
  int num_equal_to_first = 0;
  int num_equal_to_second = 0;
    
  for(int swap_index = 0; swap_index < num_gibbs_sweeps; swap_index++) {
    for(int indices_index = 0; indices_index < num_elements; indices_index++){
      data_value_index        = (int)arma::as_scalar(indices_of_split_component.at(indices_index));
      data_value              = data_points_of_state.row(data_value_index-1).t();
      for(int count_index = 0; count_index < length_launch; count_index++){
        if( (int)arma::as_scalar(assignments_launch.at(count_index)) == first_component){
          num_equal_to_first += 1;
        }
        if( (int)arma::as_scalar(assignments_launch.at(count_index)) == second_component){
          num_equal_to_second += 1;
        }
      }
      
      if( ( num_equal_to_first>1 )&( num_equal_to_second>1 ) ){
        log_prob_first        = calc_logprob_Gibbs_comp(first_precision, first_mu, component_log_probs, data_value.t(), first_component);
        log_prob_second       = calc_logprob_Gibbs_comp(second_precision, second_mu, component_log_probs, data_value.t(), second_component);
        
        log_max_log_prob = std::max(log_prob_first,log_prob_second);
        log_prob_first   = log_prob_first - log_max_log_prob;
        log_prob_second  = log_prob_second - log_max_log_prob;
        prob_first       = exp(log_prob_first)/(exp(log_prob_first)+exp(log_prob_second));
        prob_second      = exp(log_prob_second)/(exp(log_prob_first)+exp(log_prob_second));
        unif_draw        = arma::randu<double>();
        if(isnan(prob_first)){
          assignments_launch.at(data_value_index-1) = first_component;
          continue;
          
        }
        if(unif_draw<=prob_first){
          assignments_launch.at(data_value_index-1) = first_component;
          total_log_prob += log(prob_first);
        } else {
          assignments_launch.at(data_value_index-1) = second_component;
          total_log_prob += log(prob_second);
        }
      }
    }
  }
  return_items["assignments_launch"] = assignments_launch;
  return_items["total_log_prob"] = total_log_prob;
  return return_items;
}

//' Performs swap of components in the component reassignment.
//' 
//' @noRd
//' 
// [[Rcpp::export]]
arma::vec gibbs_swap_comps(const arma::mat& data_points_of_state, arma::vec& cluster_assignments, const arma::vec& regime_comp_log_probs,
                           List precisions, List mus, int assignments_maximum, int gibbs_sweeps){
  int index_of_split = 0;
  arma::vec log_prob_of_each_comp;
  arma::vec exp_distribution;
  arma::vec data_value;
  double total_density_values;
  double max_log_prob;
  arma::mat current_precision;
  arma::vec current_mu;
  
  int n_full                = data_points_of_state.n_rows;
  double sum_density_values = 0;
  double rand_value         = 0;
  
  
  for(int gibbs_index = 0; gibbs_index < gibbs_sweeps; gibbs_index++){
//    #pragma omp parallel
//    {
//    #pragma omp for
      for(int data_value_index = 0; data_value_index< n_full; data_value_index++){
        data_value                   = data_points_of_state.row(data_value_index).t();
        log_prob_of_each_comp        = arma::zeros(assignments_maximum);
        
        for(int component_index = 0; component_index < assignments_maximum; component_index++){
          current_precision = as<arma::mat>(precisions[component_index]);
          current_mu        = as<arma::vec>(mus[component_index]);
          log_prob_of_each_comp.at(component_index) = calc_logprob_Gibbs_comp(current_precision, current_mu, 
                                                                              regime_comp_log_probs,
                                                                              data_value.t(), component_index);
        }
        max_log_prob           = (double)arma::as_scalar(arma::max(log_prob_of_each_comp));
        log_prob_of_each_comp  -= arma::ones( assignments_maximum )*max_log_prob;
        exp_distribution       = exp(log_prob_of_each_comp);
        total_density_values   = accu(exp_distribution);
  
        if(isnan(total_density_values)){
          continue;
        }
        rand_value             = arma::randu<double>() * total_density_values;
        if(isnan(rand_value)){
          continue;
        }
        sum_density_values = 0;
        for(int split_index = 0; split_index < assignments_maximum; split_index++){
          sum_density_values += exp_distribution.at(split_index);
          if(rand_value <= sum_density_values){
            index_of_split     = split_index;
            break;
          }
        }
        cluster_assignments.at(data_value_index) = index_of_split+1;
      }
//    }
  }
  return cluster_assignments;
}



