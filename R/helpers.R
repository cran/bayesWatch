# Helper function for bd_gc_mcmc main method.
require(parallel)


#' Perform a Gibbs sweep on the current regime vector during the MCMC sampling.
#'
#' @param my_states vector. Current regime vector.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param transition_probabilities vector. Current estimate of unknown transition probabilities for Markov chain.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param n.cores integer. Number of cores available for this calculation.
#' @param min_regime_length integer. Gibbs sweep will not create a regime smaller than this.
#' @param verbose  logical. Will print to console if TRUE.
#'
#'
#' @noRd
#'
#'
update_states_gibbs       = function(my_states,
                                     data_points_Z,
                                     Z_timepoint_indices,
                                     previous_model_fits,
                                     transition_probabilities,
                                     hyperparameters,
                                     n.cores,
                                     min_regime_length = 1,
                                       verbose = FALSE) {
  # We are going to perform sequential Gibbs updates with two restrictions:
  ## 1) This process cannot create a singleton state (UPDATE: cannot move below min_regime_length);
  ## 2) This process will not require the drawing of a new model.
  
  # The idea here is to create a very fast update where states at which a change occurs
  # can swap to the next highest/next lowest state.
  
  # First, reject if there is only one state represented (the only thing that can handle
  # this is a split from the other update_states method):
  
  accepted                 = 0
  linger_parameter         = hyperparameters$alpha
  move_parameter           = hyperparameters$beta
  continous                = as.numeric(!hyperparameters$not.cont)
  coin_flip                = rbinom(1, 1, 0.5)
  
  if (coin_flip) {
    lower_value              = 2
    upper_value              = length(my_states) - 1
  } else {
    lower_value              = length(my_states) - 1
    upper_value              = 2
  }
  
  for (index in lower_value:upper_value) {
    # Never allow for a change in the first or the last states.  This could cause
    # a state to be removed/added.
    
    if (my_states[index - 1] != my_states[index + 1]) {
      # There is a state change, so we check to see if a state change is allowed.
      if (((my_states[index] != my_states[index + 1]) |
           (my_states[index] != my_states[index - 1]))) {
        # Now, whatever change is made will leave at least 2 observations in every state.
        first_index       = Z_timepoint_indices[[index]]$timepoint_first_index
        last_index        = Z_timepoint_indices[[index]]$timepoint_last_index
        temp_data         = data_points_Z[first_index:last_index,]
        n_value           = nrow(temp_data)
        
        upper_bound_is_equal_temp = hyperparameters$upper_bound_is_equal[first_index:last_index,]
        lower_bound_is_equal_temp = hyperparameters$lower_bound_is_equal[first_index:last_index,]
        is_missing_temp           = hyperparameters$is_missing[first_index:last_index,]
        
        if (my_states[index] == my_states[index - 1]) {
          new_state             = my_states[index + 1]
          my_states_temp        = my_states
          my_states_temp[index] = new_state
          
          if ((sum(my_states_temp == my_states[index - 1]) < min_regime_length) |
              (sum(my_states_temp == my_states[index + 1]) < min_regime_length) |
              (sum(my_states_temp == my_states[index]) < min_regime_length)) {
            if (  verbose)
              cat(
                paste(
                  "-----> Gibbs swap attempted at timepoint",
                  index,
                  "changing from regime",
                  my_states[index],
                  "to regime",
                  new_state,
                  "but would have created a regime below the min length.\n"
                )
              )
            next
          }
          
          
          # Get the latent data fit to the state ABOVE:
          temp_my_states           = my_states
          temp_my_states[index]    = temp_my_states[index + 1]
          previous_model_fits_temp = update_regime_components(
            temp_my_states[index + 1],
            my_states[index],
            index,
            previous_model_fits,
            Z_timepoint_indices,
            hyperparameters
          )
          previous_model_fits_temp = previous_model_fits_temp$previous_model_fits_item
          
          first_index         = Z_timepoint_indices[[index]]$timepoint_first_index
          last_index          = Z_timepoint_indices[[index]]$timepoint_last_index
          temp_data_above     = data_points_Z[first_index:last_index,]
          upper_bound_at_obs  = hyperparameters$upper_bound_is_equal[first_index:last_index,]
          lower_bound_at_obs  = hyperparameters$lower_bound_is_equal[first_index:last_index,]
          is_missing_at_obs   = hyperparameters$is_missing[first_index:last_index,]
          
          data_new_state      = temp_data_above
          data_old_state      = temp_data_above
          log_prob_new_model  = get_mix_log_dens_at_obs(
            index,
            temp_my_states,
            previous_model_fits_temp,
            data_new_state,
            Z_timepoint_indices
          )
          log_prob_old_model  = get_mix_log_dens_at_obs(
            index,
            my_states,
            previous_model_fits,
            data_old_state,
            Z_timepoint_indices
          )
          
          log_prob_forward    = log(1 - transition_probabilities[my_states[index -
                                                                             1]])
          log_prob_backward   = log(transition_probabilities[my_states[index -
                                                                         1]])
          
        } else {
          new_state             = my_states[index - 1]
          my_states_temp        = my_states
          my_states_temp[index] = new_state
          
          if ((sum(my_states_temp == my_states[index - 1]) < min_regime_length) |
              (sum(my_states_temp == my_states[index + 1]) < min_regime_length) |
              (sum(my_states_temp == my_states[index]) < min_regime_length)) {
            if (  verbose)
              cat(
                paste(
                  "-----> Gibbs swap attempted at timepoint",
                  index,
                  "changing from regime",
                  my_states[index],
                  "to regime",
                  new_state,
                  "but would have created a regime below the min length.\n"
                )
              )
            next
          }
          
          # Get the latent data fit to the state BELOW:
          temp_my_states           = my_states
          temp_my_states[index]    = temp_my_states[index - 1]
          previous_model_fits_temp = update_regime_components(
            temp_my_states[index - 1],
            my_states[index],
            index,
            previous_model_fits,
            Z_timepoint_indices,
            hyperparameters
          )
          previous_model_fits_temp = previous_model_fits_temp$previous_model_fits_item
          
          
          first_index              = Z_timepoint_indices[[index]]$timepoint_first_index
          last_index               = Z_timepoint_indices[[index]]$timepoint_last_index
          temp_data_below          = data_points_Z[first_index:last_index,]
          upper_bound_at_obs       = hyperparameters$upper_bound_is_equal[first_index:last_index,]
          lower_bound_at_obs       = hyperparameters$lower_bound_is_equal[first_index:last_index,]
          is_missing_at_obs        = hyperparameters$is_missing[first_index:last_index,]
          
          data_old_state      = temp_data_below
          data_new_state      = temp_data_below
          log_prob_new_model  = get_mix_log_dens_at_obs(
            index,
            temp_my_states,
            previous_model_fits_temp,
            data_new_state,
            Z_timepoint_indices
          )
          log_prob_old_model  = get_mix_log_dens_at_obs(
            index,
            my_states,
            previous_model_fits,
            data_old_state,
            Z_timepoint_indices
          )
          
          log_prob_forward   = log(transition_probabilities[my_states[index -
                                                                        1]])
          log_prob_backward  = log(1 - transition_probabilities[my_states[index -
                                                                            1]])
        }
        
        
        
        logprob_of_state_change = log_prob_new_model + log_prob_forward
        logprob_of_state_remain = log_prob_old_model + log_prob_backward
        normalization_term      = max(logprob_of_state_change, logprob_of_state_remain)
        prob_of_state_change    = exp(logprob_of_state_change - normalization_term)
        prob_of_state_remain    = exp(logprob_of_state_remain - normalization_term)
        prob_of_state_change    = prob_of_state_change / (prob_of_state_change +
                                                            prob_of_state_remain)
        rand_value              = runif(1)
        
        if ((length(prob_of_state_change) == 0) |
            (is.na(prob_of_state_change)) |
            (is.nan(prob_of_state_change))) {
          if (  verbose)
            cat(
              paste(
                "-----> Probability values are poor in the Gibbs sampler",
                "the value of prob_of_state_change is:",
                prob_of_state_change,
                "\n"
              )
            )
          stop("There were precision issues with the Gibbs Swap algorithm.")
        }
        if (!is.na(prob_of_state_change)) {
          if (rand_value <= prob_of_state_change) {
            old_state             = my_states[index]
            if (old_state != new_state) {
              if (  verbose)
                cat("gibbs swap changed the state!")
              previous_model_fits = update_regime_components(
                new_state,
                old_state,
                index,
                previous_model_fits,
                Z_timepoint_indices,
                hyperparameters
              )
              previous_model_fits = previous_model_fits$previous_model_fits_item
            } else {
              if (  verbose)
                cat(
                  "I'M PRETTY SURE THAT THIS SHOULD NEVER HAPPEN!!!!"
                )
            }
            my_states             = my_states_temp
            accepted              = 1
          } else {
            if (  verbose)
              cat("gibbs swap did not change the state!")
            accepted              = 1
          }
        }
        
      }
    }
    # If none of these if statements are triggered, then there is no change.
  }
  
  # Since no states were removed/added, there is no need to relabel the states.
  # There WILL be a need to refit, which we do directly after calling this method.
  return(
    list(
      previous_model_fits_item = previous_model_fits,
      my_states_item = my_states,
      accepted_item = accepted
    )
  )
}

#' Method to estimate a covariance matrix with missing.  Inputs missing values with column mean, then
#' performs covariance calculation with the full matrix.
#'
#' @param data_w_missing matrix. Data with missing values.
#'
#'
#' @noRd
#'
#'
empirical_cov_w_missing   = function(data_w_missing,
                                     Z_timepoint_indices,
                                     previous_states) {
  # In response to some reading that I've done on the issues of using the "pairwise.complete.obs," I am going
  # to write a better method to get a good empirical covariance matrix from data with missing values.
  # Essentially, I am going to cluster the data, impute each missing value with the average value of that value's cluster,
  # then take the covariance of the imputed data.
  
  # I originally wanted to do the clustering approach, but the kmeans w missing values has turned out to be more of
  # a meaty problem than I was expecting.  I may return to do that lift if the simple mean-imputation approach is not
  temp_full_data      = data_w_missing
  var_remove_func     = function(x) {
    var(x, na.rm = T)
  }
  mean_remove_func    = function(x) {
    mean(x, na.rm = T)
  }
  for (day_index in 1:length(previous_states)) {
    # For each day, impute values from a normal distribution with a mean & var that
    # match that day's mean and variance.
    min_index           = Z_timepoint_indices[[day_index]]$timepoint_first_index
    max_index           = Z_timepoint_indices[[day_index]]$timepoint_last_index
    temp_data           = data_w_missing[min_index:max_index, ]
    means_of_cols_temp  = apply(temp_data, 2, mean_remove_func)
    vars_of_cols_temp   = apply(temp_data, 2, var_remove_func)
    for (col_index in 1:ncol(data_w_missing)) {
      temp_col                   = temp_data[, col_index]
      normal_distn_values        = rnorm(sum(is.na(temp_col)),
                                         mean = means_of_cols_temp[col_index],
                                         sd   = sqrt(vars_of_cols_temp[col_index]))
      temp_col[is.na(temp_col)]  = normal_distn_values
      temp_data[, col_index] = temp_col
    }
    data_w_missing[min_index:max_index, ] = temp_data
    
  }
  
  return(cov(data_w_missing))
}

#' Method to assist in the parallel state splits found in the fit_regimes methods.
#'
#' @param previous_model_fits_temp rlist. The proposed parameter fit at this MCMC iteration.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param current_state integer. State at which we want to perform the split.
#' @param my_states_temp vector. Proposed regime vector.
#' @param previous_model_fits_temp rlist. The parameter fit at this MCMC iteration.
#' @param my_states vector. Current regime vector.
#'
#'
#' @noRd
#'
#'
split_parallel_helper     = function(x,
                                     previous_model_fits_temp,
                                     Z_timepoint_indices,
                                     data_points_Z,
                                     hyperparameters,
                                     first_state,
                                     my_states_temp,
                                     previous_model_fits,
                                     my_states) {
  if (x == 1) {
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state,
      my_states_temp
    )
  } else if (x == 2) {
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state + 1,
      my_states_temp
    )
  } else if (x == 3) {
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state,
      my_states
    )
    posterior_term$log_posterior  = -posterior_term$log_posterior
  } else{
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      max(my_states) + 1,
      my_states
    )
    posterior_term$log_posterior  = -posterior_term$log_posterior
  }
  return(posterior_term)
}

#' Method to assist in the parallel state merges found in the fit_regimes methods.
#'
#' @param previous_model_fits_temp rlist. The proposed parameter fit at this MCMC iteration.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param current_state integer. State at which we want to perform the merge (state before).
#' @param my_states_temp vector. Proposed regime vector.
#' @param previous_model_fits_temp rlist. The parameter fit at this MCMC iteration.
#' @param my_states vector. Current regime vector.
#'
#'
#' @noRd
#'
#'
merge_parallel_helper     = function(x,
                                     previous_model_fits_temp,
                                     Z_timepoint_indices,
                                     data_points_Z,
                                     hyperparameters,
                                     first_state,
                                     my_states_temp,
                                     previous_model_fits,
                                     my_states) {
  if (x == 1) {
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state,
      my_states_temp
    )
  } else if (x == 2) {
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      max(my_states_temp) + 1,
      my_states_temp
    )
  } else if (x == 3) {
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state,
      my_states
    )
    posterior_term$log_posterior  = -posterior_term$log_posterior
  } else{
    posterior_term                = log_Gwishart_marginals(
      previous_model_fits,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state + 1,
      my_states
    )
    posterior_term$log_posterior  = -posterior_term$log_posterior
  }
  posterior_term
}

#' Method to assist in the parallel latent data redraws found in the fit_regimes methods.
#'
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param full_data_Z matrix.  The latent data at this MCMC iteration
#' @param data_woTimeValues matrix.  The full raw data given by the user.
#' @param is_missing matrix.  Indicator matrix on whether each individual value is missing.
#' @param upper_bound_is_equal matrix.  Indicator matrix on whether each individual value achieves its upper bound.
#' @param lower_bound_is_equal matrix.  Indicator matrix on whether each individual value achieves its lower bound.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param my_states vector. Current regime vector.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param ordinal_levels vector.  Indicates the number of ordinal levels for each variable
#' @param levels_assignments vector.
#' @param discrete_levels_indicator vector.
#'
#'
#' @noRd
#'
#'
latent_data_parallel_helper = function(x,
                                       Z_timepoint_indices,
                                       full_data_Z,
                                       data_woTimeValues,
                                       is_missing,
                                       upper_bound_is_equal,
                                       lower_bound_is_equal,
                                       previous_model_fits,
                                       my_states,
                                       hyperparameters,
                                       ordinal_levels,
                                       levels_assignments,
                                       discrete_levels_indicator) {
  first_index               = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
  last_index                = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
  temp_data                 = full_data_Z[first_index:last_index, ]
  temp_raw_data             = data_woTimeValues[first_index:last_index, ]
  is_missing_temp           = is_missing[first_index:last_index, ]
  upper_bound_is_equal_temp = upper_bound_is_equal[first_index:last_index, ]
  lower_bound_is_equal_temp = lower_bound_is_equal[first_index:last_index, ]
  new_latent_data           = redraw_latent_data(
    x,
    previous_model_fits,
    hyperparameters,
    temp_data,
    temp_raw_data,
    is_missing_temp,
    upper_bound_is_equal_temp,
    lower_bound_is_equal_temp,
    ordinal_levels,
    levels_assignments,
    discrete_levels_indicator,
    1
  )
  if (anyNA(new_latent_data)) {
    stop("latent data redraw didn't work on C level -- NAs persist.")
  }
  # return(new_latent_data)
  return(new_latent_data)
  
}

#' Method to assist in the parallel split-merge algorithm on the components.
#'
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param my_states vector. Current regime vector.
#' @param full_data_Z matrix.  The latent data at this MCMC iteration
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#'
#'
#' @noRd
#'
#'
splitmerge_comps_parallel_helper = function(x,
                                            Z_timepoint_indices,
                                            my_states,
                                            full_data_Z,
                                            previous_model_fits,
                                            hyperparameters) {
  first_index               = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
  last_index                = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
  temp_data                 = full_data_Z[first_index:last_index, ]
  
  new_components_for_state  = splitmerge_gibbs_comps(
    my_states,
    previous_model_fits,
    temp_data,
    x,
    hyperparameters,
    Z_timepoint_indices,
    full_data_Z
  )
  # print(new_components_for_state)
  return(
    list(
      precisions        = new_components_for_state$precisions,
      mus               = new_components_for_state$mus,
      assigns           = new_components_for_state$assigns,
      component_probs   = new_components_for_state$comp_probs,
      component_sticks  = new_components_for_state$comp_sticks
    )
  )
}

#' Method to assist in the parallel split-merge algorithm on the components.
#'
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param full_data_Z matrix.  The latent data at this MCMC iteration
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param my_states vector. Current regime vector.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#'
#'
#' @noRd
#'
#'
gibbsswap_comps_parallel_helper = function(x,
                                           Z_timepoint_indices,
                                           full_data_Z,
                                           hyperparameters,
                                           my_states,
                                           previous_model_fits) {
  first_index               = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
  last_index                = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
  temp_data                 = full_data_Z[first_index:last_index, ]
  new_components_for_state  = gibbs_swap_comps(
    temp_data,
    previous_model_fits[[x]]$cluster_assignments,
    previous_model_fits[[x]]$component_log_probs,
    previous_model_fits[[x]]$precision,
    previous_model_fits[[x]]$mu,
    hyperparameters$component_truncation,
    2
  )
  return(list(assigns           = as.vector(new_components_for_state)))
}

#' Merge-split-swap algorithm on the regime vector as outlined in Murph et al 2023.
#'
#' @param my_states vector. Current regime vector.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param linger_parameter float. Prior parameter for Markov chain probability matrix.  Larger = less likely to change states.
#' @param move_parameter float. Prior parameter for Markov chain probability matrix.  Larger = more likely to change states.
#' @param not.cont vector. Indicator vector for which columns are discrete.
#' @param transition_probabilities vector. Current estimate of unknown transition probabilities for Markov chain.
#' @param current_graph_G matrix. Current graph structure for Gaussian Graphical Model.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param n.cores integer. Number of cores available for this calculation.
#' @param allow_mixtures logical. Whether or not a mixture model should be fit.
#' @param min_regime_length integer. Gibbs sweep will not create a regime smaller than this.
#' @param verbose logical. Will print   to console if TRUE.
#'
#'
#' @noRd
#'
#'
#'
update_states_mergesplit  = function(my_states,
                                     previous_model_fits,
                                     data_points_Z,
                                     Z_timepoint_indices,
                                     linger_parameter,
                                     move_parameter,
                                     not.cont,
                                     transition_probabilities,
                                     current_graph_G,
                                     hyperparameters,
                                     n.cores,
                                     allow_mixtures = FALSE,
                                     min_regime_length = 1,
                                     regime_selection_multiplicative_prior = NULL,
                                     split_selection_multiplicative_prior = NULL,
                                       verbose = FALSE) {
  linger_parameter         = hyperparameters$alpha
  move_parameter           = hyperparameters$beta
  mu_0                     = hyperparameters$mu_0
  lambda                   = hyperparameters$lambda
  full_data_Z              = data_points_Z
  previous_G               = current_graph_G
  p                        = ncol(data_points_Z)
  zero_means               = rep(0, times = p)
  my_states_temp           = my_states
  previous_model_fits_temp = previous_model_fits
  merge_or_split           = rbinom(1, 1, 0.5)
  
  if (merge_or_split) {
    # |------------------------ Perform a split action ------------------------|
    if ((!is.null(regime_selection_multiplicative_prior)) &
        (max(my_states) > 1)) {
      if (regime_selection_multiplicative_prior <= 1) {
        stop("The regime selection multiplicative prior must be > 1.")
      }
      
      # I start by randomly selecting a regime to split according to the regime selection prior.
      possible_regimes       = 1:(max(my_states))
      regime_distribution    = ((possible_regimes - 1) / (max(my_states) - 1)) *
        (regime_selection_multiplicative_prior - 1) + 1
      sum_density_values     = 0
      rand_value             = runif(1, min = 0, max = sum(regime_distribution))
      
      for (regime_index in 1:length(regime_distribution)) {
        sum_density_values = sum_density_values + regime_distribution[regime_index]
        if (rand_value <= sum_density_values) {
          first_state                = possible_regimes[regime_index]
          log_prob_of_regime_forward = log(regime_distribution[regime_index] / sum(regime_distribution))
          break
        }
      }
      
      # Now I need to determine the return probability, given the split, that I choose this regime to merge.
      possible_regimes             = 1:(max(my_states) + 1)
      regime_distribution          = ((possible_regimes - 1) / (max(my_states) - 1)) *
        (regime_selection_multiplicative_prior - 1) + 1
      first_state                  = possible_regimes[regime_index]
      log_prob_of_regime_backward  = log(regime_distribution[regime_index] / sum(regime_distribution))
      
    } else {
      first_state                = sample(1:(max(my_states)), 1)
      # I can choose any regime to split:
      log_prob_of_regime_forward  = log(1 / max(max(my_states)))
      # I can choose any regime but the last regime to merge
      log_prob_of_regime_backward = log(1 / max(max(my_states)))
      
    }
    log_prob_of_regime_forward  = 0
    log_prob_of_regime_backward = 0
    
    first_state           = sample(1:(max(my_states)), 1)
    if (  verbose)
      cat(
        paste("-----> starting a SPLIT on state", first_state, "\n")
      )
    
    # We will begin by automatically rejecting if the state chosen only has a single
    # observation.
    if (length(which(my_states == first_state)) <= 1) {
      if (  verbose)
        cat(
          paste(
            "-----> split automatically rejected b/c only a single instance of state",
            first_state,
            "exists\n"
          )
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          accepted_item = 0,
          transition_probabilities_item = transition_probabilities,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    }
    
    if (max(my_states) == hyperparameters$regime_truncation) {
      if (  verbose)
        cat(
          paste(
            "-----> split automatically rejected b/c reached regime truncation value of ",
            hyperparameters$regime_truncation,
            '\n'
          )
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          accepted_item = 0,
          transition_probabilities_item = transition_probabilities,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    }
    
    # Otherwise, choose a random time point in this state to act as our split.
    # Omit the first time point; introduce a new state at the chosen value and
    # set all values that happen later in time, at this state, to that value.
    indices_of_regime      = which(my_states == first_state)
    components_of_regime   = previous_model_fits[[first_state]]$cluster_assignments
    split_distribution     = get_split_distribution(
      current_graph_G,
      data_points_Z,
      Z_timepoint_indices,
      indices_of_regime,
      hyperparameters,
      n.cores,
      components_of_regime,
      hyperparameters$wishart_scale_matrix,
      hyperparameters$wishart_df,
      split_selection_multiplicative_prior,
        verbose =   verbose
    )
    
    rescale_log_values     = log(max(split_distribution))
    exp_split_distribution = exp(log(split_distribution) - rescale_log_values)
    
    if (  verbose)
      cat(
        "Rescaled and transformed split distribution:\n"
      )
    if (  verbose)
      cat(exp_split_distribution)
    if (  verbose)
      cat("\n"  )
    
    suppressWarnings({
      rand_value             = runif(1,
                                     min = 0,
                                     max = sum(exp_split_distribution))
    })
    
    if (is.na(rand_value)) {
      if (  verbose) {
        cat(
          "Something was wrong with exp split distribution"
        )
        cat(indices_of_regime)
        cat("\n"  )
        cat(split_distribution)
        cat("\n"  )
        cat(rescale_log_values)
        cat("\n"  )
        cat(exp_split_distribution)
        cat("\n"  )
        cat(rand_value)
      }
    }
    sum_density_values     = 0
    for (split_index in 1:length(split_distribution)) {
      sum_density_values = sum_density_values + exp_split_distribution[split_index]
      if (rand_value <= sum_density_values) {
        index_of_split     = indices_of_regime[split_index]
        log_prob_of_split  = log(exp_split_distribution[split_index] / sum(exp_split_distribution))
        break
      }
    }
    if (  verbose)
      cat(
        paste("-----> Index of split is:", index_of_split, '\n')
      )
    
    indices_of_state           = which(my_states == first_state)
    indices_of_state_changed   = indices_of_state[which(indices_of_state >  index_of_split)]
    indices_of_state_unchanged = indices_of_state[which(indices_of_state <= index_of_split)]
    
    my_states_temp[indices_of_state_changed] = max(my_states) + 1
    
    if ((length(which(my_states_temp == 1)) < min_regime_length)) {
      if (  verbose)
        cat(
          paste(
            "-----> split automatically rejected b/c it would create a regime less than the minimum regime length.\n"
          )
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          accepted_item = 0,
          transition_probabilities_item = transition_probabilities,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    }
    
    # Now that we know where our split is, grab the data values associated with each new group.
    # Find the first and last row of the data matrix the is included in this regime.
    Z_index_of_regime_start_changed     = Z_timepoint_indices[[min(indices_of_state_changed)]]$timepoint_first_index
    Z_index_of_regime_end_changed       = Z_timepoint_indices[[max(indices_of_state_changed)]]$timepoint_last_index
    Z_values_changed                    = data_points_Z[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    upper_bounds_changed                = hyperparameters$upper_bound_is_equal[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    lower_bounds_changed                = hyperparameters$lower_bound_is_equal[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    is_missing_changed                  = hyperparameters$is_missing[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    
    #
    Z_index_of_regime_start_unchanged   = Z_timepoint_indices[[min(indices_of_state_unchanged)]]$timepoint_first_index
    Z_index_of_regime_end_unchanged     = Z_timepoint_indices[[max(indices_of_state_unchanged)]]$timepoint_last_index
    Z_values_unchanged                  = data_points_Z[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    upper_bounds_unchanged              = hyperparameters$upper_bound_is_equal[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    lower_bounds_unchanged              = hyperparameters$lower_bound_is_equal[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    is_missing_unchanged                = hyperparameters$is_missing[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    
    # Gather all the data into one matrix as well.
    z_values_of_full_group_nochange = rbind(Z_values_unchanged, Z_values_changed)
    upper_bounds_full_group         = rbind(upper_bounds_unchanged, upper_bounds_changed)
    lower_bounds_full_group         = rbind(lower_bounds_unchanged, lower_bounds_changed)
    is_missing_full_group           = rbind(is_missing_unchanged, is_missing_changed)
    temp_n_full_group_nochange      = nrow(z_values_of_full_group_nochange)
    temp_n_changed                  = nrow(Z_values_changed)
    temp_n_unchanged                = nrow(Z_values_unchanged)
    
    # I need the new components for the split before I redraw the mixture parameters.
    log_forward_prob_component_split   = 0
    all_components                     = previous_model_fits_temp[[first_state]]$cluster_assignments
    comps_after                        = all_components[(Z_index_of_regime_end_unchanged -
                                                           Z_index_of_regime_start_unchanged + 2):length(all_components)]
    comps_before                       = all_components[1:(Z_index_of_regime_end_unchanged -
                                                             Z_index_of_regime_start_unchanged + 1)]
    
    previous_model_fits_temp[[first_state]]$cluster_assignments         = comps_before
    previous_model_fits_temp[[max(my_states_temp)]]$cluster_assignments = comps_after
    
    relabel_fits                                             = shift_states(my_states_temp, previous_model_fits_temp)
    previous_model_fits_temp                                 = relabel_fits$previous_model_fits_item
    my_states_temp                                           = relabel_fits$my_states_item
    
    redraw_states = c(max(my_states_temp), first_state)
    for (regime_index in redraw_states) {
      previous_model_fits_temp             = redraw_mixture_parameters(
        my_states_temp,
        regime_index,
        previous_model_fits_temp,
        data_points_Z,
        Z_timepoint_indices,
        linger_parameter,
        move_parameter,
        not.cont,
        current_graph_G,
        hyperparameters
      )
    }
    log_leftover_value_betaterms = calc_regimes_log_prob(my_states_temp, transition_probabilities) -
      calc_regimes_log_prob(my_states, transition_probabilities)
    
    # I now have new component assignments and new Mixture Model Parameters.  I need to calculate the return probability for the merged state.
    # I've already handled this return probability for the parameters -- I also need it for the component assignments.
    log_backward_prob_component_split = 0
    log_comp_prob_and_assignment_move = 0
    
    myfunc_5 = function(x) {
      return(
        split_parallel_helper(
          x,
          previous_model_fits_temp,
          Z_timepoint_indices,
          data_points_Z,
          hyperparameters,
          first_state,
          my_states_temp,
          previous_model_fits,
          my_states
        )
      )
    }
    
    models_temp = mclapply(1:4, myfunc_5)
    log_leftover_value_wishartterms     = 0
    laplace_nonconverge_flag            = 0
    
    for (list_index in 1:4) {
      log_leftover_value_wishartterms   = log_leftover_value_wishartterms + models_temp[[list_index]]$log_posterior
      laplace_nonconverge_flag          = max(laplace_nonconverge_flag,
                                              models_temp[[list_index]]$nonconverge_flag)
    }
    
    log_leftover_value_mu_distn     = log_mu_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state,
      my_states_temp
    ) +
      log_mu_marginals(
        previous_model_fits_temp,
        Z_timepoint_indices,
        data_points_Z,
        hyperparameters,
        first_state + 1,
        my_states_temp
      ) -
      log_mu_marginals(
        previous_model_fits,
        Z_timepoint_indices,
        data_points_Z,
        hyperparameters,
        first_state,
        my_states
      ) -
      log_mu_marginals(
        previous_model_fits,
        Z_timepoint_indices,
        data_points_Z,
        hyperparameters,
        max(my_states) + 1,
        my_states
      )
    
    pseudoprior_values              = 0
    
    # Now, either accept this split, or reject.
    M               = length(unique(my_states))
    log_likelihood_new_model = get_mixture_log_density(first_state,
                                                       previous_model_fits_temp,
                                                       Z_values_unchanged) + get_mixture_log_density(first_state + 1,
                                                                                                     previous_model_fits_temp,
                                                                                                     Z_values_changed)
    log_likelihood_old_model = get_mixture_log_density(first_state,
                                                       previous_model_fits,
                                                       z_values_of_full_group_nochange)
    MH_ratio        =  -log_prob_of_split + log_leftover_value_betaterms + log_prob_of_regime_backward - log_prob_of_regime_forward +
      log_leftover_value_wishartterms  + log_leftover_value_mu_distn -
      log_forward_prob_component_split + log_backward_prob_component_split + log_comp_prob_and_assignment_move + pseudoprior_values
    
    if (is.infinite(log_forward_prob_component_split - log_backward_prob_component_split) |
        is.na(log_forward_prob_component_split - log_backward_prob_component_split)) {
      stop("Getting some weird values for the component reassignments")
    }
    
    log_random_unif = log(runif(1))
    if (  verbose)
      cat(
        paste("-----> Metropolis-Hastings Ratio:", MH_ratio, '\n')
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from beta parameters:",
          log_leftover_value_betaterms,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from G wishart normalizing terms:",
          log_leftover_value_wishartterms,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this mixture component swap:",
          log_forward_prob_component_split - log_backward_prob_component_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "-----------> Portion of this from forward prob:",
          log_forward_prob_component_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "-----------> Portion of this from backwards prob:",-log_backward_prob_component_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "-----------> Portion of this from comp prob redraw and component vector likelihood:",
          log_comp_prob_and_assignment_move,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from proposal distribution:",-log_prob_of_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from mu distn normalizing terms:",
          log_leftover_value_mu_distn,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Value of the new likelihood:",
          log_likelihood_new_model,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Value of the old likelihood:",-log_likelihood_old_model,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Total value of likelihood ratio:",
          log_likelihood_new_model - log_likelihood_old_model,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from regime selection:",-log(M + 1) + log(M),
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste("--------> Pseudoprior terms:", pseudoprior_values, '\n')
      )
    if (  verbose)
      cat(
        paste("--------> Entire log MH ratio:", MH_ratio, '\n')
      )
    if (  verbose)
      cat(
        paste("-----> Probability of accepting SPLIT:", min(exp(MH_ratio), 1), '\n')
      )
    if (  verbose)
      cat(
        paste("-----> Log uniform random draw:", log_random_unif, '\n')
      )
    
    if (is.nan(MH_ratio)) {
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      if (  verbose)
        cat(
          "-----> completed the SPLIT (did not accept)\n"
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 0,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    } else if (laplace_nonconverge_flag) {
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      if (  verbose)
        cat(
          "-----> completed the SPLIT (did not accept because the laplace algorithm failed to converge) \n"
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 0,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    }
    
    if (log_random_unif <= MH_ratio) {
      if (  verbose)
        cat("-----> completed the SPLIT (accepted) \n")
      my_states                = my_states_temp
      previous_model_fits      = previous_model_fits_temp
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 1,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = log_leftover_value_wishartterms
        )
      )
    } else {
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      if (  verbose)
        cat(
          "-----> completed the SPLIT (did not accept)\n"
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 0,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = log_leftover_value_wishartterms
        )
      )
    }
  } else {
    # |------------------------ Perform a merge action ------------------------|
    # We will begin by automatically rejecting if there is only one state, and
    # therefore no possible merges.
    if (max(my_states) == 1) {
      if (  verbose)
        cat(
          "Cannot perform a merge b/c there is only one state! (reject automatically)\n"
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          accepted_item = 0,
          transition_probabilities_item = transition_probabilities,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    }
    # Otherwise, choose the first state randomly.
    if ((!is.null(regime_selection_multiplicative_prior)) &
        (max(my_states) > 2)) {
      if (regime_selection_multiplicative_prior <= 1) {
        stop("The regime selection multiplicative prior must be > 1.")
      }
      
      # I start by randomly selecting a regime to split according to the regime selection prior.
      possible_regimes       = 1:(max(my_states) - 1)
      regime_distribution    = ((possible_regimes - 1) / (max(possible_regimes) - 1)) *
        (regime_selection_multiplicative_prior - 1) + 1
      sum_density_values     = 0
      rand_value             = runif(1, min = 0, max = sum(regime_distribution))
      for (regime_index in 1:length(regime_distribution)) {
        sum_density_values = sum_density_values + regime_distribution[regime_index]
        if (rand_value <= sum_density_values) {
          first_state                = possible_regimes[regime_index]
          log_prob_of_regime_forward = log(regime_distribution[regime_index] / sum(regime_distribution))
          break
        }
      }
      
    } else {
      first_state                = sample(1:(max(my_states)), 1)
      # I can choose any regime to split:
      log_prob_of_regime_forward  = log(1 / max(max(my_states)))
      # I can choose any regime but the last regime to merge
      log_prob_of_regime_backward = log(1 / max(max(my_states)))
      
    }
    log_prob_of_regime_forward  = 0
    log_prob_of_regime_backward = 0
    
    first_state                              = sample(1:(max(my_states) -
                                                           1), 1)
    if (  verbose)
      cat(
        paste("-----> starting a MERGE on state", first_state, '\n')
      )
    
    # Merge this state with the state above it, making them all equal to first_state.
    indices_of_state_changed                 = which(my_states == (first_state +
                                                                     1))
    indices_of_state_unchanged               = which(my_states == (first_state))
    
    my_states_temp[min(indices_of_state_changed):max(indices_of_state_changed)] = my_states_temp[min(indices_of_state_changed):max(indices_of_state_changed)] - 1
    indices_of_merged_state                                                     = which(my_states_temp == first_state)
    
    # Now I need to calculate the split distribution to get the return probability.
    components_of_regime   = c(
      previous_model_fits[[first_state]]$cluster_assignments,
      previous_model_fits[[first_state + 1]]$cluster_assignments
    )
    split_distribution     = get_split_distribution(
      current_graph_G,
      data_points_Z,
      Z_timepoint_indices,
      indices_of_merged_state,
      hyperparameters,
      n.cores,
      components_of_regime,
      hyperparameters$wishart_scale_matrix,
      hyperparameters$wishart_df,
      split_selection_multiplicative_prior,
        verbose =   verbose
    )
    
    rescale_log_values     = log(max(split_distribution))
    exp_split_distribution = exp(log(split_distribution) - rescale_log_values)
    
    
    if (  verbose)
      cat(
        "Rescaled and transformed split distribution:\n"
      )
    if (  verbose)
      cat(exp_split_distribution)
    if (  verbose)
      cat("\n"  )
    suppressWarnings({
      rand_value             = runif(1,
                                     min = 0,
                                     max = sum(exp_split_distribution))
    })
    
    if (is.na(rand_value)) {
      if (  verbose)
        cat(
          "Something was wrong with exp split distribution\n"
        )
      if (  verbose)
        cat(split_distribution)
      if (  verbose)
        cat("\n"  )
      if (  verbose)
        cat(rescale_log_values)
      if (  verbose)
        cat("\n"  )
      if (  verbose)
        cat(exp_split_distribution)
      if (  verbose)
        cat("\n"  )
      if (  verbose)
        cat(rand_value)
      if (  verbose)
        cat("\n"  )
    }
    log_prob_of_split      = log(exp_split_distribution[length(indices_of_state_unchanged)])
    
    # I need to get the Z values for the individual groups.
    Z_index_of_regime_start_changed   = Z_timepoint_indices[[min(indices_of_state_changed)]]$timepoint_first_index
    Z_index_of_regime_end_changed     = Z_timepoint_indices[[max(indices_of_state_changed)]]$timepoint_last_index
    Z_values_changed                  = data_points_Z[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    upper_bounds_changed              = hyperparameters$upper_bound_is_equal[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    lower_bounds_changed              = hyperparameters$lower_bound_is_equal[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    is_missing_changed                = hyperparameters$is_missing[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
    Z_index_of_regime_start_unchanged = Z_timepoint_indices[[min(indices_of_state_unchanged)]]$timepoint_first_index
    Z_index_of_regime_end_unchanged   = Z_timepoint_indices[[max(indices_of_state_unchanged)]]$timepoint_last_index
    Z_values_unchanged                = data_points_Z[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    upper_bounds_unchanged            = hyperparameters$upper_bound_is_equal[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    lower_bounds_unchanged            = hyperparameters$lower_bound_is_equal[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    is_missing_unchanged              = hyperparameters$is_missing[Z_index_of_regime_start_unchanged:Z_index_of_regime_end_unchanged,]
    Z_values_of_full_group_merged     = rbind(Z_values_unchanged, Z_values_changed)
    upper_bounds_full_group           = rbind(upper_bounds_unchanged, upper_bounds_changed)
    lower_bounds_full_group           = rbind(lower_bounds_unchanged, lower_bounds_changed)
    is_missing_full_group             = rbind(is_missing_unchanged, is_missing_changed)
    
    # Information for the wishart of the two individual, initial regimes
    temp_n_changed        = nrow(Z_values_changed)
    temp_n_unchanged      = nrow(Z_values_unchanged)
    
    # Draw new parameter proposals.  The "second to last" is a random process that remains constant between merging and splitting;
    # we use the transition from the second_to_last to the new_parameter as our transition probability.
    temp_n                               = nrow(Z_values_of_full_group_merged)
    
    log_forward_prob_component_split   = 0
    previous_model_fits_temp[[first_state]]$cluster_assignments   = c(
      previous_model_fits_temp[[first_state]]$cluster_assignments,
      previous_model_fits_temp[[first_state +
                                  1]]$cluster_assignments
    )
    previous_model_fits_temp[[first_state + 1]]$cluster_assignments = NA
    
    relabel_fits                                                  = shift_states(my_states_temp, previous_model_fits_temp)
    previous_model_fits_temp                                      = relabel_fits$previous_model_fits_item
    my_states_temp                                                = relabel_fits$my_states_item
    
    previous_model_fits_temp[[max(my_states_temp) + 1]]$cluster_assignments = NA
    
    regime_states = c(first_state)
    for (regime_index in regime_states) {
      previous_model_fits_temp             = redraw_mixture_parameters(
        my_states_temp,
        regime_index,
        previous_model_fits_temp,
        data_points_Z,
        Z_timepoint_indices,
        linger_parameter,
        move_parameter,
        not.cont,
        current_graph_G,
        hyperparameters
      )
    }
    
    #   # # Whenever you split components, we will make it deterministic.
    log_backward_prob_component_split    = 0
    log_comp_prob_and_assignment_move    = 0
    
    myfunc_6 = function(x) {
      return(
        merge_parallel_helper(
          x,
          previous_model_fits_temp,
          Z_timepoint_indices,
          data_points_Z,
          hyperparameters,
          first_state,
          my_states_temp,
          previous_model_fits,
          my_states
        )
      )
    }
    models_temp = mclapply(1:4, myfunc_6)
    
    log_leftover_value_wishartterms   = 0
    laplace_nonconverge_flag          = 0
    for (list_index in 1:4) {
      log_leftover_value_wishartterms   = log_leftover_value_wishartterms + models_temp[[list_index]]$log_posterior
      laplace_nonconverge_flag          = max(laplace_nonconverge_flag,
                                              models_temp[[list_index]]$nonconverge_flag)
    }
    
    log_leftover_mu_distn           = log_mu_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      first_state,
      my_states_temp
    ) +
      log_mu_marginals(
        previous_model_fits_temp,
        Z_timepoint_indices,
        data_points_Z,
        hyperparameters,
        max(my_states_temp) + 1,
        my_states_temp
      ) -
      log_mu_marginals(
        previous_model_fits,
        Z_timepoint_indices,
        data_points_Z,
        hyperparameters,
        first_state,
        my_states
      ) -
      log_mu_marginals(
        previous_model_fits,
        Z_timepoint_indices,
        data_points_Z,
        hyperparameters,
        first_state + 1,
        my_states
      )
    
    log_leftover_value_betaterms    = calc_regimes_log_prob(my_states_temp, transition_probabilities) -
      calc_regimes_log_prob(my_states, transition_probabilities)
    
    pseudoprior_values              = 0
    
    
    # Now, either accept this split, or reject.
    M                     = length(unique(my_states))
    MH_ratio              = log_prob_of_split + log_leftover_value_betaterms + #+ log(M) - log(M-1)
      log_leftover_value_wishartterms + log_leftover_mu_distn + log_prob_of_regime_backward - log_prob_of_regime_forward - #+ log_wishart_prior_term
      log_forward_prob_component_split + log_backward_prob_component_split + log_comp_prob_and_assignment_move + pseudoprior_values
    if (is.infinite(log_forward_prob_component_split - log_backward_prob_component_split) |
        is.na(log_forward_prob_component_split - log_backward_prob_component_split)) {
      if (  verbose)
        cat(
          "We're getting some weird values for the component reassignments\n"
        )
      MH_ratio      = -Inf
    }
    
    if (  verbose)
      cat(
        paste("-----> Metropolis-Hastings Ratio:", MH_ratio, '\n')
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from beta parameters:",
          log_leftover_value_betaterms,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from G wishart parameters:",
          log_leftover_value_wishartterms,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this mixture component swap:",
          log_forward_prob_component_split - log_backward_prob_component_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "-----------> Portion of this from forward prob:",
          log_forward_prob_component_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "-----------> Portion of this from backwards prob:",-log_backward_prob_component_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "-----------> Portion of this from comp prob redraw and component vector likelihood:",
          log_comp_prob_and_assignment_move,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from proposal distribution:",
          log_prob_of_split,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from mu distn normalizing terms:",
          log_leftover_mu_distn,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion of this from regime selection:",
          log(M) - log(M - 1),
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste(
          "--------> Portion from psuedoprior:",
          pseudoprior_values,
          '\n'
        )
      )
    if (  verbose)
      cat(
        paste("--------> The full MH ratio is:", MH_ratio, '\n')
      )
    if (  verbose)
      cat(
        paste("-----> Probability of accepting MERGE:", min(exp(MH_ratio), 1), '\n')
      )
    log_random_unif       = log(runif(1))
    if (  verbose)
      cat(
        paste("-----> Log uniform random draw:", log_random_unif, '\n')
      )
    
    if (is.na(MH_ratio)) {
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      if (  verbose)
        cat(
          "-----> completed the MERGE (did not accept)\n"
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 0,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    } else if (laplace_nonconverge_flag) {
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      if (  verbose)
        cat(
          "-----> completed the MERGE (did not accept because laplace approx failed to converge)\n"
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 0,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = NA
        )
      )
    }
    
    if (log_random_unif <= MH_ratio) {
      if (  verbose)
        cat("-----> completed the MERGE (accepted) \n")
      # Accept the merge!
      # Finish by shifting the states.
      data_shifted             = shift_states(my_states_temp, previous_model_fits_temp)
      my_states                = my_states_temp
      previous_model_fits      = previous_model_fits_temp
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 1,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = log_leftover_value_wishartterms
        )
      )
    } else {
      if (  verbose)
        cat(
          "-----> completed the MERGE (did not accept)\n"
        )
      if (  verbose)
        cat(
          c(
            "[1]  -----> The current probabities are:",
            paste(sapply(
              transition_probabilities, paste, collapse = ' '
            )),
            '\n'
          )
        )
      return(
        list(
          my_states_item = my_states,
          previous_model_fits_item = previous_model_fits,
          transition_probabilities_item = transition_probabilities,
          accepted_item = 0,
          full_data_Z_item = data_points_Z,
          G_Wish_portion_of_MH = log_leftover_value_wishartterms
        )
      )
    }
  }
}

#' Relabels mixture model components.
#'
#' @param current_state integer. Regime for which we want to reassign components.
#' @param previous_model_fits rlist. Current MCMC samples parameter fits.  This includes component assignments.
#'
#' @noRd
#'
shift_components          = function(current_state, previous_model_fits) {
  # Reassign components within a regime to drop any components for which no data is assigned.
  cluster_assignments      = previous_model_fits[[current_state]]$cluster_assignments
  current_state_mus        = previous_model_fits[[current_state]]$mu
  current_state_precisions = previous_model_fits[[current_state]]$precision
  current_state_comp_probs = previous_model_fits[[current_state]]$component_log_probs
  current_state_sticks     = previous_model_fits[[current_state]]$component_sticks
  
  running_index = 0
  unique_components = 1:max(cluster_assignments)
  f = function(x) {
    sum(cluster_assignments == x)
  }
  counts_of_each    = sapply(unique_components, f)
  comps_to_iterate  = unique_components[order(counts_of_each, decreasing = TRUE)]
  
  for (actual_index in comps_to_iterate) {
    if (sum(cluster_assignments == actual_index) == 0) {
      # Then this is a state that need to be removed.
      # I think, maybe with this setup, I just do nothing here.
    } else {
      running_index                             = running_index + 1
      indices_of_component                      = which(previous_model_fits[[current_state]]$cluster_assignments == actual_index)
      if (length(current_state_mus) >= actual_index) {
        current_state_mus[[running_index]]        = current_state_mus[[actual_index]]
        current_state_precisions[[running_index]] = current_state_precisions[[actual_index]]
        current_state_comp_probs[running_index]   = current_state_comp_probs[actual_index]
        current_state_sticks[running_index]       = current_state_sticks[actual_index]
      } else {
        current_state_mus[[running_index]]        = NA
        current_state_precisions[[running_index]] = NA
      }
      
      cluster_assignments[indices_of_component] = running_index
    }
  }
  
  previous_model_fits[[current_state]]$mu                  = current_state_mus
  previous_model_fits[[current_state]]$precision           = current_state_precisions
  previous_model_fits[[current_state]]$cluster_assignments = cluster_assignments
  previous_model_fits[[current_state]]$component_log_probs = current_state_comp_probs
  previous_model_fits[[current_state]]$component_sticks    = current_state_sticks
  
  return(previous_model_fits)
}


#' There are leftover G-wishart normalizing terms when performing the Metropolis-Hastings
#' updating step in the merge-split algorithm for the regime vector.  This method calculates
#' these values.  See Murph et al 2023 for more details.
#'
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param regime_of_interest integer. Regime for which to calculate marginal.
#' @param state_vector vector. Current regime vector.
#'
#'
#' @noRd
#'
#'
log_Gwishart_marginals    = function(previous_model_fits,
                                     Z_timepoint_indices,
                                     data_points_Z,
                                     hyperparameters,
                                     regime_of_interest,
                                     state_vector) {
  # Get the 'leftover' G-Wishart normalizing terms that show up in the MH ratio of the regime-level merge-split
  # algorithm.  I used to just do this in the method, but it has become more complicated with the precense of
  # multiple components, so I moved it to its own method.
  
  # Grab data from regime_of_interest.
  indices_of_state     = which(state_vector == regime_of_interest)
  components_of_state  = previous_model_fits[[regime_of_interest]]$cluster_assignments
  mu_0                 = hyperparameters$mu_0
  p                    = hyperparameters$p
  graph_structure      = hyperparameters$G
  laplace_failure_flag = 0
  lambda               = hyperparameters$lambda
  scale_matrix         = hyperparameters$wishart_scale_matrix
  df                   = floor(hyperparameters$wishart_df)
  log_wish_prior_marg  = BDgraph::gnorm(hyperparameters$G,
                                        b = df,
                                        D = scale_matrix,
                                        iter = 1000)
  # Iterate over components.  If there is no data for a component, add in a prior term.
  if (length(indices_of_state) == 0) {
    # In this case, there are no data assigned to this regime.  So, just return log_prior terms equal to the
    # number of allowable components.
    return(
      list(
        log_posterior = hyperparameters$component_truncation * log_wish_prior_marg,
        nonconverge_flag = 0
      )
    )
  }
  
  Z_index_of_regime_start_changed   = Z_timepoint_indices[[min(indices_of_state)]]$timepoint_first_index
  Z_index_of_regime_end_changed     = Z_timepoint_indices[[max(indices_of_state)]]$timepoint_last_index
  Z_values                          = data_points_Z[Z_index_of_regime_start_changed:Z_index_of_regime_end_changed,]
  
  # Iterate over the number of possible components and add in the corresponding marginal terms.
  log_Gwishart_marginal_vals = 0
  for (comp_index in 1:hyperparameters$component_truncation) {
    indices_of_comp = which(components_of_state == comp_index)
    if (length(indices_of_comp) == 0) {
      log_Gwishart_marginal_vals = log_Gwishart_marginal_vals + log_wish_prior_marg
      # log_Gwishart_marginal_vals = log_Gwishart_marginal_vals -0.5*(df-1+p)*log(abs(det(scale_matrix))) + (df-1+p)*p*0.5*log(2) + lmvgamma(0.5*(df-1+p),p)
      next
    } else if (length(indices_of_comp) == 1) {
      # print("getting the marginal for a single obs")
      temp_Z_values = Z_values[indices_of_comp,]
      temp_n        = 1
      mu            = matrix(Z_values, nrow = p, ncol = 1)
      xbar          = matrix(rep(mu, times = temp_n), temp_n, p, byrow =
                               T)
      nS2           = matrix(0, p, p)
      
      posterior_scale_matrix = nS2 + scale_matrix + (temp_n * lambda / (temp_n +
                                                                          lambda)) * (mu - mu_0) %*% t(mu - mu_0)
    } else {
      # print("getting the marginal for a multiple obs")
      temp_Z_values = Z_values[indices_of_comp,]
      temp_n        = nrow(temp_Z_values)
      mu            = apply(temp_Z_values, 2, mean)
      xbar          = matrix(rep(mu, times = temp_n), temp_n, p, byrow =
                               T)
      nS2           = t(temp_Z_values - xbar) %*% (temp_Z_values - xbar)
      
      posterior_scale_matrix = nS2 + scale_matrix + (temp_n * lambda / (temp_n +
                                                                          lambda)) * (mu - mu_0) %*% t(mu - mu_0)
    }
    
    # Note that delta and n get added in the C++ code.
    if ((df + temp_n) < 30) {
      log_Gwishart_marginal_vals = log_Gwishart_marginal_vals + BDgraph::gnorm(
        hyperparameters$G,
        b = df + temp_n,
        D = posterior_scale_matrix,
        iter = 1000
      )
    } else {
      posterior_list             = log_normalizing_g_wishart_posterior_laplace(as.matrix(unclass(graph_structure)),
                                                                               posterior_scale_matrix,
                                                                               df,
                                                                               temp_n,
                                                                               hyperparameters$p)
      laplace_failure_flag       = laplace_failure_flag + posterior_list[["nonconverge_flag"]]
      log_Gwishart_marginal_vals = log_Gwishart_marginal_vals + posterior_list[["log_posterior"]]
    }
  }
  return(list(
    log_posterior = log_Gwishart_marginal_vals,
    nonconverge_flag = (laplace_failure_flag > 0)
  ))
}


#' There are leftover Normal-Inverse-Wishart terms when performing the Metropolis-Hastings
#' updating step in the merge-split algorithm for the regime vector.  This method calculates
#' these values.  See Murph et al 2023 for more details.
#'
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param regime_of_interest integer. Regime for which to calculate marginal.
#' @param state_vector vector. Current regime vector.
#'
#'
#' @noRd
#'
#'
log_mu_marginals          = function(previous_model_fits,
                                     Z_timepoint_indices,
                                     data_points_Z,
                                     hyperparameters,
                                     regime_of_interest,
                                     state_vector) {
  # Get the 'leftover' G-Wishart normalizing terms that show up in the MH ratio of the regime-level merge-split
  # algorithm.  I used to just do this in the method, but it has become more complicated with the presense of
  # multiple components, so I moved it to its own method.
  
  # Grab data from regime_of_interest.
  indices_of_state    = which(state_vector == regime_of_interest)
  components_of_state = previous_model_fits[[regime_of_interest]]$cluster_assignments
  p                   = hyperparameters$p
  lambda              = hyperparameters$lambda
  df                  = floor(hyperparameters$wishart_df)
  
  # Iterate over components.  If there is no data for a component, add in a prior term.
  if (length(indices_of_state) == 0) {
    # In this case, there are no data assigned to this regime.  So, just return log_prior terms equal to the
    # number of allowable components.
    return(-hyperparameters$component_truncation * 0.5 * p * log(lambda))
  }
  
  # Iterate over the number of possible components and add in the corresponding marginal terms.
  log_mu_marginal_vals = 0
  for (comp_index in 1:hyperparameters$component_truncation) {
    indices_of_comp = which(components_of_state == comp_index)
    if (length(indices_of_comp) == 0) {
      log_mu_marginal_vals = log_mu_marginal_vals - 0.5 * p * log(lambda)
      next
    } else if (length(indices_of_comp) == 1) {
      temp_n        = 1
      log_mu_marginal_vals = log_mu_marginal_vals - 0.5 * p * log(lambda +
                                                                    1)
      
    } else {
      temp_n        = length(indices_of_comp)
      log_mu_marginal_vals = log_mu_marginal_vals - 0.5 * p * log(lambda +
                                                                    temp_n)
    }
  }
  return(log_mu_marginal_vals)
}


#' The algorithms for resampling the regime vector require a
#' relabeling at the end.  This method does that relabeling.
#'
#' @param my_states vector. Current regime vector.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#'
#'
#' @noRd
#'
#'
shift_states              = function(my_states, previous_model_fits) {
  current_state            = 1
  previous_model_fits_temp = previous_model_fits
  my_states_temp           = rep(0, times = length(my_states))
  for (index in 2:(length(my_states))) {
    # Fill in the first state, if we are just starting.
    if ((index - 1) == 1) {
      previous_model_fits_temp[[current_state]][["precision"]]           = previous_model_fits[[my_states[[index -
                                                                                                             1]]]][["precision"]]
      previous_model_fits_temp[[current_state]][["mu"]]                  = previous_model_fits[[my_states[[index -
                                                                                                             1]]]][["mu"]]
      previous_model_fits_temp[[current_state]][["cluster_assignments"]] = previous_model_fits[[my_states[[index -
                                                                                                             1]]]][["cluster_assignments"]]
      previous_model_fits_temp[[current_state]][["component_log_probs"]] = previous_model_fits[[my_states[[index -
                                                                                                             1]]]][["component_log_probs"]]
      previous_model_fits_temp[[current_state]][["component_sticks"]] = previous_model_fits[[my_states[[index -
                                                                                                          1]]]][["component_sticks"]]
      my_states_temp[index - 1]                                            = 1
    }
    
    # Check to see if there has been a change of state.  If so, increase the current state by one and fill in.
    if (my_states[index - 1] != my_states[index]) {
      current_state                                                      = current_state + 1
      previous_model_fits_temp[[current_state]][["precision"]]           = previous_model_fits[[my_states[[index]]]][["precision"]]
      previous_model_fits_temp[[current_state]][["mu"]]                  = previous_model_fits[[my_states[[index]]]][["mu"]]
      previous_model_fits_temp[[current_state]][["cluster_assignments"]] = previous_model_fits[[my_states[[index]]]][["cluster_assignments"]]
      previous_model_fits_temp[[current_state]][["component_log_probs"]] = previous_model_fits[[my_states[[index]]]][["component_log_probs"]]
      previous_model_fits_temp[[current_state]][["component_sticks"]] = previous_model_fits[[my_states[[index]]]][["component_sticks"]]
    }
    my_states_temp[index]                                                = current_state
  }
  return(
    list(
      my_states_item = my_states_temp,
      previous_model_fits_item = previous_model_fits_temp
    )
  )
}


#' Re sampling the unknown transition probabilities for the Markov chain.
#'
#' @param my_states vector. Current regime vector.
#' @param my_alpha Prior parameter for Markov chain probability matrix.  Larger = less likely to change states.
#' @param my_beta float. Prior parameter for Markov chain probability matrix.  Larger = more likely to change states.
#' @param n.cores integer. Number of cores available for this calculation.
#'
#'
#' @noRd
#'
#'
redraw_transition_probs   = function(my_states, my_alpha, my_beta, n.cores) {
  probs_temp = c()
  for (x in 1:length(my_states)) {
    num_observed_of_state = sum(my_states == x)
    if (x < max(my_states)) {
      probs_temp = c(probs_temp,
                     rbeta(1, my_alpha + num_observed_of_state - 1, my_beta + 1))
    } else if (x == max(my_states)) {
      probs_temp = c(probs_temp,
                     rbeta(1, my_alpha + num_observed_of_state - 1, my_beta))
    } else {
      probs_temp = c(probs_temp, rbeta(1, my_alpha, my_beta))
    }
  }
  return(probs_temp)
}


#' Resample hyperparameters according to a Bayesian hierarchical framework.
#'
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param transition_probs vector. Current estimate of unknown transition probabilities for Markov chain.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param my_states vector. Current regime vector.
#' @param verbose logical. Will print   to console if TRUE.
#'
#'
#' @noRd
#'
#'
redraw_hyperparameters    = function(hyperparameters,
                                     transition_probs,
                                     previous_model_fits,
                                     my_states,
                                       verbose = FALSE) {
  rate_of_sampling_distn <- NULL
  # | ---------------------- Gather Data ---------------------------- |
  threshold      = 1e-8
  alpha          = hyperparameters$alpha
  beta           = hyperparameters$beta
  # Beta Hyperparameters for alpha, beta.  Assumes a gamma prior on each param.
  shape_beta  = hyperparameters$beta_hyperparameter
  rate_beta   = 1
  shape_alpha = hyperparameters$alpha_hyperparameter
  rate_alpha  = 1
  # Normal Hyperparameters for m
  p       = hyperparameters$p
  m_0     = rep(0, times = p)
  prec_0  = diag(p)
  # Gamma hyperparameters for lambda
  rate_0  = 10
  shape_0 = hyperparameters$lambda_hyperparameter_shape * 10
  
  # G-Wishart hyperparameters for Wishart scale matrix
  ## This update is NEARLY conjugate
  scale_hyperprior_matrix  = hyperparameters$hyperprior_scale_matrix
  df_hyperprior            = hyperparameters$hyperprior_b
  
  # Gamma hyperparameters for degrees of freedom
  spread_of_df_random_walk = 1
  rate_prior_for_wish_df   = 10
  shape_prior_for_wish_df  = hyperparameters$original_wishart_df * 10
  
  # | --------------------- Redraw Sticks from Truncated Dirichlet Process ------------------------- |
  temp_component_sticks = rep(0, times = hyperparameters$component_truncation)
  for (i in 1:hyperparameters$regime_truncation) {
    cluster_assignments = previous_model_fits[[i]]$cluster_assignments
    if (anyNA(cluster_assignments)) {
      # This means that there are no data assigned to this regime.
      temp_component_sticks   = rbeta(hyperparameters$component_truncation,
                                      1,
                                      hyperparameters$dirichlet_prior)
    } else {
      for (stick_index in 1:hyperparameters$component_truncation) {
        if (stick_index %in% cluster_assignments) {
          temp_component_sticks[stick_index] = rbeta(
            1,
            1 + sum(cluster_assignments == stick_index),
            hyperparameters$dirichlet_prior + sum(cluster_assignments > stick_index)
          )
        } else {
          temp_component_sticks[stick_index] = rbeta(1, 1, hyperparameters$dirichlet_prior)
        }
      }
    }
    previous_model_fits[[i]][["component_log_probs"]] = sticks_to_log_probs(temp_component_sticks)
    previous_model_fits[[i]][["component_sticks"]]    = temp_component_sticks
  }
  
  # Clean and examine input data.
  if ((alpha <= 0) | (beta <= 0)) {
    return(NA)
  }
  transition_probabilities = transition_probs[which(!is.na(transition_probs))]
  
  # | --------------------- Redraw Alpha, Beta ------------------------- |
  # Propose new Alpha, Beta from prior.  Calculate the MH ratio.
  log_prob_mass  = 0
  perturb_alpha  = rnorm(1, 0, 0.1)
  perturb_beta   = rnorm(1, 0, 0.1)
  new_alpha      = exp(log(hyperparameters$alpha) + perturb_alpha)
  new_beta       = exp(log(hyperparameters$beta) + perturb_beta)
  
  for (value in transition_probabilities) {
    log_prob_mass = log_prob_mass + dbeta(value, new_alpha, new_beta, log = TRUE) -
      dbeta(value, alpha, beta, log = TRUE)
    
  }
  log_prob_mass = log_prob_mass + dgamma(new_alpha, shape_alpha, rate = rate_alpha, log = TRUE) +
    dgamma(new_beta, shape_beta, rate = rate_beta, log = TRUE) -
    dgamma(alpha, shape_alpha, rate = rate_alpha, log = TRUE) -
    dgamma(beta, shape_beta, rate = rate_beta, log = TRUE) +
    log(new_alpha) + log(new_beta) - log(alpha) - log(beta)
  
  
  if (is.na(log_prob_mass)) {
    #
    if (  verbose)
      cat("Got an NA for Log Probability Mass.\n")
    if (  verbose)
      cat(new_alpha)
    if (  verbose)
      cat("\n"  )
    if (  verbose)
      cat(new_beta)
    if (  verbose)
      cat("\n"  )
    if (  verbose)
      cat(transition_probabilities)
    if (  verbose)
      cat("\n"  )
    log_prob_mass = -Inf
  }
  
  # Accept or Reject based on the MH ratio.
  rand_value   = log(runif(1))
  if (rand_value <= log_prob_mass) {
    # ACCEPT
    hyperparameters$alpha = new_alpha
    hyperparameters$beta  = new_beta
    accepted   = 1
  } else {
    # REJECT
    accepted   = 0
  }
  
  # | ------------------------ Redraw m ------------------------------- |
  ## Conjugate update.  We put a simple normal distribution on m.
  different_states       = unique(my_states)
  components_in_state    = hyperparameters$component_truncation
  
  sum_precision_matrices     = matrix(0, p, p)
  sum_precision_times_mu     = matrix(0, nrow = p)
  sum_prec_matries_wo_lambda = matrix(0, p, p)
  sum_log_det_prec_matrices  = rep(0, times = hyperparameters$regime_truncation)
  
  max_states       = max(my_states)
  count_components = 0
  for (index in 1:max(my_states)) {
    precisions_in_state                 = previous_model_fits[[index]]$precision
    mus_in_state                        = previous_model_fits[[index]]$mu
    comps_in_state                      = previous_model_fits[[index]]$cluster_assignments
    for (comp_index in 1:unique(comps_in_state)) {
      count_components                    = count_components + 1
      sum_precision_matrices              = sum_precision_matrices + hyperparameters$lambda * precisions_in_state[[comp_index]]
      sum_prec_matries_wo_lambda          = sum_prec_matries_wo_lambda + precisions_in_state[[comp_index]]
      sum_log_det_prec_matrices[index]    = sum_log_det_prec_matrices[index] + log(abs(det(precisions_in_state[[comp_index]])))
      sum_precision_times_mu              = sum_precision_times_mu + hyperparameters$lambda *
        precisions_in_state[[comp_index]] %*% mus_in_state[[comp_index]]
    }
  }
  
  epsilon = 1e-10
  if (abs(rcond(prec_0 + sum_precision_matrices)) > epsilon) {
    hyperprior_mean = solve(prec_0 + sum_precision_matrices) %*% (prec_0 %*%
                                                                    m_0 + sum_precision_times_mu)
    new_mus = rmvn_Rcpp(as.double(hyperprior_mean),
                        as.double((hyperparameters$lambda) * (sum_precision_matrices +
                                                                prec_0)
                        ),
                        as.integer(hyperparameters$p))
    
    if (!anyNA(new_mus)) {
      hyperparameters$mu_0     = new_mus
    }
  }
  
  # | --------------------- Redraw Scale Matrix D ------------------------- |
  # We removed this for this version of the algorithm.  See the paper for details on possible
  # extensions in this method.
  
  # | --------------------- Redraw Degrees of Freedom b ------------------------- |
  old_df                 = hyperparameters$wishart_df
  new_df                 = rnorm(1, mean = old_df, sd = spread_of_df_random_walk)
  
  if (!is.nan(new_df)) {
    if (new_df >= 3) {
      regime_scale_matrix = hyperparameters$wishart_scale_matrix
      
      if ((new_df > 50) | (old_df > 50)) {
        posterior_term_1      = log_normalizing_g_wishart_posterior_laplace(as.matrix(unclass(hyperparameters$G)),
                                                                            regime_scale_matrix,
                                                                            old_df,
                                                                            0,
                                                                            hyperparameters$p)
        posterior_term_2      = log_normalizing_g_wishart_posterior_laplace(as.matrix(unclass(hyperparameters$G)),
                                                                            regime_scale_matrix,
                                                                            new_df,
                                                                            0,
                                                                            hyperparameters$p)
        MH_ratio              = posterior_term_1[["log_posterior"]] - posterior_term_2[["log_posterior"]]
        
        laplace_failure_flag  = max(posterior_term_1[["nonconverge_flag"]], posterior_term_1[["nonconverge_flag"]])
        
      } else {
        MH_ratio              = BDgraph::gnorm(as.matrix(unclass(hyperparameters$G)),
                                               b = old_df,
                                               regime_scale_matrix,
                                               iter = 500) -
          BDgraph::gnorm(as.matrix(unclass(hyperparameters$G)),
                         b = new_df,
                         regime_scale_matrix,
                         iter = 500)
        laplace_failure_flag  = 0
      }
      
      MH_ratio   = MH_ratio * count_components
      MH_ratio   = MH_ratio - new_df * (-0.5 * sum(sum_log_det_prec_matrices) + old_df * (-0.5 *
                                                                                            sum(sum_log_det_prec_matrices)))
      MH_ratio   = MH_ratio + dgamma(new_df,
                                     shape = shape_prior_for_wish_df,
                                     rate = rate_prior_for_wish_df,
                                     log = TRUE) -
        dgamma(old_df,
               shape = shape_prior_for_wish_df,
               rate = rate_prior_for_wish_df,
               log = TRUE)
      
      if ((log(runif(1)) <= MH_ratio) |
          (laplace_failure_flag == 0)) {
        hyperparameters$wishart_df = new_df
      } else if (laplace_failure_flag) {
        if (  verbose)
          cat(
            "MH ratio failed automatically since the flag was triggered\n"
          )
      }
    }
  } else{
    if (  verbose)
      cat("new df was nan!!\n"  )
    if (  verbose)
      cat(
        "rate of the gamma sampling distribution was likley invalid.  It was:\n"
      )
    if (  verbose)
      cat(rate_of_sampling_distn)
    if (  verbose)
      cat("\n"  )
  }
  
  # | ---------------------- Redraw lambda ---------------------------- |
  ## Using Gamma prior, this is also a conjugate update.
  trace_term             = 0
  for (index in 1:max(my_states)) {
    precisions_in_state    = previous_model_fits[[index]]$precision
    mus_in_state           = previous_model_fits[[index]]$mu
    for (comp_index in unique(components_in_state)) {
      center_term_temp     = mus_in_state[[comp_index]] - hyperparameters$mu_0
      trace_term           = trace_term  + 0.5 * t(center_term_temp) %*% precisions_in_state[[comp_index]] %*% center_term_temp
    }
  }
  posterior_rate         = trace_term + rate_0
  posterior_shape        = shape_0 + 0.5 * p * count_components
  new_lambda             = rgamma(1, shape = posterior_shape, rate = posterior_rate)
  if (!is.na(new_lambda)) {
    if (  verbose)
      cat(paste("NEW LAMBDA IS:", new_lambda, '\n'))
    if (  verbose)
      cat(
        paste("posterior rate:", posterior_rate, '\n')
      )
    if (  verbose)
      cat(
        paste("posterior shape:", posterior_shape, '\n')
      )
    hyperparameters$lambda = new_lambda
  }
  
  # | ------------------------- Done! ------------------------------- |
  return(
    list(
      previous_model_fits_item = previous_model_fits,
      hyperparameters_item = hyperparameters,
      accepted_item = accepted
    )
  )
}


#' Draw Lambda and mu from a NI-G-Wishart distribution according to Murph et al 2023.
#'
#' @param current_G matrix. Graph structure of G-Wishart.
#' @param b integer. DF of NIW.
#' @param nS2 matrix. n*empirical covariance matrix.
#' @param n integer. Number of observations.
#' @param y_bar vector. Column average of data.
#' @param mu_0 vector. NIW hyperparameter. Center of mu distribution.
#' @param lambda float. NIW hyperparameter.  Controls effect of nS2 on mu.
#' @param scale_matrix matrix. NIW hyperparameter.
#'
#'
#' @noRd
#'
#'
draw_mvn_parameters       = function(current_G,
                                     b,
                                     nS2,
                                     n,
                                     y_bar,
                                     mu_0,
                                     lambda,
                                     scale_matrix) {
  # If needed, I could speed things up with C++ here.  Maybe pass the data in instead of nS2 and n, then do all
  # the matrix algebra in C/FORTRAN.
  # | ----------- Gather Parameters -------------- |
  p              = ncol(nS2)
  y_bar_mat      = matrix(y_bar - mu_0, p, 1)
  scale_matrix_n = scale_matrix + nS2 + ((n * lambda) / (n + lambda)) * (y_bar_mat %*% t(y_bar_mat))
  mu_n           = (n * y_bar + mu_0 * lambda) / (n + lambda)
  threshold      = 1e-8
  
  # | ----------- Draw Precision -------------- |
  result = rgwish_Rcpp(
    as.double(current_G),
    as.double(scale_matrix_n),
    as.integer(b + n),
    as.integer(p),
    as.double(threshold)
  )
  flag       = result[['failed']]
  result     = result[['K']]
  
  new_K  = matrix(result, p, p)
  
  # | ----------- Draw New Mu ----------------- |
  new_mu = rmvn_Rcpp(as.double(mu_n),
                     as.double((lambda + n) * new_K), as.integer(p))
  
  new_mu = as.matrix(new_mu, p, 1)
  
  
  return(list(final_mu = new_mu, final_K = new_K))
}


#' Process to find a split point for the merge-split algorithm on the regime vector.
#' Probability of proposal is relative to the likelihood.
#'
#' @param current_G matrix. Current graph structure for Gaussian Graphical Model.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param indices_of_regime vector. Data instances associated with current regime.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param n.cores integer. Number of cores available for this calculation.
#' @param components_of_regime vector. Mixture component assignments for current regime.
#' @param scale_matrix_of_regime matrix. Hyperparameter for NI-G-W
#' @param df_of_regime float. DF of NI-G-W.
#' @param verbose logical. Will print   to console if TRUE.
#'
#'
#' @noRd
#'
#'
get_split_distribution    = function(current_G,
                                     data_points_Z,
                                     Z_timepoint_indices,
                                     indices_of_regime,
                                     hyperparameters,
                                     n.cores,
                                     components_of_regime,
                                     scale_matrix_of_regime,
                                     df_of_regime,
                                     split_selection_multiplicative_prior = NULL,
                                       verbose = FALSE) {
  # This method assumes that the number of indices in the regime is at least 2.
  # I'll likely do at least 3, since the split point for 2 is deterministic.
  # This value can be tuned for a precision vs. runtime tradeoff.
  mu_dataframe               = NULL
  mu_0                       = hyperparameters$mu_0
  b                          = df_of_regime
  p                          = hyperparameters$p
  lambda                     = hyperparameters$lambda
  scale_matrix_D             = scale_matrix_of_regime
  discretization_val         = 5
  split_distribution         = rep(0, times = (length(indices_of_regime) -
                                                 1))
  number_of_breaks           = max(floor((length(
    indices_of_regime
  ) - 1) / discretization_val), 2)
  
  # Find the first and last row of the data matrix the is included in this regime.
  Z_index_of_regime_start    = Z_timepoint_indices[[min(indices_of_regime)]]$timepoint_first_index
  Z_index_of_regime_end      = Z_timepoint_indices[[max(indices_of_regime)]]$timepoint_last_index
  
  # We will want to save the points at which we draw parameters and evaluate the likelihood.
  list_of_evaluation_points  = c()
  list_of_evaluation_indices = c()
  
  # Start by filling in the vector with likelihood values that I actually evaluate.
  for (index in 1:number_of_breaks) {
    # Let the first index evaluated be the very first index always.
    if (index == 1) {
      evaluation_point        = min(indices_of_regime)
      evaluation_index        = 1
      
      # Let the last index evaluated be the very last index always.
    } else if (index == number_of_breaks) {
      evaluation_point        = indices_of_regime[length(indices_of_regime) -
                                                    1]
      evaluation_index        = length(indices_of_regime) - 1
      
      # Otherwise, choose the index evaluated at random.
    } else {
      index_within_chunck     = sample(0:(discretization_val - 1), 1)
      evaluation_index        = index * discretization_val + index_within_chunck
      evaluation_point        = indices_of_regime[evaluation_index]
    }
    # Save at which indices we evaluate the likelihood.
    list_of_evaluation_indices = c(list_of_evaluation_indices, evaluation_index)
    
    # Split data, draw parameters, evaluate the likelihood.
    Z_index_of_eval                       = Z_timepoint_indices[[evaluation_point]]$timepoint_last_index
    data_before_full                      = data_points_Z[Z_index_of_regime_start:Z_index_of_eval,]
    likelihood_value                      = 0
    components_before                     = components_of_regime[1:(Z_index_of_eval - Z_index_of_regime_start + 1)]
    
    for (component_index in unique(components_before)) {
      data_before          = data_before_full[which(components_before == component_index),]
      n_before             = nrow(data_before)
      if (is.null(n_before)) {
        n_before                  = 1
        x_bar_before              = data_before
        x_bar_before_mat          = data_before
        nS2_before                = matrix(0, hyperparameters$p, hyperparameters$p)
      } else {
        x_bar_before         = apply(data_before, 2, mean)
        x_bar_before_mat     = matrix(rep(x_bar_before, times = n_before), n_before, p, byrow = T)
        nS2_before           = t(data_before - x_bar_before_mat) %*% (data_before -
                                                                        x_bar_before_mat)
      }
      
      parameters_before = draw_mvn_parameters(current_G,
                                              b,
                                              nS2_before,
                                              n_before,
                                              x_bar_before,
                                              mu_0,
                                              lambda,
                                              scale_matrix_D)
      
      if (is.na(max(abs(parameters_before$final_K)))) {
        likelihood_value = -1e200
        if (  verbose)
          cat(
            "SPLIT DISTRIBUTION SUGGESTED OUTLANDISH PRECISION VALUE(S)\n"
          )
      } else if ((max(abs(parameters_before$final_K)) > 1e100)) {
        likelihood_value = -1e200
        if (  verbose)
          cat(
            "SPLIT DISTRIBUTION SUGGESTED OUTLANDISH PRECISION VALUE(S)\n"
          )
      } else {
        precision   = parameters_before$final_K
        mu          = parameters_before$final_mu
        data_before       = matrix(data_before, n_before, hyperparameters$p)
        
        likelihood_before =  log_dmvnrm_arma_regular(data_before,
                                                     parameters_before$final_mu,
                                                     parameters_before$final_K)
      }
    }
    
    
    data_after_full      = data_points_Z[(Z_index_of_eval + 1):Z_index_of_regime_end, ]
    components_after     = components_of_regime[(Z_index_of_eval - Z_index_of_regime_start + 2):(length(components_of_regime))]
    
    for (component_index in unique(components_after)) {
      data_after          = data_after_full[which(components_after == component_index),]
      n_after             = nrow(data_after)
      if (is.null(n_after)) {
        n_after = 1
        x_bar_after               = data_after
        x_bar_after_mat           = data_after
        nS2_after                 = matrix(0, hyperparameters$p, hyperparameters$p)
      } else {
        x_bar_after         = apply(data_after, 2, mean)
        x_bar_after_mat     = matrix(rep(x_bar_after, times = n_after), n_after, p, byrow = T)
        nS2_after           = t(data_after - x_bar_after_mat) %*% (data_after -
                                                                     x_bar_after_mat)
      }
      
      parameters_after = draw_mvn_parameters(current_G,
                                             b,
                                             nS2_after,
                                             n_after,
                                             x_bar_after,
                                             mu_0,
                                             lambda,
                                             scale_matrix_D)
      
      if (is.na(max(abs(parameters_after$final_K)))) {
        likelihood_value = -1e200
        if (  verbose)
          cat(
            "SPLIT DISTRIBUTION SUGGESTED OUTLANDISH PRECISION VALUE(S)\n"
          )
      } else if ((max(abs(parameters_after$final_K)) > 1e100)) {
        likelihood_value = -1e200
        if (  verbose)
          cat(
            "SPLIT DISTRIBUTION SUGGESTED OUTLANDISH PRECISION VALUE(S) \n"
          )
      } else {
        precision        = parameters_after$final_K
        mu               = parameters_after$final_mu
        data_after       = matrix(data_after, n_after, hyperparameters$p)
        likelihood_after = log_dmvnrm_arma_regular(data_after,
                                                   parameters_after$final_mu,
                                                   parameters_after$final_K)
      }
    }
    
    likelihood_value = likelihood_before + likelihood_after
    # Now, linearly interpolate the values in between.
    if (index > 1) {
      index_before     = list_of_evaluation_indices[index - 1]
      index_difference = evaluation_index - index_before
    } else {
      index_difference = 1
    }
    
    if ((index > 1) &
        (length(split_distribution) > 1) & (index_difference > 1)) {
      split_distribution[evaluation_index] = likelihood_value
      index_before = list_of_evaluation_indices[index - 1]
      
      previous_value = split_distribution[index_before]
      slope     = (likelihood_value - previous_value) /
        (evaluation_index - index_before)
      intercept = likelihood_value - (slope * evaluation_index)
      split_distribution[(index_before + 1):(evaluation_index - 1)] = ((index_before +
                                                                          1):(evaluation_index - 1)) * slope + intercept
      
    } else if (((index == 1) & (length(split_distribution) > 1)) |
               ((index_difference == 1) &
                (length(split_distribution) > 1))) {
      split_distribution[evaluation_index] = likelihood_value
      
    } else if (length(split_distribution) == 1) {
      split_distribution = likelihood_value
    }
  }
  
  exp_x_values       = 10 * (0:(length(split_distribution) - 1))
  exponential_values = unlist(lapply(exp_x_values, function(x) {
    exp(-(1 / length(split_distribution)) * x)
  }))
  split_distribution[order(-split_distribution)] = exponential_values
  
  # At this point, I rescale things based on the split_selection_multiplicative_prior.
  if ((!is.null(split_selection_multiplicative_prior)) &
      (length(split_distribution) > 1)) {
    if (split_selection_multiplicative_prior < 1) {
      stop("split_selection_multiplicative_prior must be >= 1.")
    }
    possible_split_points = 1:length(split_distribution)
    multipliers           = (possible_split_points - 1) / (max(possible_split_points) -
                                                             1) * (split_selection_multiplicative_prior - 1) + 1
    split_distribution    = split_distribution * multipliers
  }
  
  return(split_distribution)
}


#' Resample the mixture components according to a merge-split algorithm followed by a
#' Gibbs sweep.
#'
#' @param my_states vector. Current regime vector.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param data_points_of_state matrix. Data instances in current state.
#' @param current_state integer. The current regime for which we are resampling components.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param verbose logical. Will print   to console if TRUE.
#'
#'
#' @noRd
#'
#'
splitmerge_gibbs_comps    = function(my_states,
                                     previous_model_fits,
                                     data_points_of_state,
                                     current_state,
                                     hyperparameters,
                                     Z_timepoint_indices,
                                     data_points_Z,
                                       verbose = FALSE) {
  # Method to perform the split-merge algorithm for the regime components, according to STORLIE
  num_of_gibbs_sweeps     = 2
  cluster_assignments     = previous_model_fits[[current_state]]$cluster_assignments
  precisions              = previous_model_fits[[current_state]]$precision
  mus                     = previous_model_fits[[current_state]]$mu
  uniform_draw_components = sample(1:length(cluster_assignments), 2, replace =
                                     F)
  n_full                  = nrow(data_points_of_state)
  assignments_maximum     = hyperparameters$component_truncation
  max_achieved            = (max(cluster_assignments) == assignments_maximum)
  move_parameter          = hyperparameters$beta
  linger_parameter        = hyperparameters$alpha
  not.cont                = hyperparameters$not.cont
  current_graph_G         = hyperparameters$G
  
  if ((cluster_assignments[uniform_draw_components[1]] == cluster_assignments[uniform_draw_components[2]]) &
      (!max_achieved)) {
    # If the uniform draws are equal, we perform a split action.
    ####################################################################################################
    ## GET THE LAUNCH STATE VECTOR
    # I will start by acquiring a launch vector for the split proposal.
    assignments_launch   = cluster_assignments
    data_point_of_first  = data_points_of_state[uniform_draw_components[1],]
    data_point_of_second = data_points_of_state[uniform_draw_components[2],]
    for (launch_index in 1:length(assignments_launch)) {
      current_data_point = data_points_of_state[launch_index,]
      first_norm_value   = norm(current_data_point - data_point_of_first, '2')
      second_norm_value  = norm(current_data_point - data_point_of_second, '2')
      if (is.na(first_norm_value) | is.na(second_norm_value)) {
        stop("Issue chosing the launch vector for component splits")
      }
      if (!(first_norm_value <= second_norm_value)) {
        # Based on the distance a given value is from the two randomly chosen data points, we'll either keep
        # a data point in its already assigned component, or put it in a completely new component.
        assignments_launch[launch_index] = max(cluster_assignments) + 1
      }
    }
    # I think that I need to redraw the parameter values here, otherwise the Gibbs step just reverts back to
    # the original vector.  I'm worried that the MH ratio will become very complicated here, though.
    previous_model_fits_temp                                      = previous_model_fits
    previous_model_fits_temp[[current_state]]$cluster_assignments = assignments_launch
    previous_model_fits_temp                                      = redraw_mixture_parameters(
      my_states,
      current_state,
      previous_model_fits_temp,
      data_points_Z,
      Z_timepoint_indices,
      linger_parameter,
      move_parameter,
      not.cont,
      current_graph_G,
      hyperparameters
    )
    
    # Now that I have a valid launch, I can perform multiple Gibbs swaps to get my
    # second-to-last proposal vector.
    first_state                 = cluster_assignments[uniform_draw_components[1]]
    second_state                = max(cluster_assignments) + 1
    indices_of_split_comp       = which(cluster_assignments == first_state)
    assignments_launch          = as.vector(
      gibbs_swap_btwn_two(
        previous_model_fits_temp[[current_state]]$precision[[first_state]],
        previous_model_fits_temp[[current_state]]$precision[[second_state]],
        previous_model_fits_temp[[current_state]]$mu[[first_state]],
        previous_model_fits_temp[[current_state]]$mu[[second_state]],
        previous_model_fits_temp[[current_state]]$component_log_probs,
        indices_of_split_comp,
        data_points_of_state,
        assignments_launch,
        first_state,
        second_state,
        2
      )[["assignments_launch"]]
    )
    ####################################################################################################
    ## GET PROPOSAL STATE VECTOR AND PROPOSAL PROBABILITY
    previous_model_fits_temp[[current_state]]$cluster_assignments = assignments_launch
    previous_model_fits_temp                                      = redraw_mixture_parameters(
      my_states,
      current_state,
      previous_model_fits_temp,
      data_points_Z,
      Z_timepoint_indices,
      linger_parameter,
      move_parameter,
      not.cont,
      current_graph_G,
      hyperparameters
    )
    # I am going to need to calculate the mismatch in the MH caused by doing this.
    total_log_prob                  = 0
    swap_return_items               = gibbs_swap_btwn_two(
      previous_model_fits_temp[[current_state]]$precision[[first_state]],
      previous_model_fits_temp[[current_state]]$precision[[second_state]],
      previous_model_fits_temp[[current_state]]$mu[[first_state]],
      previous_model_fits_temp[[current_state]]$mu[[second_state]],
      previous_model_fits_temp[[current_state]]$component_log_probs,
      indices_of_split_comp,
      data_points_of_state,
      assignments_launch,
      first_state,
      second_state,
      1
    )
    assignments_launch              = as.vector(swap_return_items[['assignments_launch']])
    total_log_prob                  = as.vector(swap_return_items[['total_log_prob']])
    
    ####################################################################################################
    ## CALCULATE THE MH STEP AND FINISH
    # Recall that this is for a split.  Gather the necessary data values.
    MH_prob                           = 0
    
    # Calculate the MH log ratio
    previous_model_fits_temp[[current_state]]$cluster_assignments = assignments_launch
    previous_model_fits_temp                                      = redraw_mixture_parameters(
      my_states,
      current_state,
      previous_model_fits_temp,
      data_points_Z,
      Z_timepoint_indices,
      linger_parameter,
      move_parameter,
      not.cont,
      current_graph_G,
      hyperparameters
    )
    marginal_term_1                 = log_Gwishart_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      current_state,
      my_states
    )
    marginal_term_2                 = log_Gwishart_marginals(
      previous_model_fits,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      current_state,
      my_states
    )
    log_leftover_value_wishartterms = marginal_term_1[["log_posterior"]] - marginal_term_2[["log_posterior"]]
    laplace_nonconvergence_flag     = max(marginal_term_1[["nonconverge_flag"]], marginal_term_2[["nonconverge_flag"]])
    
    # I need the mismatch marginals from the G-Wishart distribution, similar to what I do with the
    # merge-split algorithm on the states
    MH_prob = MH_prob + log_leftover_value_wishartterms
    
    # The total_log_prob is from merge->split, which is the FORWARD move here, so I divide by it.
    # By assumption, the merge happens with prob 1.
    MH_prob = MH_prob - total_log_prob
    log_unif = log(runif(1))
    if (is.na(MH_prob)) {
      # Catch this case first, don't accept.
    } else if (laplace_nonconvergence_flag) {
      
    } else if (MH_prob >= log_unif) {
      # Accept the new state vector.
      previous_model_fits = previous_model_fits_temp
    }
  } else {
    # If the uniform draws are not equal, we perform a merge action.
    ####################################################################################################
    ## GET THE LAUNCH STATE VECTOR
    # I will start by acquiring a launch vector for the reverse-merge proposal.
    assignments_launch          = cluster_assignments
    first_state                 = cluster_assignments[uniform_draw_components[1]]
    second_state                = cluster_assignments[uniform_draw_components[2]]
    # Pretend that I have the merged state
    assignments_launch[which(assignments_launch == second_state)] = first_state
    
    # Now, do precisely what we did for the split move, to get the launch vector.
    assignments_launch    = cluster_assignments
    data_point_of_first   = data_points_of_state[uniform_draw_components[1],]
    data_point_of_second  = data_points_of_state[uniform_draw_components[2],]
    for (launch_index in 1:length(assignments_launch)) {
      current_data_point  = data_points_of_state[launch_index,]
      first_norm_value    = norm(current_data_point - data_point_of_first, '2')
      second_norm_value   = norm(current_data_point - data_point_of_second, '2')
      if (is.na(first_norm_value) | is.na(second_norm_value)) {
        next
      } else if (!(first_norm_value <= second_norm_value)) {
        assignments_launch[launch_index] = second_state
      }
    }
    
    # I think that I need to redraw the parameter values here, otherwise the Gibbs step just reverts back to
    # the original vector.  I'm worried that the MH ratio will become very complicated here, though.
    previous_model_fits_temp                                      = previous_model_fits
    previous_model_fits_temp[[current_state]]$cluster_assignments = assignments_launch
    previous_model_fits_temp                                      = redraw_mixture_parameters(
      my_states,
      current_state,
      previous_model_fits_temp,
      data_points_Z,
      Z_timepoint_indices,
      linger_parameter,
      move_parameter,
      not.cont,
      current_graph_G,
      hyperparameters
    )
    
    # Now that I have a valid launch, I can perform multiple Gibbs swaps to get my
    # second-to-last proposal vector.
    indices_of_split_comp       = which(
      (cluster_assignments == cluster_assignments[uniform_draw_components[1]]) |
        (cluster_assignments == cluster_assignments[uniform_draw_components[2]])
    )
    assignments_launch          = as.vector(
      gibbs_swap_btwn_two(
        previous_model_fits_temp[[current_state]]$precision[[first_state]],
        previous_model_fits_temp[[current_state]]$precision[[second_state]],
        previous_model_fits_temp[[current_state]]$mu[[first_state]],
        previous_model_fits_temp[[current_state]]$mu[[second_state]],
        previous_model_fits_temp[[current_state]]$component_log_probs,
        indices_of_split_comp,
        data_points_of_state,
        assignments_launch,
        first_state,
        second_state,
        2
      )[["assignments_launch"]]
    )
    ####################################################################################################
    ## GET PROPOSAL STATE VECTOR AND PROPOSAL PROBABILITY
    # Now I have my second-to-last proposal vector.  What remains is to calculate the probability of
    # moving from it to the split vector with which I started.
    # That is, from assignments_launch -> cluster_assignments.
    swap_return_items               = gibbs_swap_btwn_two(
      previous_model_fits_temp[[current_state]]$precision[[first_state]],
      previous_model_fits_temp[[current_state]]$precision[[second_state]],
      previous_model_fits_temp[[current_state]]$mu[[first_state]],
      previous_model_fits_temp[[current_state]]$mu[[second_state]],
      previous_model_fits_temp[[current_state]]$component_log_probs,
      indices_of_split_comp,
      data_points_of_state,
      assignments_launch,
      first_state,
      second_state,
      1
    )
    assignments_launch              = as.vector(swap_return_items[['assignments_launch']])
    total_log_prob                  = as.vector(swap_return_items[['total_log_prob']])
    ####################################################################################################
    ## CALCULATE THE MH STEP AND FINISH
    # Recall that this is for a merge.  Gather the necessary data values.
    # Now that I have the proposal probability, I can reset the launch vector.
    assignments_launch          = cluster_assignments
    # Pretend that I have the merged state
    assignments_launch[which(assignments_launch == second_state)] = first_state
    
    previous_model_fits_temp[[current_state]]$cluster_assignments = assignments_launch
    previous_model_fits_temp                                      = redraw_mixture_parameters(
      my_states,
      current_state,
      previous_model_fits_temp,
      data_points_Z,
      Z_timepoint_indices,
      linger_parameter,
      move_parameter,
      not.cont,
      current_graph_G,
      hyperparameters
    )
    
    MH_prob                         = 0
    posterior_term_1                = log_Gwishart_marginals(
      previous_model_fits_temp,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      current_state,
      my_states
    )
    posterior_term_2                = log_Gwishart_marginals(
      previous_model_fits,
      Z_timepoint_indices,
      data_points_Z,
      hyperparameters,
      current_state,
      my_states
    )
    log_leftover_value_wishartterms = posterior_term_1[["log_posterior"]] - posterior_term_2[["log_posterior"]]
    laplace_nonconverge_flag        = max(posterior_term_1[["nonconverge_flag"]], posterior_term_2[["nonconverge_flag"]])
    
    
    # Calculate the MH log ratio
    MH_prob  = MH_prob + log_leftover_value_wishartterms
    
    MH_prob  = MH_prob + total_log_prob
    log_unif = log(runif(1))
    if (is.na(MH_prob)) {
      # Catch this case first, don't accept.
    } else if (laplace_nonconverge_flag) {
    } else if (MH_prob >= log_unif) {
      previous_model_fits = previous_model_fits_temp
    } else {
    }
    
  }
  
  ## PERFORM THE FINAL GIBBS SWEEP
  temp_component_sticks = rep(0, times = hyperparameters$component_truncation)
  i                     = current_state
  cluster_assignments   = previous_model_fits[[i]]$cluster_assignments
  if (anyNA(cluster_assignments)) {
    # This means that there are no data assigned to this regime.
    temp_component_sticks   = rbeta(hyperparameters$component_truncation,
                                    1,
                                    hyperparameters$dirichlet_prior)
  } else {
    for (stick_index in 1:hyperparameters$component_truncation) {
      if (stick_index %in% cluster_assignments) {
        temp_component_sticks[stick_index] = rbeta(
          1,
          1 + sum(cluster_assignments == stick_index),
          hyperparameters$dirichlet_prior + sum(cluster_assignments > stick_index)
        )
      } else {
        temp_component_sticks[stick_index] = rbeta(1, 1, hyperparameters$dirichlet_prior)
      }
    }
  }
  previous_model_fits[[i]][["component_log_probs"]] = sticks_to_log_probs(temp_component_sticks)
  previous_model_fits[[i]][["component_sticks"]]    = temp_component_sticks
  previous_model_fits[[current_state]]$cluster_assignments = as.vector(
    gibbs_swap_comps(
      data_points_of_state,
      previous_model_fits[[current_state]]$cluster_assignments,
      previous_model_fits[[current_state]]$component_log_probs,
      previous_model_fits[[current_state]]$precision,
      previous_model_fits[[current_state]]$mu,
      hyperparameters$component_truncation,
      3
    )
  )
  
  return(
    list(
      precisions  = previous_model_fits[[current_state]]$precision,
      mus         = previous_model_fits[[current_state]]$mu,
      assigns     = previous_model_fits[[current_state]]$cluster_assignments,
      comp_probs  = previous_model_fits[[current_state]]$component_log_probs,
      comp_sticks = previous_model_fits[[current_state]]$component_sticks
    )
  )
}


#' Redraw the Lambda and my parameters according to a mixture NI-G-W model according
#' to Murph et al 2023.
#'
#' @param my_states vector. Current regime vector.
#' @param state_to_redraw integer. State's parameter values to redraw.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param linger_parameter float. Prior parameter for Markov chain probability matrix.  Larger = less likely to change states.
#' @param move_parameter float. Prior parameter for Markov chain probability matrix.  Larger = more likely to change states.
#' @param not.cont vector. Indicator vector for which columns are discrete.
#' @param current_graph_G matrix. Current graph structure for Gaussian Graphical Model.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#'
#'
#'
#' @noRd
#'
#'
redraw_mixture_parameters = function(my_states,
                                     state_to_redraw,
                                     previous_model_fits,
                                     data_points_Z,
                                     Z_timepoint_indices,
                                     linger_parameter,
                                     move_parameter,
                                     not.cont,
                                     current_graph_G,
                                     hyperparameters) {
  # This is a method for redrawing the parameters, given the components within the regime.
  
  cluster_assignments       = previous_model_fits[[state_to_redraw]]$cluster_assignments
  cluster_precisions        = previous_model_fits[[state_to_redraw]]$precision
  cluster_mus               = previous_model_fits[[state_to_redraw]]$mu
  component_truncation      = hyperparameters$component_truncation
  scale_matrix              = hyperparameters$wishart_scale_matrix
  df                        = hyperparameters$wishart_df
  indices_of_state          = which(my_states == state_to_redraw)
  
  if (length(indices_of_state) > 0) {
    Z_index_of_regime_start   = Z_timepoint_indices[[min(indices_of_state)]]$timepoint_first_index
    Z_index_of_regime_end     = Z_timepoint_indices[[max(indices_of_state)]]$timepoint_last_index
    Z_values_of_regime        = data_points_Z[Z_index_of_regime_start:Z_index_of_regime_end,]
    upper_bound_is_equal      = hyperparameters$upper_bound_is_equal[Z_index_of_regime_start:Z_index_of_regime_end,]
    lower_bound_is_equal      = hyperparameters$lower_bound_is_equal[Z_index_of_regime_start:Z_index_of_regime_end,]
    is_missing                = hyperparameters$is_missing[Z_index_of_regime_start:Z_index_of_regime_end,]
  } else {
    Z_values_of_regime        = NA
  }
  
  p                           = hyperparameters$p
  
  for (component_index in 1:(component_truncation + 1)) {
    if (any(is.na(cluster_assignments)) |
        any(is.null(cluster_assignments)) |
        any(is.na(Z_values_of_regime))) {
      # In this case, there are no assignments to this component.  We just draw
      # precisions and mus from our prior.
      result              = rgwish_Rcpp(
        as.double(hyperparameters$G),
        as.double(scale_matrix),
        as.integer(df),
        as.integer(hyperparameters$p),
        as.double(1e-8)
      )
      flag                = result[['failed']]
      result              = result[['K']]
      starting_precision  = matrix(result, hyperparameters$p, hyperparameters$p)
      precision           = starting_precision
      
      new_mu              = rmvn_Rcpp(
        as.double(hyperparameters$mu_0),
        as.double((hyperparameters$lambda) *
                    starting_precision),
        as.integer(hyperparameters$p)
      )
      
      mu                  = as.matrix(new_mu, p, 1)
      
      cluster_precisions[[component_index]]                      = precision
      cluster_mus[[component_index]]                             = mu
      previous_model_fits[[state_to_redraw]]$cluster_assignments = NA
      
      
    } else if (sum(cluster_assignments == component_index) == 0) {
      # I'm doing the same thing as last time, since this condition also needs the prior draw,
      # but I can't verify it this way unless I have non-null values for cluster_assignments.
      result              = rgwish_Rcpp(
        as.double(hyperparameters$G),
        as.double(scale_matrix),
        as.integer(df),
        as.integer(hyperparameters$p),
        as.double(1e-8)
      )
      flag                = result[['failed']]
      result              = result[['K']]
      starting_precision  = matrix(result, hyperparameters$p, hyperparameters$p)
      precision           = starting_precision
      
      new_mu              = rmvn_Rcpp(
        as.double(hyperparameters$mu_0),
        as.double((hyperparameters$lambda) *
                    starting_precision),
        as.integer(hyperparameters$p)
      )
      
      mu                  = as.matrix(new_mu, p, 1)
      
      cluster_precisions[[component_index]]                      = precision
      cluster_mus[[component_index]]                             = mu
      
    } else {
      # Otherwise, draw the precision and mu from the mvn posterior.
      indices_of_comp_in_regime = which(cluster_assignments == component_index)
      
      Z_values                  = Z_values_of_regime[indices_of_comp_in_regime,]
      upper_bound_is_equal_temp = upper_bound_is_equal[indices_of_comp_in_regime,]
      lower_bound_is_equal_temp = lower_bound_is_equal[indices_of_comp_in_regime,]
      is_missing_temp           = is_missing[indices_of_comp_in_regime,]
      
      temp_n                    = nrow(Z_values)
      if (is.null(temp_n)) {
        temp_n = 1
        upper_bound_is_equal_temp = t(as.matrix(upper_bound_is_equal_temp))
        lower_bound_is_equal_temp = t(as.matrix(lower_bound_is_equal_temp))
        is_missing_temp           = t(as.matrix(is_missing_temp))
        mu                        = Z_values
        x_bar_vector              = Z_values
        nS2                       = matrix(0, hyperparameters$p, hyperparameters$p)
      } else {
        mu                        = apply(Z_values, 2, mean)
        x_bar_vector              = matrix(rep(mu, times = temp_n),
                                           temp_n ,
                                           hyperparameters$p,
                                           byrow = T)
        nS2                       = t(Z_values - x_bar_vector) %*% (Z_values -
                                                                      x_bar_vector)
      }
      
      
      new_parameter_values      = draw_mvn_parameters(
        current_graph_G,
        df,
        nS2,
        temp_n,
        mu,
        hyperparameters$mu_0,
        hyperparameters$lambda,
        scale_matrix
      )
      cluster_precisions[[component_index]] = new_parameter_values$final_K
      cluster_mus[[component_index]]        = new_parameter_values$final_mu
    }
    
  }
  
  previous_model_fits[[state_to_redraw]]$precision = cluster_precisions
  previous_model_fits[[state_to_redraw]]$mu        = cluster_mus
  return(previous_model_fits)
}


#' Redraws the Lambda and Mu parameters according to a new graph structure.
#' Details in Murph et al 2023.
#'
#' @param my_states vector. Current regime vector.
#' @param state_to_redraw integer. State for which parameters are to be redrawn.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param data_points_Z matrix.  The latent data at this MCMC iteration
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param current_graph_G matrix. Current graph structure for Gaussian Graphical Model.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param new_proposed_G matrix. New graph structure proposed for Gaussian Graphical Model.
#' @param selected_edge_i integer. Row index of changed element in G.
#' @param selected_edge_j integer. Col index of changed element in G.
#' @param g_sampling_distribution matrix. Probability distribution of edge inclusion.
#'
#'
#'
#' @noRd
#'
#'
redraw_G_with_mixture     = function(my_states,
                                     state_to_redraw,
                                     previous_model_fits,
                                     data_points_Z,
                                     Z_timepoint_indices,
                                     current_graph_G,
                                     hyperparameters,
                                     new_proposed_G,
                                     selected_edge_i,
                                     selected_edge_j,
                                     g_sampling_distribution) {
  upper_bound_is_equal <- lower_bound_is_equal <- is_missing <- NULL
  n.cores_eachchild         = 1
  new_G_proposed            = new_proposed_G
  cluster_assignments       = previous_model_fits[[state_to_redraw]]$cluster_assignments
  cluster_precisions        = previous_model_fits[[state_to_redraw]]$precision
  cluster_mus               = previous_model_fits[[state_to_redraw]]$mu
  component_truncation      = hyperparameters$component_truncation
  p                         = hyperparameters$p
  indices_of_state          = which(my_states == state_to_redraw)
  spread_parameter_sd2      = 0.5
  g.prior                   = g_sampling_distribution[selected_edge_i, selected_edge_j]
  
  scale_matrix          = hyperparameters$wishart_scale_matrix
  df                    = hyperparameters$wishart_df
  mean_hyperparameter   = hyperparameters$mu_0
  lambda_hyperparameter = hyperparameters$lambda
  total_log_MH          = 0
  new_scale_mat_proposal = scale_matrix
  
  if (length(indices_of_state) > 0) {
    Z_index_of_regime_start   = Z_timepoint_indices[[min(indices_of_state)]]$timepoint_first_index
    Z_index_of_regime_end     = Z_timepoint_indices[[max(indices_of_state)]]$timepoint_last_index
    Z_values_of_regime        = data_points_Z[Z_index_of_regime_start:Z_index_of_regime_end,]
  } else {
    Z_values_of_regime        = NA
  }
  
  for (component_index in 1:component_truncation) {
    if (any(is.na(cluster_assignments)) |
        any(is.null(cluster_assignments)) |
        any(is.na(Z_values_of_regime))) {
      # In this case, there are no assignments to this component.  We just draw
      # precisions and mus from our prior.  This may be cheating, but it also seems practical.
      result               = rgwish_Rcpp(
        as.double(new_proposed_G),
        as.double(new_scale_mat_proposal),
        as.integer(df),
        as.integer(hyperparameters$p),
        as.double(1e-8)
      )
      flag                = result[['failed']]
      result              = result[['K']]
      starting_precision  = matrix(result, hyperparameters$p, hyperparameters$p)
      precision           = starting_precision
      
      new_mu              = rmvn_Rcpp(
        as.double(hyperparameters$mu_0),
        as.double((hyperparameters$lambda) *
                    starting_precision),
        as.integer(hyperparameters$p)
      )
      
      mu                  = as.matrix(new_mu, p, 1)
      
      cluster_precisions[[component_index]]                      = precision
      cluster_mus[[component_index]]                             = mu
      previous_model_fits[[state_to_redraw]]$cluster_assignments = NA
      
    } else if (sum(cluster_assignments == component_index) == 0) {
      result              = rgwish_Rcpp(
        as.double(new_proposed_G),
        as.double(new_scale_mat_proposal),
        as.integer(df),
        as.integer(hyperparameters$p),
        as.double(1e-8)
      )
      flag                = result[['failed']]
      result              = result[['K']]
      starting_precision  = matrix(result, hyperparameters$p, hyperparameters$p)
      precision           = starting_precision
      
      new_mu              = rmvn_Rcpp(
        as.double(hyperparameters$mu_0),
        as.double((hyperparameters$lambda) *
                    starting_precision),
        as.integer(hyperparameters$p)
      )
      
      mu                  = as.matrix(new_mu, p, 1)
      
      cluster_precisions[[component_index]]                      = precision
      cluster_mus[[component_index]]                             = mu
      
    } else {
      # Otherwise, draw the precision and mu from the mvn posterior.
      indices_of_comp_in_regime = which(cluster_assignments == component_index)
      
      Z_values                  = Z_values_of_regime[indices_of_comp_in_regime,]
      upper_bound_is_equal_temp = upper_bound_is_equal[indices_of_comp_in_regime,]
      lower_bound_is_equal_temp = lower_bound_is_equal[indices_of_comp_in_regime,]
      is_missing_temp           = is_missing[indices_of_comp_in_regime,]
      
      temp_n                    = nrow(Z_values)
      if (is.null(temp_n)) {
        temp_n                    = 1
        upper_bound_is_equal_temp = t(as.matrix(upper_bound_is_equal_temp))
        lower_bound_is_equal_temp = t(as.matrix(lower_bound_is_equal_temp))
        is_missing_temp           = t(as.matrix(is_missing_temp))
        mu                        = Z_values
        x_bar_vector              = Z_values
        nS2                       = matrix(0, hyperparameters$p, hyperparameters$p)
      } else {
        mu                        = apply(Z_values, 2, mean)
        x_bar_vector              = matrix(rep(mu, times = temp_n),
                                           temp_n ,
                                           hyperparameters$p,
                                           byrow = T)
        nS2                       = t(Z_values - x_bar_vector) %*% (Z_values -
                                                                      x_bar_vector)
      }
      
      # Sample the half-jump lambda_0 (I do it out here b/c of a FORTRAN conflict within Armadillo code)
      lambda_0       = rgwish_Rcpp(
        as.double(new_G_proposed),
        as.double(new_scale_mat_proposal),
        as.integer(df),
        as.integer(p),
        as.double(1e-8)
      )
      flag           = lambda_0[['failed']]
      lambda_0       = lambda_0[['K']]
      lambda_0       = matrix(lambda_0, p, p)
      current_lambda = cluster_precisions[[component_index]]
      current_mu     = cluster_mus[[component_index]]
      edge_updated_i = selected_edge_i - 1
      edge_updated_j = selected_edge_j - 1
      current_G      = current_graph_G
      
      output = DRJ_MCMC_singlestep(
        current_lambda,
        lambda_0,
        current_G,
        p,
        n.cores_eachchild,
        edge_updated_i,
        edge_updated_j,
        new_scale_mat_proposal,
        temp_n,
        mu,
        nS2,
        df,
        spread_parameter_sd2,
        mean_hyperparameter,
        lambda_hyperparameter,
        g.prior
      )
      
      total_log_MH                          = total_log_MH + output$log_MH_ratio
      cluster_precisions[[component_index]] = output$new_lambda
      cluster_mus[[component_index]]        = output$new_mu
    }
    
  }
  
  previous_model_fits[[state_to_redraw]]$precision = cluster_precisions
  previous_model_fits[[state_to_redraw]]$mu        = cluster_mus
  return(list(
    previous_model_fits_item = previous_model_fits,
    log_MH_ratio            = total_log_MH
  ))
}


#' Get log density of data according to the likelihood, assuming a mixture distribution.
#'
#' @param current_state integer. State for which we need to calculate a mixture log density.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param data_points_of_state matrix. Data instances in current state.
#'
#'
#'
#' @noRd
#'
#'
get_mixture_log_density   = function(current_state,
                                     previous_model_fits,
                                     data_points_of_state) {
  # This method calculates the log normal density value, given the component assignments.
  current_density_value   = 0
  cluster_assignments     = previous_model_fits[[current_state]]$cluster_assignments
  cluster_precisions      = previous_model_fits[[current_state]]$precision
  cluster_mus             = previous_model_fits[[current_state]]$mu
  cluster_probs           = previous_model_fits[[current_state]]$cluster_probs
  all_indices             = 1:nrow(data_points_of_state)
  for (cluster_index in 1:max(cluster_assignments)) {
    current_precision     = cluster_precisions[[cluster_index]]
    current_mu            = cluster_mus[[cluster_index]]
    indicies_of_cluster   = which(cluster_assignments == cluster_index)
    data_of_cluster       = data_points_of_state[indicies_of_cluster,]
    if (length(indicies_of_cluster) == 1) {
      data_of_cluster = t(as.matrix(data_of_cluster))
    }
    if (is.null(data_of_cluster)) {
      stop("no data in the cluster")
    }
    if (is.null(current_mu)) {
      stop("no mean vector")
    }
    current_density_value = current_density_value + log_dmvnrm_arma_regular(data_of_cluster, current_mu, current_precision)
  }
  return(current_density_value)
}


#' Gets the log mixture likelihood of a single data instance.
#'
#' @param index_of_observation integer. Index of full data corresponding to single data instance.
#' @param my_states vector. Current regime vector.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param data_of_observation vector. Data of a single observation.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#'
#'
#' @noRd
#'
#'
get_mix_log_dens_at_obs   = function(index_of_observation,
                                     my_states,
                                     previous_model_fits,
                                     data_of_observation,
                                     Z_timepoint_indices) {
  # This method calculates the log normal density value, given the component assignments, for a single data point, NOT for
  # the whole state.
  first_index             = Z_timepoint_indices[[index_of_observation]]$timepoint_first_index
  last_index              = Z_timepoint_indices[[index_of_observation]]$timepoint_last_index
  first_index_of_state    = Z_timepoint_indices[[min(which(my_states == my_states[index_of_observation]))]]$timepoint_first_index
  current_density_value   = 0
  cluster_assignments     = previous_model_fits[[my_states[index_of_observation]]]$cluster_assignments
  cluster_assignments     = cluster_assignments[(first_index - first_index_of_state +
                                                   1):(last_index - first_index_of_state + 1)]
  cluster_precisions      = previous_model_fits[[my_states[index_of_observation]]]$precision
  cluster_mus             = previous_model_fits[[my_states[index_of_observation]]]$mu
  
  for (cluster_index in unique(cluster_assignments)) {
    current_precision     = cluster_precisions[[cluster_index]]
    current_mu            = cluster_mus[[cluster_index]]
    cluster_indices       = which(cluster_assignments == cluster_index)
    data_of_cluster       = data_of_observation[cluster_indices,]
    if (length(cluster_indices) == 1) {
      data_of_cluster     = t(as.matrix(data_of_cluster))
    }
    if (is.null(data_of_cluster)) {
      stop("no data in the cluster (at obs)")
    }
    if (is.null(current_mu)) {
      stop("no mu vector")
    }
    if (is.nan(log_dmvnrm_arma_regular(data_of_cluster, current_mu, current_precision))) {
      stop("issue with log normal density evaluation")
    }
    current_density_value = current_density_value + log_dmvnrm_arma_regular(data_of_cluster, current_mu, current_precision)
  }
  return(current_density_value)
}


#' Method to redraw estimated latent data for a single state.
#'
#' @param state_to_redraw integer. State whose latent data we will redraw in this method.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param latent_data_of_state matrix. The latent data at this MCMC iteration corresponding to given state.
#' @param raw_data_from_state matrix. The raw data at this MCMC iteration corresponding to given state.
#' @param is_missing_state matrix. Indicators as to when elements of raw_data_from_state are missing.
#' @param upper_bound_is_equal_state matrix. Indicators as to when elements of raw_data_from_state achieve their upper bounds.
#' @param lower_bound_is_equal_state matrix. Indicators as to when elements of raw_data_from_state achieve their lower bounds.
#' @param n_cores_each_child integer. Number of cores to use for paralellization at the C/C++ level.
#'
#' @noRd
#'
#'
redraw_latent_data        = function(state_to_redraw,
                                     previous_model_fits,
                                     hyperparameters,
                                     latent_data_of_state,
                                     raw_data_from_state,
                                     is_missing_state,
                                     upper_bound_is_equal_state,
                                     lower_bound_is_equal_state,
                                     ordinal_levels,
                                     levels_assignments,
                                     discrete_levels_indicator,
                                     n_cores_each_child = 1) {
  # Method to redraw latent data FOR A SINGLE STATE ONLY.  Assumes Normal Mixture Model.
  components        = previous_model_fits[[state_to_redraw]]$cluster_assignments
  precisions        = previous_model_fits[[state_to_redraw]]$precision
  mus               = previous_model_fits[[state_to_redraw]]$mu
  continuous        = as.numeric(!hyperparameters$not.cont)
  
  for (component_index in 1:max(components)) {
    indices_of_data           = which(components == component_index)
    temp_data                 = latent_data_of_state[indices_of_data, ]
    temp_data_raw             = raw_data_from_state[indices_of_data, ]
    n_value                   = nrow(temp_data)
    
    if (is.null(n_value)) {
      n_value = 1
      temp_data = t(as.matrix(temp_data))
      temp_data_raw = t(as.matrix(temp_data_raw))
      upper_bound_is_equal_temp = t(as.matrix(upper_bound_is_equal_state[indices_of_data,]))
      lower_bound_is_equal_temp = t(as.matrix(lower_bound_is_equal_state[indices_of_data,]))
      is_missing_temp           = t(as.matrix(is_missing_state[indices_of_data,]))
    } else {
      upper_bound_is_equal_temp = upper_bound_is_equal_state[indices_of_data,]
      lower_bound_is_equal_temp = lower_bound_is_equal_state[indices_of_data,]
      is_missing_temp           = is_missing_state[indices_of_data,]
    }
    prec_of_component         = precisions[[component_index]]
    mu_of_component           = mus[[component_index]]
    
    redraw_output             = redraw_Z_arma(
      temp_data,
      prec_of_component,
      mu_of_component,
      hyperparameters$p,
      hyperparameters$lower_bounds,
      hyperparameters$upper_bounds,
      lower_bound_is_equal_temp,
      upper_bound_is_equal_temp,
      is_missing_temp,
      continuous,
      temp_data_raw,
      ordinal_levels,
      levels_assignments,
      discrete_levels_indicator,
      n_cores_each_child
    )
    new_data_for_comp         = matrix(redraw_output, n_value, hyperparameters$p)
    
    # Update the data for that component.
    latent_data_of_state[indices_of_data,] = new_data_for_comp
  }
  return(latent_data_of_state)
}


#' Swap components from one regime to another whenever Gibbs sampling reassigns it.
#'
#' @param state_value_gained integer. State that is gaining additional components.
#' @param state_value_lost integer. State that is loosing additional components.
#' @param value_index integer. State value that is changing its regime.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#'
#'
#' @noRd
#'
#'
update_regime_components  = function(state_value_gained,
                                     state_value_lost,
                                     value_index,
                                     previous_model_fits,
                                     Z_timepoint_indices,
                                     hyperparameters) {
  ## I need to be able to swap one time-point-instance's component data from one regime to another.
  components_gained_state     = previous_model_fits[[state_value_gained]]$cluster_assignments
  precisions_gained_state     = previous_model_fits[[state_value_gained]]$precision
  mus_gained_state            = previous_model_fits[[state_value_gained]]$mu
  comp_log_probs_gained_state = previous_model_fits[[state_value_gained]]$component_log_probs
  components_lost_state       = previous_model_fits[[state_value_lost]]$cluster_assignments
  precisions_lost_state       = previous_model_fits[[state_value_lost]]$precision
  mus_lost_state              = previous_model_fits[[state_value_lost]]$mu
  comp_log_probs_lost_state   = previous_model_fits[[state_value_lost]]$component_log_probs
  index_first                 = Z_timepoint_indices[[value_index]]$timepoint_first_index
  index_last                  = Z_timepoint_indices[[value_index]]$timepoint_last_index
  
  if (state_value_gained < state_value_lost) {
    # Add all the data for the timepoint AFTER the existing data values
    # Take the VALUE from the beginning of the state_value_lost and add to the
    # end of state_value_gained.
    indices_of_component      = 1:(index_last - index_first + 1)
    components_to_add         = components_lost_state[indices_of_component]
    new_component_assignments = components_to_add
    components_lost_state     = components_lost_state[-indices_of_component]
    available_components      = 1:hyperparameters$component_truncation
    log_prob_forward          = 0
    log_prob_backwards        = 0
    
    # Drop from the 'available_components' any places where the prob is way too small.
    prob_too_small       = log(1e-8)
    available_components = available_components[comp_log_probs_gained_state >
                                                  prob_too_small]
    for (comp_index in sort(unique(components_to_add))) {
      if (length(available_components) == 0) {
        new_component_assignments[which(components_to_add == comp_index)] = 1
        next
      }
      prob_distn              = rep(0, times = length(available_components))
      log_prob_backwards      =  log_prob_backwards + comp_log_probs_lost_state[comp_index] *
        sum(components_to_add == comp_index)
      for (new_comp_index in 1:length(available_components)) {
        prob_distn[new_comp_index] = comp_log_probs_gained_state[available_components[new_comp_index]] *
          (sum(components_to_add == comp_index))
      }
      # Determine to which component I want to assign this group of new values.
      sum_density_values     = 0
      if (is.na(sum(prob_distn)) | (sum(prob_distn) == 0)) {
        log_prob_forward       = -Inf
      } else {
        suppressWarnings({
          rand_value = runif(1, min = 0, max = sum(prob_distn))
        })
        
        for (split_index in 1:length(available_components)) {
          sum_density_values = sum_density_values + exp(prob_distn[split_index])
          if (is.na(rand_value)) {
            log_prob_forward       = -Inf
            break
          } else if (rand_value <= sum_density_values) {
            log_prob_forward                                                = log_prob_forward + prob_distn[split_index]
            new_component_assignments[which(components_to_add == comp_index)] = available_components[split_index]
            available_components                                            = available_components[-split_index]
            break
          }
        }
      }
      
    }
    
    # Otherwise, there is no need to touch the MVN parameters.
    components_gained_state = c(components_gained_state, new_component_assignments)
    
  } else {
    # Add all the data for the timepoint BEFORE the existing data values.
    # Otherwise, there is no need to touch the MVN parameters.
    indices_of_component      = (length(components_lost_state) - (index_last -
                                                                    index_first + 1) + 1):length(components_lost_state)
    components_to_add         = components_lost_state[indices_of_component]
    new_component_assignments = components_to_add
    components_lost_state     = components_lost_state[-indices_of_component]
    available_components      = 1:hyperparameters$component_truncation
    log_prob_forward          = 0
    log_prob_backwards        = 0
    
    # Drop from the 'available_components' any places where the prob is way too small.
    prob_too_small = 5e-5
    available_components = available_components[comp_log_probs_gained_state >
                                                  log(prob_too_small)]
    for (comp_index in sort(unique(components_to_add))) {
      if (length(available_components) == 0) {
        new_component_assignments[which(components_to_add == comp_index)] = 1
        next
      }
      log_prob_backwards           =  log_prob_backwards + comp_log_probs_lost_state[comp_index] *
        sum(components_to_add == comp_index)
      prob_distn                   = rep(0, times = length(available_components))
      for (new_comp_index in 1:length(available_components)) {
        prob_distn[new_comp_index] = comp_log_probs_gained_state[available_components[new_comp_index]] *
          (sum(components_to_add == comp_index))
      }
      # Determine to which component I want to assign this group of new values.
      sum_density_values     = 0
      if (is.na(sum(prob_distn)) | (sum(prob_distn) == 0)) {
        log_prob_forward       = -Inf
      } else {
        suppressWarnings({
          rand_value = runif(1, min = 0, max = sum(prob_distn))
        })
        
        for (split_index in 1:length(available_components)) {
          sum_density_values = sum_density_values + exp(prob_distn[split_index])
          if (is.na(rand_value)) {
            log_prob_forward       = -Inf
            break
          } else if (rand_value <= sum_density_values) {
            log_prob_forward                                                = log_prob_forward + prob_distn[split_index]
            new_component_assignments[which(components_to_add == comp_index)] = available_components[split_index]
            available_components                                            = available_components[-split_index]
            break
          }
        }
      }
    }
    # Otherwise, there is no need to touch the MVN parameters.
    components_gained_state = c(new_component_assignments, components_gained_state)
    
  }
  previous_model_fits[[state_value_gained]]$cluster_assignments = components_gained_state
  previous_model_fits[[state_value_lost]]$cluster_assignments   = components_lost_state
  
  # I believe that the update will go away in the MH ratio.  The proposal
  # distribution should be the same as the likelihood evaluation.
  return(list(previous_model_fits_item = previous_model_fits))#,
  
}


#' Splits components of a state following a regime split in the mergesplit algorithm.
#'
#' @param state_split integer. State that is being split.
#' @param my_states vector. Current regime vector.
#' @param index_of_split integer. Data index where the state split occurs.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param Z_timepoint_indices rlist. Indices of data_points_Z corresponded to different days.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param data_of_state_changed matrix. Latent data of the state that is changing regimes.
#'
#'
#' @noRd
#'
#'
split_regime_components   = function(state_split,
                                     my_states,
                                     index_of_split,
                                     previous_model_fits,
                                     Z_timepoint_indices,
                                     hyperparameters,
                                     data_of_state_changed) {
  ## I need to be able to swap one time-point-instance's component data from one regime to another.
  ## Note that here the VALUE_INDEX is an element of the regime that happens earlier in time.
  components_lost_state     = previous_model_fits[[state_split]]$cluster_assignments
  precisions_gained_state   = previous_model_fits[[max(my_states) + 1]]$precision
  mus_gained_state          = previous_model_fits[[max(my_states) + 1]]$mu
  which_states              = which(my_states == state_split)
  index_first               = Z_timepoint_indices[[min(which_states)]]$timepoint_first_index
  index_last                = Z_timepoint_indices[[index_of_split]]$timepoint_last_index
  indices_of_component      = (index_last - index_first + 2):length(components_lost_state)
  components_gained_state   = components_lost_state[indices_of_component]
  components_to_add         = components_gained_state
  new_component_assignments = components_to_add
  components_lost_state     = components_lost_state[-indices_of_component]
  
  comp_log_probs_gained_state         = previous_model_fits[[max(my_states) +
                                                               1]]$component_log_probs
  comp_log_probs_lost_state           = previous_model_fits[[state_split]]$component_log_probs
  
  available_components    = 1:hyperparameters$component_truncation
  log_prob_forward        = 0
  
  # Drop from the 'available_components' any places where the prob is way too small.
  unique_components = sort(unique(components_to_add))
  f = function(x) {
    sum(components_to_add == x)
  }
  counts_of_each    = sapply(unique_components, f)
  comps_to_iterate  = unique_components[order(counts_of_each, decreasing = TRUE)]
  
  for (comp_index in comps_to_iterate) {
    if (length(available_components) == 0) {
      new_component_assignments[which(components_to_add == comp_index)] = 1
      next
    }
    
    prob_distn              = rep(0, times = length(available_components))
    prob_distn_just_density = rep(0, times = length(available_components))
    for (new_comp_index in 1:length(available_components)) {
      indices_of_this_component = which(components_to_add == comp_index)
      if (length(indices_of_this_component) > 1) {
        data_of_this_component                  = data_of_state_changed[indices_of_this_component,]
        current_mu                              = mus_gained_state[[available_components[new_comp_index]]]
        current_precision                       = precisions_gained_state[[available_components[new_comp_index]]]
        prob_distn[new_comp_index]              = log_dmvnrm_arma_regular(data_of_this_component,
                                                                          current_mu,
                                                                          current_precision)
        prob_distn_just_density[new_comp_index] = prob_distn[new_comp_index]
      } else if (length(indices_of_this_component) == 1) {
        data_of_this_component                  = matrix(data_of_state_changed[indices_of_this_component,], nrow =
                                                           1)
        current_mu                              = mus_gained_state[[available_components[new_comp_index]]]
        current_precision                       = precisions_gained_state[[available_components[new_comp_index]]]
        prob_distn[new_comp_index]              = log_dmvnrm_arma_regular(data_of_this_component,
                                                                          current_mu,
                                                                          current_precision)
        prob_distn_just_density[new_comp_index] = prob_distn[new_comp_index]
      }
      
      prob_distn[new_comp_index] = prob_distn[new_comp_index] +
        comp_log_probs_gained_state[available_components[new_comp_index]] *
        (sum(components_to_add == comp_index))
    }
    orig_prob_distn = prob_distn
    prob_distn      = prob_distn - max(prob_distn)
    prob_distn      = prob_distn - log(sum(exp(prob_distn)))
    
    # Determine to which component I want to assign this group of new values.
    sum_density_values     = 0
    suppressWarnings({
      rand_value = runif(1, min = 0, max = sum(exp(prob_distn)))
    })
    
    for (split_index in 1:length(available_components)) {
      sum_density_values = sum_density_values + exp(prob_distn[split_index])
      if (is.na(rand_value)) {
        log_prob_forward = Inf
        break
      } else if (rand_value <= sum_density_values) {
        log_prob_forward       = log_prob_forward + prob_distn[split_index]
        break
      }
    }
    new_component_assignments[which(components_to_add == comp_index)] = available_components[split_index]
    available_components                                            = available_components[-split_index]
  }
  
  previous_model_fits[[max(my_states) + 1]]$cluster_assignments = new_component_assignments
  previous_model_fits[[state_split]]$cluster_assignments      = components_lost_state
  return(
    list(
      previous_model_fits_item                = previous_model_fits,
      log_forward_prob_item                   = log_prob_forward,
      original_merged_assignments_state2_item = components_gained_state
    )
  )
}


#' Adding components to a state whenever a merge occurs in the merge-split algorithm for regimes.
#'
#' @param state_merged integer. Index of regime that is being merged with the regime above.
#' @param previous_model_fits rlist. The parameter fit at this MCMC iteration.
#' @param hyperparameters rlist. Various hyperparameters according to Bayesian setup.
#' @param data_in_second_state matrix. Latent data of the state that is joining state_merged.
#'
#'
#' @noRd
#'
#'
merge_regime_components   = function(state_merged,
                                     previous_model_fits,
                                     hyperparameters,
                                     data_in_second_state) {
  ## I need to be able to swap one time-point-instance's component data from one regime to another.
  ## Note that here the VALUE_INDEX is an element of the regime that happens earlier in time.
  components_lost_state     = previous_model_fits[[state_merged + 1]]$cluster_assignments
  precisions_lost_state     = previous_model_fits[[state_merged + 1]]$precision
  mus_lost_state            = previous_model_fits[[state_merged + 1]]$mu
  components_gained_state   = previous_model_fits[[state_merged]]$cluster_assignments
  precisions_gained_state   = previous_model_fits[[state_merged]]$precision
  mus_gained_state          = previous_model_fits[[state_merged]]$mu
  
  comp_log_probs_gained_state   = previous_model_fits[[state_merged]]$component_log_probs
  comp_log_probs_lost_state     = previous_model_fits[[state_merged + 1]]$component_log_probs
  
  components_to_add         = components_lost_state
  new_component_assignments = components_to_add
  
  available_components    = 1:hyperparameters$component_truncation
  log_prob_forward        = 0
  
  unique_components = sort(unique(components_to_add))
  f = function(x) {
    sum(components_to_add == x)
  }
  counts_of_each    = sapply(unique_components, f)
  comps_to_iterate  = unique_components[order(counts_of_each, decreasing = TRUE)]
  
  for (comp_index in comps_to_iterate) {
    if (length(available_components) == 0) {
      new_component_assignments[which(components_to_add == comp_index)] = 1
    }
    
    prob_distn              = rep(0, times = length(available_components))
    prob_distn_just_density = rep(0, times = length(available_components))
    for (new_comp_index in 1:length(available_components)) {
      # Add in log prob mass from the density:
      indices_of_this_component = which(components_to_add == comp_index)
      if (length(indices_of_this_component) > 1) {
        data_of_this_component                  = data_in_second_state[indices_of_this_component,]
        current_mu                              = mus_gained_state[[available_components[new_comp_index]]]
        current_precision                       = precisions_gained_state[[available_components[new_comp_index]]]
        prob_distn[new_comp_index]              = log_dmvnrm_arma_regular(data_of_this_component,
                                                                          current_mu,
                                                                          current_precision)
        prob_distn_just_density[new_comp_index] = prob_distn[new_comp_index]
      } else if (length(indices_of_this_component) == 1) {
        data_of_this_component                  = matrix(data_in_second_state[indices_of_this_component,], nrow =
                                                           1)
        current_mu                              = mus_gained_state[[available_components[new_comp_index]]]
        current_precision                       = precisions_gained_state[[available_components[new_comp_index]]]
        prob_distn[new_comp_index]              = log_dmvnrm_arma_regular(data_of_this_component,
                                                                          current_mu,
                                                                          current_precision)
        prob_distn_just_density[new_comp_index] = prob_distn[new_comp_index]
      }
      
      # Add in the log prob mass from the component probabilities:
      prob_distn[new_comp_index] = prob_distn[new_comp_index] +
        comp_log_probs_gained_state[available_components[new_comp_index]] *
        (sum(components_to_add == comp_index))
    }
    orig_prob_distn = prob_distn
    prob_distn      = prob_distn - max(prob_distn)
    prob_distn      = prob_distn - log(sum(exp(prob_distn)))
    
    # Determine to which component I want to assign this group of new values.
    sum_density_values     = 0
    suppressWarnings({
      rand_value = runif(1, min = 0, max = sum(exp(prob_distn)))
    })
    for (split_index in 1:length(available_components)) {
      sum_density_values = sum_density_values + exp(prob_distn[split_index])
      if (is.na(rand_value)) {
        log_prob_forward = Inf
        break
      } else if (rand_value <= sum_density_values) {
        log_prob_forward       = log_prob_forward + prob_distn[split_index]
        break
      }
    }
    new_component_assignments[which(components_to_add == comp_index)] = available_components[split_index]
    available_components                                              = available_components[-split_index]
  }
  
  previous_model_fits[[state_merged]]$cluster_assignments     = c(components_gained_state, new_component_assignments)
  previous_model_fits[[state_merged + 1]]$cluster_assignments = NA
  
  return(
    list(
      previous_model_fits_item              = previous_model_fits,
      log_forward_prob_item                 = log_prob_forward,
      original_split_components_state2_item = components_lost_state
    )
  )
}


#' Converts the sticks in a stick-breaking process into component inclusion probabilities.
#'
#' @param component_sticks vector. Sticks according to a stick-breaking process.  Dirichlet process prior.
#'
#'
#' @noRd
#'
#'
sticks_to_log_probs       = function(component_sticks) {
  n               = length(component_sticks)
  component_probs = rep(0, times = n)
  
  for (stick_index in 1:n) {
    component_probs[stick_index] = log(component_sticks[stick_index])
    for (backward_index in 1:(stick_index - 1)) {
      if (stick_index != 1) {
        component_probs[stick_index] = component_probs[stick_index] + log(1 - component_sticks[backward_index])
      }
    }
  }
  return(component_probs)
}


#' Gets log probability of regime vector according to the Markov process.
#'
#' @param state_vector vector. Current estimate regimes.
#' @param transition_probabilities vector. Estimate for transition probabilities according to the Markov Chain.
#'
#'
#' @noRd
#'
#'
calc_regimes_log_prob     = function(state_vector, transition_probabilities) {
  log_prob = 0
  for (regime_index in 1:(length(state_vector) - 1)) {
    if (state_vector[regime_index] != state_vector[regime_index + 1]) {
      log_prob = log_prob + log(1 - transition_probabilities[state_vector[regime_index]])
    } else {
      log_prob = log_prob + log(transition_probabilities[state_vector[regime_index]])
    }
  }
  return(log_prob)
}
