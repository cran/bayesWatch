# Script file for fitting a Hidden Markov Model with Dirichlet Prior
## Author: Alexander Murph
require(parallel)
require(Rcpp)
require(Matrix)
require(CholWishart)
require(ggplot2)
require(Hotelling)
library(ess)


#' Print function for a bayesWatch object.  Prints only the posterior change-point probabilities.
#'
#'
#' @param x bayesWatch object. Fit from bayesWatch main method.
#' @param ... Additional plotting arguments.
#'
#' @usage \method{print}{bayesWatch}(x) <- NULL
#' @exportS3Method
print.bayesWatch = function(x, ...)
{
  cat("\n     bayesWatch object\n")
  cat("----------------------------------\n")
  print(x$changepoint_probabilities)
}


#' Fit a bayesWatch object.
#'
#' @description Main method of package.  MCMC sampling for change-point probabilities with fault detection
#' according to the model by Murph et al. 2023.  Creates a bayesWatch object for analysis of change-points.
#'
#' @param data_woTimeValues matrix. Raw data matrix without datetime stamps.
#' @param time_of_observations vector. Datetime stamps for every data instance in data_woTimeValues.
#' @param time_points vector. Time points that mark each 'day' of time. Range should include every datetime in time_of_observations.
#' @param not.cont vector. Indicator variable as to which columns are discrete.
#' @param iterations integer. Number of MCMC samples to take (including burn-in).
#' @param burnin integer. Number of burn-in samples. iterations > burnin necessarily.
#' @param lower_bounds vector. Lower bounds for each data column.
#' @param upper_bounds vector. Upper bounds for each data column.
#' @param ordinal_indicators vector. Discrete values, one for each column, indicating which variables are ordinal.
#' @param list_of_ordinal_levels vector. Discrete values, one for each column, indicating which variables are part of the same ordinal group.
#' @param categorical_indicators vector.  Each nominal d categorical variable must be broken down into d different indicator variables.  This vector marks which variables are such indicators.
#' @param previous_states vector. Starting regime vector, if known, of the same length as the number of 'days' in time_points.
#' @param previous_model_fits rlist. Starting parameter fits corresponding to regime vector previous_states.
#' @param linger_parameter float. Prior parameter for Markov chain probability matrix.  Larger = less likely to change states.
#' @param move_parameter float. Prior parameter for Markov chain probability matrix.  Larger = more likely to change states.
#' @param g.prior float in (0,1). Prior probability on edge inclusion for graph structure G.
#' @param set_G matrix. Starting graph structure, if known.
#' @param wishart_df_initial integer (>= 3).  Starting DF for G-Wishart prior.
#' @param lambda float. Parameter for NI-G-W prior, controls affect of precision sample on the center sample.
#' @param g_sampling_distribution matrix. Prior probability on edge inclusion if not uniform across G.
#' @param n.cores integer. Number of cores available for parallelization.
#' @param scaleMatrix matrix. Parameter for NI-G-W prior.
#' @param allow_for_mixture_models logical. Whether or not method should fix mixture distributions to regimes.
#' @param dirichlet_prior float. Parameter for the dirichlet process for fitting components in the mixture model.
#' @param component_truncation integer. Maximum component allowed.  Should be sufficiently large.
#' @param regime_truncation integer. Maximum regime allowed. Should be sufficiently large.
#' @param hyperprior_b integer. Hyperprior on Wishart distribution fit to the scaleMatrix.
#' @param model_params_save_every integer. How frequently to save model fits for the fault detection method.
#' @param simulation_iter integer. Used for simulation studies.  Deprecated value at package launch.
#' @param T2_window_size integer. Length of sliding window for Hotelling T2 pre-step.  Used when an initial value for previous_states is not provided.
#' @param determining_p_cutoff logical. Method for estimating the probability cutoff on the posterior distribution for determining change-points.  Deprecated at package launch date.
#' @param variable_names vector. Vector of names of columnsof data_woTimeValues.
#' @param prob_cutoff float. Changepoints are determined (for fault detection process) if posterior probability exceeds this value.
#' @param model_log_type character vector. The type of log (used to distinguish logfiles).
#' @param regime_selection_multiplicative_prior float. Must be >=1.  Gives additional probability to the most recent day for the selection of a new split point.
#' @param split_selection_multiplicative_prior float. 
#' @param is_initial_fit logical. True when there is no previously fit bayesWatch object fed through the algorithm..
#' @param verbose logical. Prints verbose model output for debugging when TRUE.  It is highly recommended that you pipe this to a text file.
#'
#' @return bayesWatch object. A model fit for the analysis of posterior change-points and fault detection.
#' @import parallel
#' @export
#'
#' @examples
#' \donttest{
#' library(bayesWatch)
#' data("full_data")
#' data("day_of_observations")
#' data("day_dts")
#'
#' x       = bayeswatch(full_data, day_of_observations, day_dts,
#'                    iterations = 500, g.prior = 1, linger_parameter = 20, n.cores=3,
#'                    wishart_df_initial = 3, hyperprior_b = 3, lambda = 5)
#'
#' print(x)
#' plot(x)
#' detect_faults(x)
#' }
bayeswatch = function(data_woTimeValues,
                      time_of_observations,
                      time_points,
                      variable_names = 1:ncol(data_woTimeValues),
                      not.cont = NULL,
                      iterations = 100,
                      burnin = floor(iterations / 2),
                      lower_bounds = NULL,
                      upper_bounds = NULL,
                      ordinal_indicators = NULL,
                      list_of_ordinal_levels = NULL,
                      categorical_indicators = NULL,
                      previous_states = NULL,
                      previous_model_fits = NULL,
                      linger_parameter = 500,
                      move_parameter = 100,
                      g.prior = 0.2,
                      set_G = NULL,
                      wishart_df_initial = 1500,
                      lambda = 1500,
                      g_sampling_distribution = NULL,
                      n.cores = 1,
                      scaleMatrix = NULL,
                      allow_for_mixture_models = FALSE,
                      dirichlet_prior = 0.001,
                      component_truncation = 7,
                      regime_truncation = 15,
                      hyperprior_b = 20,
                      model_params_save_every = 5,
                      simulation_iter = NULL,
                      T2_window_size = 3,
                      determining_p_cutoff = FALSE,
                      prob_cutoff = 0.5,
                      model_log_type = "NoModelSpecified",
                      regime_selection_multiplicative_prior = 2,
                      split_selection_multiplicative_prior = 2,
                      is_initial_fit = TRUE,
                      verbose = FALSE) {
  scale_matrix = NULL
  # If the user wants a verbose output of the entire model, this is saved to a log file (it is very very long).
  if (verbose) {
    cat("It is recommended that you pipe this output to a text file.\n")
  }
  
  # If you are performing fault detection, logs of the models drawn must be saved.  The graph
  # is made at the end of this method so that these logs need not be saved past this method.
  model_saves_list  = list()
  model_saves_names = c()
  model_save_count = 1
  if (!is.null(names(data_woTimeValues))) {
    variable_names = names(data_woTimeValues)
  } else {
    variable_names = as.character(1:ncol(data_woTimeValues))
  }
  
  original_p         = ncol(data_woTimeValues)
  has_missing_values = apply(data_woTimeValues, 2, anyNA)
  names_of_columns   = variable_names
  for (miss_col_index in which(has_missing_values)) {
    rows_with_missing         = is.na(data_woTimeValues[, miss_col_index])
    new_col_name              = paste(names_of_columns[miss_col_index], "NAs", sep =
                                        "")
    data_woTimeValues         = cbind(data_woTimeValues, as.numeric(rows_with_missing))
    variable_names            = c(variable_names, new_col_name)
    not.cont                  = c(not.cont, 1)
  }
  
  # As a sanity check, component truncation should be set to 1 when we do not allow for mixture models:
  if (!allow_for_mixture_models) {
    component_truncation = 1
  }
  # | --------------------- Establish Initial Hyperparameters -------------------------- |
  p = ncol(data_woTimeValues)
  if (verbose)
    cat("inside the fit function \n")
  # Sort days and observations
  data_woTimeValues    = data_woTimeValues[order(time_of_observations),]
  time_of_observations = time_of_observations[order(time_of_observations)]
  time_points          = time_points[order(time_points)]
  data_set_list        = list()
  Z_timepoint_indices  = list(list(
    timepoint_first_index = NA,
    timepoint_last_index = NA
  ))
  negative_count       = 0
  for (i in 1:(length(time_points) - 1)) {
    if (i == 1) {
      if (length(which(time_of_observations <= time_points[i + 1 - negative_count])) < p) {
        if (verbose)
          cat(
            paste(
              "There are not enough observations for time point,",
              i - negative_count,
              "\n"
            )
          )
        index_to_remove = i - negative_count + 1
        time_points     = time_points[-index_to_remove]
        negative_count  = negative_count + 1
        next
      }
      if (verbose)
        cat(
          paste(
            "Number of timepoints in day",
            1,
            "is:",
            length(which(
              time_of_observations <= time_points[i + 1 - negative_count]
            )),
            '\n'
          )
        )
      
      Z_timepoint_indices[[i - negative_count]]$timepoint_first_index = 1
      Z_timepoint_indices[[i - negative_count]]$timepoint_last_index  = max(which(time_of_observations <= time_points[i +
                                                                                                                        1 - negative_count]))
    } else {
      indices_of_timepoint                           = which((time_of_observations > time_points[i -
                                                                                                   negative_count]) &
                                                               (time_of_observations <= time_points[i +
                                                                                                      1 - negative_count]))
      if (length(indices_of_timepoint) < p) {
        if (verbose)
          cat(
            paste(
              "There are not enough observations for time point,",
              i - negative_count,
              '\n'
            )
          )
        index_to_remove = i - negative_count + 1
        time_points         = time_points[-index_to_remove]
        negative_count  = negative_count + 1
        next
      }
      if (verbose)
        cat(
          paste(
            "Number of timepoints in day",
            i - negative_count,
            "is:",
            length(indices_of_timepoint),
            '\n'
          )
        )
      Z_timepoint_indices[[i - negative_count]]                       = list(timepoint_first_index = NA,
                                                                             timepoint_last_index = NA)
      Z_timepoint_indices[[i - negative_count]]$timepoint_first_index = min(indices_of_timepoint)
      Z_timepoint_indices[[i - negative_count]]$timepoint_last_index  = max(indices_of_timepoint)
    }
    data_set_list[[i - negative_count]] = data_woTimeValues[Z_timepoint_indices[[i -
                                                                                   negative_count]]$timepoint_first_index:Z_timepoint_indices[[i - negative_count]]$timepoint_last_index,]
  }
  
  
  if (is.null(previous_states)) {
    time_since_last_change   = 1
    current_regime           = 1
    previous_states          = rep(1, times = (length(time_points) - 1))
  }
  
  my_states = previous_states
  if (verbose)
    cat(
      c("[1]  -----> current states are:",
        my_states,
        "\n")
    )
  
  p              = ncol(data_woTimeValues)
  
  if (is.null(scaleMatrix)) {
    # I'm following the formulation of Lenkowski and Dobra 2011 for the scale matrix.  In contrast
    # to the Wishart Distribution from Wikipedia, the average of their distribution is df * solve(scaleMatrix).
    # Thus, to make the average of this Wishart distribution to be equal to the observed precision matrix, I'll
    # require the following scale matrix:
    sum_empirical_prec    = matrix(0, p, p)
    # IMPORTANT NOTE: in the formulation of Lenkowski and Dobra, the average of the distribution is D^{-1} * (df.prior + p - 1).  Thus, when the
    #                 following scale matrix is inverted, and multiplied by (df.prior + p - 1), what should remain is the empirical average precision.
    regime_index        = 1
    min_index           = Z_timepoint_indices[[min(which(previous_states ==
                                                           regime_index))]]$timepoint_first_index
    max_index           = Z_timepoint_indices[[max(which(previous_states ==
                                                           regime_index))]]$timepoint_last_index
    temp_data           = data_woTimeValues[min_index:max_index,]
    covMatrix           = cov(temp_data, use = "pairwise.complete.obs")
    
    if ((sum(eigen(covMatrix)$values <= 0) > 0)) {
      if (verbose)
        cat(
          "NOTE: We were not able to get a good initial guess for the precision matrix."
        )
      covMatrix         = empirical_cov_w_missing(temp_data, Z_timepoint_indices, previous_states)
    }
    
    if (sum(eigen(covMatrix)$values <= 0) > 0) {
      if (verbose)
        cat(
          "Was unable to get a starting covariance matrix for this data.  Starting with an identity matrix."
        )
      my_func_1             = function(x) {
        sd(x, na.rm = T) ^ 2
      }
      covMatrix             = diag(apply(data_woTimeValues, 2, my_func_1))
    }
    
  }
  if (verbose)
    cat("Empirical Covariance Estimate: \n")
  if (verbose)
    cat(covMatrix)
  if (verbose)
    cat("\n")
  hyperparameters        = list(
    hyperprior_b = hyperprior_b,
    alpha = linger_parameter,
    beta = move_parameter,
    alpha_hyperparameter = linger_parameter,
    beta_hyperparameter = move_parameter,
    mu_0 = NA,
    lambda = lambda,
    p = NA,
    hyperprior_scale_matrix = 2 * diag(p),
    #solve(covMatrix)
    wishart_scale_matrix = covMatrix * (wishart_df_initial +
                                          p - 1),
    #*component_truncation*regime_truncation,
    wishart_df = wishart_df_initial,
    original_wishart_df = wishart_df_initial,
    log_wishart_prior_term = NA,
    dirichlet_prior = dirichlet_prior,
    component_truncation = component_truncation,
    regime_truncation = regime_truncation,
    g.prior = g.prior,
    lambda_hyperparameter_shape = lambda
  )
  # | ------------------------ Fill in Missing Inputs ---------------------------------- |
  data_woTimeValues[is.infinite(data_woTimeValues)] = NA
  model_states                                      = list()
  transition_probabilities                          = rep(NA, times = nrow(data_woTimeValues))
  
  # Here I'm assuming that the first values in time_of_observations is the earliest
  # available time, and the last is the last available time.
  iter = iterations
  if (is.null(lower_bounds)) {
    min_value    = min(data_woTimeValues, na.rm = TRUE) - 1e100
    lower_bounds = rep(min_value, times = ncol(data_woTimeValues))
  }
  if (is.null(upper_bounds)) {
    max_value    = max(data_woTimeValues, na.rm = TRUE) + 1e100
    upper_bounds = rep(max_value, times = ncol(data_woTimeValues))
  }
  if (is.null(not.cont)) {
    not.cont     = rep(0, times = ncol(data_woTimeValues))
  }
  
  if (is.null(ordinal_indicators)) {
    ordinal_indicators     = rep(0, times = ncol(data_woTimeValues))
    list_of_ordinal_levels = list()
  }
  
  if (is.null(categorical_indicators)) {
    categorical_indicators = rep(0, times = ncol(data_woTimeValues))
  }
  
  if (length(not.cont) < ncol(data_woTimeValues)) {
    not.cont = c(not.cont, rep(1, times = (
      ncol(data_woTimeValues) - length(not.cont)
    )))
  }
  
  is_continuous            = !not.cont
  hyperparameters$not.cont = not.cont
  
  if (sum(apply(data_woTimeValues, 2, class) == "numeric") != ncol(data_woTimeValues)) {
    stop("Data must be all numeric.  Makes sure that you've seperated the POSIXct objects.")
  }
  if (!all((
    class(time_of_observations) %in% c("POSIXct", "POSIXlt" , "POSIXt")
  ))) {
    stop("time_of_observations must contain only POSIXt objects.")
  }
  
  # | --------------------------- Gathering Data Values -------------------------------------- |
  if (wishart_df_initial < 3)
    stop(" 'hyperprior.df' must be >= 3.")
  if (iter < burnin)
    stop(" Number of iteration must be more than number of burn-in.")
  burnin = floor(burnin)
  
  # Gather indicator information used for the latent variable draws
  n                              = nrow(data_woTimeValues)
  p                              = ncol(data_woTimeValues)
  hyperparameters$p              = p
  lower_bound_is_equal           = data_woTimeValues
  upper_bound_is_equal           = data_woTimeValues
  is_missing                     = data_woTimeValues
  mean_wo_missing                = function(x) {
    mean(x, na.rm = T)
  }
  mu_0                           = rep(0, times = p)
  
  # If I don't have the sampling distribution over the graph structure set, start it with a
  # uniform distribution.
  if (verbose)
    cat("establishing the starting graph structure.\n")
  
  if ((is.null(set_G)) & (g.prior < 1)) {
    # If G is not set, I'll draw from the distribution over all possible graphs determined by g.prior,
    # the probability of a given edge being present.
    full_graph_size = 0.5 * p * (p - 1)
    data.sim.mixed_group1 = BDgraph::graph.sim(
      p = p,
      graph = "scale-free",
      size = floor(g.prior * full_graph_size)
    )
    previous_G       = data.sim.mixed_group1
    
    temp_link_list = BDgraph::adj2link(previous_G)
    previous_G_link_list = ess::make_null_graph(nodes = as.character(1:p))
    for (p_index in 1:p) {
      previous_G_link_list[[p_index]] = as.character(temp_link_list[which(temp_link_list[, 1] ==
                                                                            p_index), 2])
    }
    count = 1
    
    # This code checks to see that we are starting with a decomposable graph.
    while ((count < 500) &
           (!ess::is_decomposable(previous_G_link_list))) {
      data.sim.mixed_group1 = BDgraph::graph.sim(
        p = p,
        graph = "scale-free",
        size = floor(g.prior * full_graph_size)
      )
      previous_G       = data.sim.mixed_group1
      
      temp_link_list = BDgraph::adj2link(previous_G)
      previous_G_link_list = ess::make_null_graph(nodes = as.character(1:p))
      for (p_index in 1:p) {
        previous_G_link_list[[p_index]] = as.character(temp_link_list[which(temp_link_list[, 1] ==
                                                                              p_index), 2])
      }
      
      count = count + 1
    }
    if (count == 500) {
      stop(
        "We were unable to sample a starting decomposable graph structure.  Provide your own or consider using the full model."
      )
    }
    
    # This code tries to get a graph with a size closer to our g.prior amount.
    count = 1
    while ((count < 500) &
           (((sum(previous_G) / 2) / full_graph_size) <= g.prior)) {
      previous_G_temp                   = previous_G
      # Randomly choose a link to add:
      previous_G[lower.tri(previous_G)] = 1
      previous_G                        = previous_G + diag(p)
      locations_of_zero                 = which(previous_G == 0)
      random_location                   = sample(locations_of_zero, 1)
      previous_G[random_location]       = 1
      previous_G                        = previous_G - diag(diag(previous_G))
      previous_G[lower.tri(previous_G)] = 0
      previous_G[lower.tri(previous_G)] = t(previous_G)[lower.tri(previous_G)]
      
      temp_link_list = BDgraph::adj2link(previous_G)
      previous_G_link_list = ess::make_null_graph(nodes = as.character(1:p))
      for (p_index in 1:p) {
        previous_G_link_list[[p_index]] = as.character(temp_link_list[which(temp_link_list[, 1] ==
                                                                              p_index), 2])
      }
      
      if (!ess::is_decomposable(previous_G_link_list)) {
        previous_G = previous_G_temp
      }
      count = count + 1
    }
    
    
    
  } else if ((is.null(set_G)) & (g.prior == 1)) {
    previous_G = matrix(1, p, p) - diag(p)
  } else {
    previous_G    = set_G
  }
  hyperparameters$G = previous_G
  
  if (is.null(g_sampling_distribution)) {
    g_sampling_distribution = matrix(g.prior, nrow = p, ncol = p)
    g_sampling_distribution = (1 - g_sampling_distribution) * (previous_G) + g_sampling_distribution *
      (1 - previous_G)
  }
  
  
  for (index in 1:p) {
    lower_bound_is_equal[, index] = (data_woTimeValues[, index] == lower_bounds[index])
    upper_bound_is_equal[, index] = (data_woTimeValues[, index] == upper_bounds[index])
    is_missing[, index]           = is.na(data_woTimeValues[, index])
    if (sum(is.na(data_woTimeValues[, index])) == nrow(data_woTimeValues)) {
      stop("Drop the data column with all missing values!")
    }
    mu_0[index]                  = mean(data_woTimeValues[, index], na.rm = TRUE)
  }
  lower_bound_is_equal[is.na(lower_bound_is_equal)] = 0
  upper_bound_is_equal[is.na(upper_bound_is_equal)] = 0
  
  hyperparameters$mu_0                   = mu_0
  hyperparameters$is_missing             = is_missing
  hyperparameters$lower_bound_is_equal   = lower_bound_is_equal
  hyperparameters$upper_bound_is_equal   = upper_bound_is_equal
  hyperparameters$lower_bounds           = lower_bounds
  hyperparameters$upper_bounds           = upper_bounds
  my_states                              = previous_states
  
  # | --------------------------- Initial Model Fits -------------------------------------- |
  transition_probabilities_temp = redraw_transition_probs(previous_states,
                                                          linger_parameter,
                                                          move_parameter,
                                                          n.cores - 1)
  transition_probabilities      = c(transition_probabilities_temp,
                                    rep(NA, times = (
                                      length(my_states)  - length(transition_probabilities_temp)
                                    )))
  if (verbose)
    cat("starting to fit the models for each state \n")
  
  if (verbose)
    cat(
      "data has been subset -- setting up parallelization \n"
    )
  # n.cores_eachparent = max(min(regime_truncation, floor(n.cores/4)), 1)
  # n.cores_eachchild  = max(4, floor(n.cores/n.cores_eachparent))
  n.cores_eachparent = max(n.cores - 1, 1)
  n.cores_eachchild  = 1
  if (verbose)
    cat(
      paste(
        "Setting up",
        n.cores_eachparent,
        "children, each can fork",
        n.cores_eachchild,
        "additional children.\n"
      )
    )
  
  if (is.null(previous_model_fits)) {
    # If we don't have model fits from a previous run, we need to fit models
    # based off of the states.
    df_prior_on_wish_start = wishart_df_initial
    previous_model_fits    = lapply(1:(regime_truncation + 1),
                                    function(x) {
                                      list(
                                        precision = list(),
                                        mu = list(),
                                        cluster_assignments = NA,
                                        component_log_probs = NA,
                                        component_sticks = NA
                                      )
                                    })
    models_temp = list()
    for (x in 1:(regime_truncation + 1)) {
      linger_parameter         = hyperparameters$alpha
      move_parameter           = hyperparameters$beta
      not.cont                 = hyperparameters$not.cont
      
      if ((sum(hyperparameters$G) / 2) >= (0.5 * p * (p - 1))) {
        tryCatch({
          result              = stats::rWishart(
            1,
            (hyperparameters$hyperprior_b + hyperparameters$p - 1),
            solve(hyperparameters$hyperprior_scale_matrix)
          )[, , 1]
          new_scale           = result
        }, error = function(cond) {
          result = scale_matrix
          flag   = 1
          message("Here's the original error message:")
          message(cond)
        })
        
        
      } else {
        result                   = rgwish_Rcpp(
          as.double(previous_G),
          as.double(hyperparameters$hyperprior_scale_matrix),
          as.integer(hyperparameters$hyperprior_b),
          as.integer(hyperparameters$p),
          as.double(1e-8)
        )
        result                   = result[['K']]
        new_scale                = matrix(result, hyperparameters$p, hyperparameters$p)
      }
      
      previous_model_fits      = lapply(1:(regime_truncation + 1),
                                        function(x) {
                                          list(
                                            precision = list(),
                                            mu = list(),
                                            cluster_assignments = NA,
                                            component_log_probs = NA,
                                            component_sticks = NA
                                          )
                                        })
      
      
      
      if (sum(my_states == x) > 0) {
        first_index = Z_timepoint_indices[[min(which(my_states == x))]]$timepoint_first_index
        last_index  = Z_timepoint_indices[[max(which(my_states == x))]]$timepoint_last_index
        indicies    = c(first_index:last_index)
        temp_data   = data_woTimeValues[first_index:last_index, ]
      }
      
      # Start by drawing a prior for everything.
      previous_model_fits = redraw_mixture_parameters(
        my_states,
        x,
        previous_model_fits,
        data_woTimeValues,
        Z_timepoint_indices,
        linger_parameter,
        move_parameter,
        not.cont,
        previous_G,
        hyperparameters
      )
      if (sum(my_states == x) > 0) {
        previous_model_fits[[x]]$cluster_assignments = rep(1, times = nrow(temp_data))
      }
      
      models_temp[[x]] = list(
        precision = previous_model_fits[[x]]$precision,
        mu = previous_model_fits[[x]]$mu,
        cluster_assignments = previous_model_fits[[x]]$cluster_assignments
      )
    }
    
    if (verbose)
      cat("initial model fits complete.")
    
    # Now we fill in the results from our setup into previous_model_fits,
    # and do our data subset now to save on time:
    # max_value = c()
    for (i in 1:regime_truncation) {
      previous_model_fits[[i]][["precision"]]           = models_temp[[i]]$precision
      previous_model_fits[[i]][["mu"]]                  = models_temp[[i]]$mu
      previous_model_fits[[i]][["cluster_assignments"]] = models_temp[[i]]$cluster_assignments
      
      # Here, I am going to fill in each of the component assignment 'sticks' from the beta (dirichlet process) prior.
      temp_component_sticks = rep(0, times = component_truncation)
      if (anyNA(previous_model_fits[[i]][["cluster_assignments"]])) {
        # This means that there are no data assigned to this regime.
        previous_model_fits[[i]][["component_log_probs"]]     = rbeta(component_truncation,
                                                                      1,
                                                                      hyperparameters$dirichlet_prior)
      } else {
        for (stick_index in 1:component_truncation) {
          if (stick_index %in% previous_model_fits[[i]][["cluster_assignments"]]) {
            temp_component_sticks[stick_index] = rbeta(
              1,
              1 + sum(previous_model_fits[[i]][["cluster_assignments"]] == stick_index),
              hyperparameters$dirichlet_prior +
                sum(previous_model_fits[[i]][["cluster_assignments"]] > stick_index)
            )
          } else {
            temp_component_sticks[stick_index] = rbeta(1, 1, hyperparameters$dirichlet_prior)
          }
        }
      }
      previous_model_fits[[i]][["component_log_probs"]] = sticks_to_log_probs(temp_component_sticks)
      previous_model_fits[[i]][["component_sticks"]]    = temp_component_sticks
    }
    
  } else {
    if ((ncol(previous_model_fits[[1]]$precision[[1]]) != p) |
        (nrow(previous_model_fits[[1]]$precision[[1]]) != p) |
        (length(previous_model_fits[[1]]$mu[[1]]) != p)) {
      stop(
        "The dimension of the model fit parameters does not match the dimension of the input data.  You may need to consider a new initial_fit."
      )
    }
  }
  
  
  full_data_Z              = data_woTimeValues
  hyperparameters$raw_data = data_woTimeValues
  hyperparameters$G        = previous_G
  
  # | -------------------------------Starting MCMC-------------------------------------- |
  # If the models were not built, now they have been.  We can begin the mcmc
  # updates.
  model_states                    = list()
  transition_probabilities_states = list()
  mergesplit_accepts              = c()
  DRJ_accepts                     = c()
  gibbswap_accepts                = c()
  alphabeta_accepts               = c()
  cholesky_failed_prec            = c()
  G_Wish_portion_trace            = c()
  D_prior_dens_trace              = c()
  my_flag                         = 0
  
  mergesplit_accepted   = 0
  DRJ_accepted          = 0
  gibbswap_accepted     = 0
  alphabeta_accepted    = 0
  record_count          = 1
  
  if (verbose)
    cat("Beginning MCMC Chain \n")
  if (verbose)
    cat(
      paste(
        "-----> alpha:",
        hyperparameters$alpha,
        ", beta:",
        hyperparameters$beta,
        ", lambda:",
        hyperparameters$lambda,
        '\n'
      )
    )
  for (iter in 1:iterations) {
    if (verbose)
      cat(paste("Starting Iteration", iter, "\n"))
    # | --------------------------- Redraw the Latent Data ------------------------------------- |
    if (verbose)
      cat(
        "-----> setting up first-level parallelization to draw from latent data\n"
      )
    if (verbose)
      cat(
        paste("-----> Setting up", n.cores, "cores.\n")
      )
    
    myfunc_2       = function(x_num) {
      return(
        latent_data_parallel_helper(
          x_num,
          Z_timepoint_indices,
          full_data_Z,
          data_woTimeValues,
          is_missing,
          upper_bound_is_equal,
          lower_bound_is_equal,
          previous_model_fits,
          my_states,
          hyperparameters,
          list_of_ordinal_levels,
          ordinal_indicators,
          categorical_indicators
        )
      )
    }
    data_subsets_z = mclapply(1:max(my_states), myfunc_2)
    
    full_data_Z = NULL
    for (i in 1:max(my_states)) {
      full_data_Z = rbind(full_data_Z, data_subsets_z[[i]])
    }
    if (verbose)
      cat(
        paste("-----> Redraw latent data complete.\n")
      )
    
    # | ------------------------------ Update States: Merge-Split --------------------------------- |
    if (verbose)
      cat(
        paste("-----> Redraw hyperparameters before Merge-Split.\n")
      )
    new_hyperparameters                         = redraw_hyperparameters(
      hyperparameters,
      transition_probabilities,
      previous_model_fits,
      my_states,
      verbose = verbose
    )
    hyperparameters                             = new_hyperparameters$hyperparameters_item
    previous_model_fits                         = new_hyperparameters$previous_model_fits_item
    symmetric_scale                             = hyperparameters$wishart_scale_matrix
    symmetric_scale[lower.tri(symmetric_scale)] = 0
    symmetric_scale                             = symmetric_scale + t(symmetric_scale) - diag(diag(symmetric_scale))
    D_prior_dens_trace                          = c(
      D_prior_dens_trace,
      CholWishart::dWishart(
        symmetric_scale,
        hyperparameters$wishart_df,
        diag(hyperparameters$p),
        log = TRUE
      )
    )
    if (verbose)
      cat(
        paste("-----> Hyperparameters redraw completed.\n")
      )
    
    if (iter < floor(iterations / 5)) {
      launching = TRUE
    } else {
      launching = FALSE
    }
    
    if (verbose)
      cat(
        paste("-----> Starting the Merge-Split Algorithm.\n")
      )
    states_update            = update_states_mergesplit(
      my_states,
      previous_model_fits,
      full_data_Z,
      Z_timepoint_indices,
      linger_parameter,
      move_parameter,
      not.cont,
      transition_probabilities,
      previous_G,
      hyperparameters,
      n.cores,
      allow_mixtures = TRUE,
      min_regime_length = 1,
      regime_selection_multiplicative_prior = regime_selection_multiplicative_prior,
      split_selection_multiplicative_prior = split_selection_multiplicative_prior,
      verbose = verbose
    )
    if (verbose)
      cat(
        paste("-----> Merge-Split Algorithm complete.\n")
      )
    previous_model_fits      = states_update$previous_model_fits_item
    my_states                = states_update$my_states_item
    transition_probabilities = states_update$transition_probabilities_item
    mergesplit_accepted      = states_update$accepted_item
    full_data_Z              = states_update$full_data_Z_item
    G_Wish_portion_trace     = c(G_Wish_portion_trace,
                                 states_update$G_Wish_portion_of_MH)
    
    # | --------------------------- Update States: Gibbs Sampler ------------------------------ |
    if (verbose)
      cat(paste("-----> Starting the Gibbs Sweep.\n"))
    
    states_update            = update_states_gibbs(
      my_states,
      full_data_Z,
      Z_timepoint_indices,
      previous_model_fits,
      transition_probabilities,
      hyperparameters,
      n.cores - 1,
      min_regime_length = 2,
      verbose = verbose
    )
    
    previous_model_fits      = states_update$previous_model_fits_item
    my_states                = states_update$my_states_item
    transition_probabilities = redraw_transition_probs(my_states,
                                                       hyperparameters$alpha,
                                                       hyperparameters$beta,
                                                       n.cores - 1)
    gibbswap_accepted        = states_update$accepted_item
    if (verbose)
      cat(paste("-----> Gibbs Sweep complete.\n"))
    
    my_states_all_equal = unique(my_states)
    values_between      = 1:max(my_states)
    if (all.equal(my_states_all_equal, values_between) != TRUE) {
      stop("Gibbs update caused an error with state relabeling.")
    }
    
    
    # | ------------------ Merge-Split-Gibbs the Mixture Components ----------------------- |
    # I need to make sure that the component probabilities are well-behaved here,
    # so I'll update the component probabilities using a conjugate update.
    if (((allow_for_mixture_models) &
         (iter %% 5 == 1)) |
        ((allow_for_mixture_models) &
         (iter < floor(iterations / 4)))) {
      if (verbose)
        cat(
          "about to merge-split-gibbs the mixture components \n"
        )
      myfunc_3               = function(x) {
        return(
          splitmerge_comps_parallel_helper(
            x,
            Z_timepoint_indices,
            my_states,
            full_data_Z,
            previous_model_fits,
            hyperparameters
          )
        )
      }
      new_mixture_components = mclapply(1:max(my_states), myfunc_3)
      
      for (i in 1:max(my_states)) {
        previous_model_fits[[i]]$precision           = new_mixture_components[[i]]$precisions
        previous_model_fits[[i]]$mu                  = new_mixture_components[[i]]$mus
        previous_model_fits[[i]]$cluster_assignments = new_mixture_components[[i]]$assigns
        previous_model_fits[[i]]$component_log_probs = new_mixture_components[[i]]$component_probs
        previous_model_fits[[i]]$component_sticks    = new_mixture_components[[i]]$component_sticks
        previous_model_fits                          = shift_components(i, previous_model_fits)
        if (verbose)
          cat(
            paste("cluster assignments for regime", i, "are:\n")
          )
        if (min(previous_model_fits[[i]]$cluster_assignments) > 1) {
          stop("Component shift failed!!")
        }
      }
    } else if (allow_for_mixture_models) {
      # This is the case where (in between the regular split merge steps that i do every 5 iterations) I
      # just do a quick gibbs swap.
      myfunc_4 = function(x) {
        return(
          gibbsswap_comps_parallel_helper(
            x,
            Z_timepoint_indices,
            full_data_Z,
            hyperparameters,
            my_states,
            previous_model_fits
          )
        )
      }
      
      new_mixture_components = mclapply(1:max(my_states), myfunc_4)
      for (i in 1:max(my_states)) {
        previous_model_fits[[i]]$cluster_assignments = new_mixture_components[[i]]$assigns
        previous_model_fits                          = shift_components(i, previous_model_fits)
      }
    }
    
    # | -------------------------------Double Reversible Jump Graph Resampling---------------------------------------- |
    # First, resample everything based on the posteriors.
    if (verbose)
      cat(
        paste("-----> Starting to refit the model parameters.\n")
      )
    
    # Randomly select an edge based on the sampling distribution on the graph structure.
    if ((g.prior < 1) & (iter %% 5 == 0)) {
      total_prob_mass = 0
      suppressWarnings({
        rand_selection  = runif(1,
                                min = 0,
                                max = sum(g_sampling_distribution))
      })
      
      if (is.na(sum(g_sampling_distribution))) {
        if (verbose)
          cat(
            "Something is wrong with the g sampling distribution! \n"
          )
        if (verbose)
          cat('\n')
        stop("Check the BDgraph package access.")
      }
      done = FALSE
      for (g_index_i in 1:(p - 1)) {
        for (g_index_j in (g_index_i + 1):p) {
          total_prob_mass = total_prob_mass + g_sampling_distribution[g_index_i, g_index_j]
          if (total_prob_mass >= rand_selection) {
            selected_edge_i = g_index_i
            selected_edge_j = g_index_j
            done = TRUE
            break
          }
        }
        if (done) {
          break
        }
      }
      if (!done) {
        selected_edge_i = 2
        selected_edge_j = 3
      }
      new_G_proposed  = previous_G
      if (previous_G[selected_edge_i, selected_edge_j]) {
        new_G_proposed[selected_edge_i, selected_edge_j] = 0
        new_G_proposed[selected_edge_j, selected_edge_i] = 0
      } else {
        new_G_proposed[selected_edge_i, selected_edge_j] = 1
        new_G_proposed[selected_edge_j, selected_edge_i] = 1
      }
      
      # Next, we fix this graph G and sample new precision matrices K1,...,Km for each regime.
      # n.cores_eachparent = max(min(max(my_states), n.cores - 1), 1)
      # n.cores_eachchild  = max(floor((n.cores - n.cores_eachparent) / n.cores_eachparent), 1)
      if (verbose)
        cat(
          "-----> Setting up to draw new proposal precision matrices \n"
        )
      
      ## First, we'll just draw new parameter values for every state, which we will accept deterministically by conjugate update.
      temp_previous_G_link_list = ess::make_null_graph(nodes = as.character(1:p))
      adj_new_G_proposed = BDgraph::adj2link(new_G_proposed)
      for (p_index in 1:p) {
        temp_previous_G_link_list[[p_index]] = as.character(adj_new_G_proposed[which(adj_new_G_proposed[, 1] ==
                                                                                       p_index), 2])
      }
      
      ## Next, we update our graph structure G using a Double Reversible Jump
      if ((ess::is_decomposable(temp_previous_G_link_list))) {
        ## Next, we update our graph structure G using a Double Reversible Jump
        
        tryCatch({
          models_temp = NULL
          models_temp = list()
          for (x in 1:max(my_states)) {
            # library(palliative.changepoints)
            
            state_to_redraw = x
            current_graph_G = previous_G
            redraw_G_output = redraw_G_with_mixture(
              my_states,
              state_to_redraw,
              previous_model_fits,
              full_data_Z,
              Z_timepoint_indices,
              current_graph_G,
              hyperparameters,
              new_G_proposed,
              selected_edge_i,
              selected_edge_j,
              g_sampling_distribution
            )
            
            models_temp[[x]] = list(
              new_previous_model_fits_item = redraw_G_output$previous_model_fits_item,
              MH_value_item                = redraw_G_output$log_MH_ratio
            )
            
          }
          #     # We now evaluate the metropolis-hastings ratio.
          if (length(models_temp[[1]]$MH_value_item) != 0) {
            cholesky_failed_prec     = c(cholesky_failed_prec, 0)
            log_MHratio_all_models   = 0
            for (index in 1:max(my_states)) {
              log_MHratio_all_models = log_MHratio_all_models + models_temp[[index]]$MH_value_item
            }
            log_unif                             = log(runif(1))
            
            if (is.na(log_MHratio_all_models)) {
              if (verbose)
                cat(
                  "-----> DRJ model fits complete (did not accept).  Missing values! \n"
                )
              DRJ_accepted                             = 0
            } else if ((log_unif <= log_MHratio_all_models) &
                       (!is.infinite(log_MHratio_all_models))) {
              if (verbose)
                cat(
                  "-----> DRJ model fits complete (accepted)."
                )
              
              # Here, we ACCEPT
              for (index in 1:max(my_states)) {
                temp_previous_model_fits               = models_temp[[index]]$new_previous_model_fits_item
                previous_model_fits[[index]]           = temp_previous_model_fits[[index]]
              }
              previous_G                               = new_G_proposed
              #
              DRJ_accepted                             = 1
              
              if (verbose)
                cat(
                  paste(
                    "-----> DRJ model fits DRJ-MH ratio:",
                    log_MHratio_all_models
                  )
                )
            } else {
              # Note that I'm rejecting all infinite MH ratio values.
              if (verbose)
                cat(
                  "-----> DRJ model fits complete (did not accept)."
                )
              
              DRJ_accepted                             = 0
              if (verbose)
                cat(
                  paste(
                    "-----> DRJ model fits DRJ-MH ratio:",
                    log_MHratio_all_models
                  )
                )
            }
          } else {
            if (verbose)
              cat(
                "-----> Cholesky Failed!  DRJ model fits incomplete. (did not accept)."
              )
            cholesky_failed_prec     = c(cholesky_failed_prec, 0)
            log_MHratio_all_models   = NA
          }
        },
        error = function(cond) {
          if (verbose)
            cat(
              "-----> DRJ model fits complete (did not accept due to error)."
            )
          log_MHratio_all_models = -Inf
          DRJ_accepted                             = 0
          if (verbose)
            cat(
              paste(
                "-----> DRJ model fits DRJ-MH ratio:",
                log_MHratio_all_models
              )
            )
        })
        
      } else {
        nondecomposable_count = nondecomposable_count + 1
        if (verbose)
          cat("proposed graph was nondecomposable!")
      }
      
    }
    hyperparameters$G                      = previous_G
    
    # | ----------------------------- Print Hyperparameters ---------------------------------------- |
    if (verbose)
      cat(
        paste(
          "-----> alpha:",
          hyperparameters$alpha,
          ", beta:",
          hyperparameters$beta,
          ", lambda:",
          hyperparameters$lambda,
          '\n'
        )
      )
    if (verbose)
      cat(
        c("[1]  -----> current mu_0 is:", hyperparameters$mu_0),
        "\n"
      )
    if (length(my_states) > 30) {
      if (verbose)
        cat(
          c("[1]  -----> current states are:",
            paste(
              sapply(my_states, paste, collapse = ' ')
            ),
            "\n")
        )
    } else {
      if (verbose)
        cat(
          c("[1]  -----> current states are:",
            paste(
              sapply(my_states, paste, collapse = ' ')
            ),
            "\n")
        )
    }
    # | ----------------------------- Record MCMC Information ---------------------------------------- |
    mergesplit_accepts                      = c(mergesplit_accepts, mergesplit_accepted)
    if (verbose){
      cat("MERGSPLIT ACCEPTS \n")
      cat(mergesplit_accepts)
      cat("\n")
    }
    DRJ_accepts                             = c(DRJ_accepts, DRJ_accepted)
    if (verbose){
      cat("NUMBER OF DRJ ACCEPTS IS:\n")
    }
    if (verbose)
      cat(sum(DRJ_accepts))
    if (verbose)
      cat("\n")
    gibbswap_accepts                        = c(gibbswap_accepts, gibbswap_accepted)
    alphabeta_accepts                       = c(alphabeta_accepts, alphabeta_accepted)
    transition_probabilities_states[[iter]] = transition_probabilities
    
    if (verbose)
      cat(paste("saving data at location", iter, "\n"))
    if ((iter %% model_params_save_every == 0) & (iter > burnin)) {
      # I want to save model parameters associated with different state vectors for later analysis.
      # I am only going to save the model parameters for the final two regimes.
      if (length(unique(my_states)) == 1) {
        temp_list = previous_model_fits
      } else {
        temp_list = previous_model_fits[(max(my_states) - 1):max(my_states)]
      }
      
      # I think that the only information that I really want is the location of the LAST changepoint.
      # I want every model that puts a changepoint here so that I can see if they are all putting one
      # there for the same reason. I don't need an exact match.
      temp_states = rep(0, times = (length(my_states) - 1))
      for (j in 1:(length(my_states) - 1)) {
        temp_states[j] = my_states[j + 1] - my_states[j]
      }
      if (sum(temp_states) > 1) {
        indicies_of_changepoint                                               = which(temp_states ==
                                                                                        1)
        temp_states                                                           = rep(0, times = (length(my_states) -
                                                                                                  1))
        temp_states[indicies_of_changepoint[length(indicies_of_changepoint)]] = 1
      }
      data_name = paste(paste(temp_states, collapse = ""),
                        "_",
                        iter,
                        '.RData',
                        sep = "")
      
      
      # Usually, I only need to save the model information for the last two regimes (fault detection)
      # However, when I am searching for a reasonable probability cutoff on the posterior, I need them all.
      if (model_params_save_every < iterations) {
        model_saves_list[[model_save_count]] = temp_list
        model_saves_names = c(model_saves_names, data_name)
        model_save_count  = model_save_count + 1
      }
    }
    
    if (determining_p_cutoff) {
      data_name = paste(
        "Cuttoff_Model_Logs/",
        paste(my_states, collapse = "n"),
        "_",
        simulation_iter,
        "_",
        iter,
        '.RData',
        sep = ""
      )
      # saveRDS(previous_model_fits, file = data_name)
    }
    
    model_states[[iter]]         = my_states
    record_count                 = record_count + 1
    
    five_percent = floor(iterations * 0.05)
    if (iter %% five_percent == 0) {
      percentile_location = floor(iter / five_percent) * 5
      if (percentile_location == 100) {
        cat(paste(percentile_location, "\n", sep = ""))
      } else{
        cat(paste(percentile_location, "->", sep = ""))
      }
      
    } else if (iter == 1) {
      cat("MCMC chain running...\n")
    }
  }
  hyperparameters$current_G = previous_G
  
  # Grab the changepoint probabilities.
  my_states = model_states
  states_df = data.frame(matrix(my_states[[1]], nrow = 1))
  for (i in 2:length(my_states)) {
    states_df = rbind(states_df, data.frame(matrix(my_states[[i]], nrow = 1)))
  }
  states_df_burnt = states_df[floor(2 * nrow(states_df) / 4):nrow(states_df), ]
  prop_mat = matrix(0, nrow = 1, ncol = (ncol(states_df_burnt) - 1))
  number_of_segments = nrow(states_df)
  for (j in 1:(ncol(states_df_burnt) - 1)) {
    temp_difference = states_df_burnt[, j + 1] - states_df_burnt[, j]
    prop_mat[1, j] = mean(temp_difference == 1)
  }
  time_points        = matrix(1:ncol(prop_mat),
                              ncol = ncol(prop_mat),
                              nrow = 1)
  prop_mat           = data.frame(t(rbind(time_points, prop_mat)))
  names(prop_mat)    = c("Time-point", "Change-point probability")
  
  mylist = list(
    changepoint_probabilities = prop_mat,
    models = previous_model_fits,
    states = model_states,
    transition_probs = transition_probabilities_states,
    mergesplit_accepts = mergesplit_accepts,
    DRJ_accepts = DRJ_accepts,
    gibbswap_accepts = gibbswap_accepts,
    alphabeta_accepts = alphabeta_accepts,
    G_Wish_portion_trace = G_Wish_portion_trace,
    D_prior_dens_trace = D_prior_dens_trace,
    cholesky_failed_prec = cholesky_failed_prec,
    hyperparameters = hyperparameters,
    latent_data = full_data_Z,
    raw_data = data_woTimeValues,
    Z_timepoint_indices = Z_timepoint_indices,
    variable_names = variable_names,
    prob_cutoff = prob_cutoff
  )
  
  class(mylist) = "bayesWatch"
  
  if (model_params_save_every < iterations) {
    names(model_saves_list) = model_saves_names
    differences_of_marginals = determine_marginals(mylist, model_saves_list, prob_cutoff)
    fault_plot = graph_posterior_distributions(differences_of_marginals,
                                               mylist,
                                               model_saves_list,
                                               prob_cutoff,
                                               f_divergence = "Hellinger")
  }
  
  mylist$fault_graph = fault_plot
  
  return(mylist)
}


#' Create an estimate on posterior distribution of change-points.
#'
#' @description  Given a bayesWatch object and a probability cutoff, finds change-points.
#'
#' @param regime_fit_object bayesWatch object. Fit with the bayesWatch method.
#' @param prob_cutoff float in (0,1). Posterior probabilities above this cutoff will be considered changepoints.
#'
#' @return vector. Indicator values corresponding to change-point locations.
#' @export
#'
get_point_estimate = function(regime_fit_object, prob_cutoff) {
  ### Grab the change-point probabilities
  yy = regime_fit_object
  MVN_model_mergesplit_accepts = yy$mergesplit_accepts
  my_states = yy$states
  states_df = data.frame(matrix(my_states[[1]], nrow = 1))
  for (i in 2:length(my_states)) {
    states_df = rbind(states_df, data.frame(matrix(my_states[[i]], nrow = 1)))
  }
  states_df_burnt = states_df[floor(2 * nrow(states_df) / 4):nrow(states_df), ]
  prop_mat = matrix(0, nrow = 1, ncol = (ncol(states_df_burnt) - 1))
  number_of_segments = nrow(states_df)
  for (j in 1:(ncol(states_df_burnt) - 1)) {
    temp_difference = states_df_burnt[, j + 1] - states_df_burnt[, j]
    prop_mat[1, j] = mean(temp_difference == 1)
  }
  changepoint_probs   = prop_mat
  changepoints        = changepoint_probs >= prob_cutoff
  changepoints        = data.frame(matrix(
    changepoints,
    nrow = 1,
    ncol = length(changepoints)
  ))
  return(changepoints)
}


#' Print function for a bayesWatch object.  Prints only the posterior change-point probabilities.
#'
#'
#' @param x bayesWatch object. Fit from bayesWatch main method.
#' @param ... Additional plotting arguments.
#'
#' @usage \method{plot}{bayesWatch}(x) <- value
#' @exportS3Method
plot.bayesWatch = function(x, ...) {
  time_point <- prob_value <- NULL
  regime_fit_object   = x
  prob_cutoff         = regime_fit_object$prob_cutoff
  changepoint_probs   = regime_fit_object$changepoint_probabilities[, 2]
  changepoints_data   = data.frame(prob_value = as.vector(changepoint_probs),
                                   time_point = as.factor(1:length(as.vector(
                                     as.vector(changepoint_probs)
                                   ))))
  
  ggplot2::ggplot(changepoints_data, aes(x = time_point, y = prob_value)) + ggplot2::geom_bar(stat = "identity", fill =
                                                                                                "blue") +
    ggplot2::geom_hline(yintercept = prob_cutoff, color = "red") +
    ggplot2::annotate("text", label = "Probability Cutoff Value")
}


#' Determine the cause of a change-point.
#'
#' @description Prints out fault detection graphics given a bayesWatch object. This method can only be run
#' if fault detection was run on the bayesWatch fit (if model_params_save_every < iterations).
#'
#' @param regime_fit_object bayesWatch object.  Fit with main method of package.
#'
#' @return ggplot object. Fault detection graphs.
#' @export
#'
#'
detect_faults = function(regime_fit_object) {
  if (is.null(regime_fit_object$fault_graph)) {
    stop("This regime fit was run without fault detection.")
  } else {
    regime_fit_object$fault_graph
    return(regime_fit_object$fault_graph)
  }
}
