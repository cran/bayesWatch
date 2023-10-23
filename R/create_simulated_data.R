# Create simulated data with artificial change-point imposed.
library(CholWishart)
library(MASS)

#' Create simulated data used in the bayesWatch package.  This method copies much of the
#' code used in the original paper (Murph et al 2023).  Much of the possible simulated examples
#' are disabled here (we consider only a mean change example).
#'
#' @noRd
create_simulated_data = function(){
  job_num = 7
  parameters_number         = 8
  p                         = 5
  orig_p = p
  discrete_vals = c(1,5)
  # discrete_vals = c()
  has_discretes = 1
  
  mean_shift                = 0.001
  var_shift                 = 0.6 
  prob_cutoff               = 0.1
  half_dist_between_modes   = mean_shift
  percent_missing_index     = parameters_number
  base_data_size            = 50
  
  # Presently, the maximum value for num_of_days is 51.  This can be easily updated by updating the
  # day_dts value below.
  num_of_days               = 10
  
  possible_missings_percent = c(20,0,20,20,20,0,20,20)
  missing_nonmissing        = c(1,0,1,1,1,0,1,1)
  distribution_change       = c('miss','mean','cov','cov','cov','mean','miss', 'mean')
  modal_type                = c('unimodal','unimodal','unimodal','unimodal','bimodal','bimodal','unimodal','unimodal')
  MAR_list                  = c(0,1,1,1,1,1,0,1)
  MAR                       = MAR_list[parameters_number]
  
  should_include_missings   = 1
  
  orig_job_num = job_num
  port_number = (job_num-1)*31 + 11001
  
  not.cont = rep(0, times = p)
  
  graph_structure       = matrix(1,p,p) - diag(p)
  
  df      = p + 50
  S1      = CholWishart::rInvWishart(1, df, diag(p)/df )[,,1]
  orig_S2 = S1
  # S1      = S1* 0.25
  
  S2      = S1
  # S2[3,3] = S2[3,3]*2
  # S2[4,4] = S2[4,4]*4
  S3      = S1
  mean1   = rep(0, p)
  mean2   = mean1 
  mean3   = mean1
  
  mean2[3:4] = mean1[3:4] + mean_shift
  mean2[4]   = mean1[4] + 2*mean_shift
  data_type = "MeanIsDifferent"
  
  day_dts = seq(as.POSIXct("2020-03-01 01:00:00", tz="CST6CDT"),as.POSIXct("2020-04-22 01:00:00", tz="CST6CDT"),by="1 day")
  day_dts = day_dts[1:(num_of_days+1)]
  
  full_data            = NULL
  count                = 0
  day_of_observations  = NULL
  data_set_list        = list()
  for(day in 1:(length(day_dts)-1)){
    count=count+1
    if(is.null(full_data)){
      
      if(modal_type[parameters_number]=='unimodal'){
        # data.sim.mixed_group1 <- bdgraph.sim( n = 2*base_data_size, p = p,
        #                                       graph = "fixed", 
        #                                       sigma=S1, mean = mean1)
        # full_data             = data.sim.mixed_group1$data
        
        data.sim.mixed_group1 = MASS::mvrnorm(2*base_data_size, mu = mean1, Sigma = S1)
        full_data             = data.sim.mixed_group1
        
        
        
        temp_day_log          = full_data
        day_of_observations   = rep(day_dts[day+1], times = nrow(data.sim.mixed_group1))#$data))
        
      } else {
        # data.sim.mixed_group1 <- bdgraph.sim( n = base_data_size, p = p,
        #                                       graph = "fixed", 
        #                                       sigma=S1, mean = mean1-half_dist_between_modes)
        # full_data = data.sim.mixed_group1$data
        
        data.sim.mixed_group1 = MASS::mvrnorm(base_data_size, mu = mean1-half_dist_between_modes, Sigma = S1)
        full_data = data.sim.mixed_group1
        
        day_of_observations = rep(day_dts[day+1], times = nrow(data.sim.mixed_group1))#$data))
        temp_day_log = full_data
        
        # data.sim.mixed_group3 <- bdgraph.sim( n = base_data_size, p = p,
        #                                       graph = "fixed", 
        #                                       sigma=S1, mean = mean1+half_dist_between_modes)
        # full_data = rbind(full_data, data.sim.mixed_group3$data)
        data.sim.mixed_group3 = MASS::mvrnorm(base_data_size, mu = mean1+half_dist_between_modes, Sigma = S1)
        full_data             = rbind(full_data, data.sim.mixed_group3)
        
        day_of_observations =  c(day_of_observations,rep(day_dts[day+1], times = nrow(data.sim.mixed_group3)))#$data)))
        
        temp_day_log = rbind(temp_day_log, data.sim.mixed_group3)#$data)
      }
      
    }else if(count < 15){
      
      if(modal_type[parameters_number]=='unimodal'){
        # data.sim.mixed_group1 <- bdgraph.sim( n = 2*base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S1, mean = mean1)
        # full_data = rbind(full_data, data.sim.mixed_group1$data)
        
        data.sim.mixed_group1 = MASS::mvrnorm(2*base_data_size, mu = mean1, Sigma = S1)
        full_data = rbind(full_data, data.sim.mixed_group1)
        
        day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group1)))#$data)))
        
        
        temp_day_log = data.sim.mixed_group1#$data
        
      } else {
        # data.sim.mixed_group1 <- bdgraph.sim( n = base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S1, mean = mean1-half_dist_between_modes)
        
        data.sim.mixed_group1 = MASS::mvrnorm(base_data_size, mu = mean1-half_dist_between_modes, Sigma = S1)
        
        full_data = rbind(full_data, data.sim.mixed_group1)#$data)
        temp_day_log = data.sim.mixed_group1#$data
        day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group1)))#$data)))
        # data.sim.mixed_group3 <- bdgraph.sim( n = base_data_size, p = p,
        #                                       graph = "fixed", 
        #                                       sigma=S1, mean = mean1+half_dist_between_modes)
        
        data.sim.mixed_group3 = MASS::mvrnorm(base_data_size, mu = mean1+half_dist_between_modes, Sigma = S1)
        
        full_data = rbind(full_data, data.sim.mixed_group3)#$data)
        day_of_observations =  c(day_of_observations,rep(day_dts[day+1], times = nrow(data.sim.mixed_group3)))#$data)))
        
        temp_day_log = rbind(temp_day_log, data.sim.mixed_group3)#$data)
        
      }
      
    } else if(count < 35) {
      
      
      if(modal_type[parameters_number]=='unimodal'){
        # 
        # if(data_type == "LinearDrift"){
        #   if(distribution_change[parameters_number] == 'mean'){
        #     mean1 = mean1 + mean_shift
        #     mean2 = mean1
        #   } else {
        #     print(paste("count is currently",count))
        #     alpha = (count-15)/21
        #     S2 = solve(solve(S1)*(1-alpha) + solve(orig_S2)*(alpha))
        #   }
        # }
        # data.sim.mixed_group2 = bdgraph.sim( n = 2*base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S2, mean = mean2)
        
        data.sim.mixed_group2 = MASS::mvrnorm(2*base_data_size, mu = mean2, Sigma = S2)
        
        full_data             = rbind(full_data, data.sim.mixed_group2)#$data)
        day_of_observations   = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group2)))#$data)))
        
        temp_day_log = data.sim.mixed_group2#$data
        
      } else {
        if(data_type == "MeanIsDifferent"){
          bimodal_shift = mean_shift
        } else {
          bimodal_shift = 0
        }
        
        
        # data.sim.mixed_group2 <- bdgraph.sim( n = base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S2, mean = mean2-half_dist_between_modes - bimodal_shift)
        
        data.sim.mixed_group2 = MASS::mvrnorm(base_data_size, mu = mean2-half_dist_between_modes - bimodal_shift, Sigma = S2)
        
        full_data = rbind(full_data, data.sim.mixed_group2)#$data)
        temp_day_log = data.sim.mixed_group2#$data
        day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group2)))#$data)))
        
        
        # data.sim.mixed_group4 <- bdgraph.sim( n = base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S2, mean = mean2+half_dist_between_modes + bimodal_shift)
        
        data.sim.mixed_group4 = MASS::mvrnorm(base_data_size, mu = mean2+half_dist_between_modes+bimodal_shift, Sigma = S2)
        
        full_data = rbind(full_data, data.sim.mixed_group4)#$data)
        day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group4)))#$data)))
        
        
        temp_day_log = rbind(temp_day_log, data.sim.mixed_group4)#$data)
      }
      # print(paste("means for day", day, "is:"))
      # print(paste(apply(temp_day_log, 2, mean)))
    } else {
      
      
      if(modal_type[parameters_number]=='unimodal'){
        
        # if(data_type == "LinearDrift"){
        #   if((distribution_change[parameters_number] == 'mean')&&(count==35)){
        #     mean1 = mean1 + 0.2*mean_shift
        #     mean3 = mean1
        #   } else if (distribution_change[parameters_number] == 'mean') {
        #     mean3 = mean1
        #   } else {
        #     S3 = orig_S2
        #   }
        # }
        # data.sim.mixed_group3 <- bdgraph.sim( n = 2*base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S3, mean = mean3)
        
        data.sim.mixed_group3 = MASS::mvrnorm(2*base_data_size, mu = mean3, Sigma = S3)
        
        full_data = rbind(full_data, data.sim.mixed_group3)#$data)
        day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group3)))#$data)))
        
        temp_day_log = data.sim.mixed_group3#$data
        
        
      } else {
        if(data_type == "MeanIsDifferent"){
          bimodal_shift = mean_shift*2
        } else {
          bimodal_shift = 0
        }
        # data.sim.mixed_group3 <- bdgraph.sim( n = base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S3, mean = mean3-half_dist_between_modes - bimodal_shift)
        
        data.sim.mixed_group3 = MASS::mvrnorm(base_data_size, mu = mean3-half_dist_between_modes - bimodal_shift, Sigma = S3)
        
        
        full_data = rbind(full_data, data.sim.mixed_group3)#$data)
        temp_day_log = data.sim.mixed_group3#$data
        day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group3)))#$data)))
        
        # data.sim.mixed_group4 <- bdgraph.sim( n = base_data_size, p = p, 
        #                                       graph = "fixed",
        #                                       sigma=S3, mean = mean3+half_dist_between_modes + bimodal_shift)
        
        data.sim.mixed_group4 = MASS::mvrnorm(base_data_size, mu = mean3+half_dist_between_modes + bimodal_shift, Sigma = S3)
        
        full_data = rbind(full_data, data.sim.mixed_group4)#$data)
        day_of_observations = c(day_of_observations, rep(day_dts[day+1], times = nrow(data.sim.mixed_group4)))#$data)))
        
        
        temp_day_log = rbind(temp_day_log, data.sim.mixed_group4)#$data)
        
        # print(paste("means for day", day, "is:"))
        # print(paste(apply(temp_day_log, 2, mean)))
      }
    }
  }
  
  all_indices = 1:p
  not.cont = as.numeric(1:p %in% discrete_vals)
  
  for(discrete_index in all_indices[which(not.cont == 1)]){
    while_loop_count = 0
    while(TRUE){
      while_loop_count = while_loop_count + 1
      cutoff_value = runif(1, -3,3)
      if(sum((full_data[, discrete_index] >= cutoff_value)) >= (nrow(full_data)-70) ){
      } else if(sum((full_data[, discrete_index] < cutoff_value)) >= (nrow(full_data)-70) ){
      } else if( (mean((full_data[, discrete_index] >= cutoff_value)) > 0.3) & (discrete_index > 10) ){
      } else {
        full_data[,discrete_index] = (full_data[,discrete_index] >= cutoff_value)
        break
      }
    }
    if(while_loop_count>1000){
      stop("WE COULD NOT CREATE THE DISCRETE VARIABLES")
    }
  }
  
  # Mean center non-discrete rows. 
  for(discrete_index in all_indices[which(not.cont == 0)]){
    full_data[,discrete_index] = (full_data[,discrete_index] - 
                                    mean(full_data[,discrete_index],na.rm=T))/sd(full_data[,discrete_index],na.rm=T)
  }
  
  # Add in some missing values.
  count = 0
  missing_nums              = floor((possible_missings_percent[percent_missing_index]/100)*nrow(full_data))
  if(MAR){
    ## Put in missingness, following a simple Missing At Random structure.
    ## I assume that the discrete values are always non-missing.
    continuous_col_numbers = 1:ncol(full_data)
    continuous_col_numbers = continuous_col_numbers[which(not.cont==0)]
    list_of_cols_w_missing = c()
    if(should_include_missings){
      while(count < missing_nums){
        random_i               = sample(1:nrow(full_data), size = 1)
        random_j               = sample(1:ncol(full_data), size = 1)
        list_of_cols_w_missing = c(list_of_cols_w_missing, random_j)
        if(!is.na(full_data[random_i,random_j])){
          full_data[random_i,random_j] = NA
          count = count + 1
        }
      }
      list_of_cols_w_missing = unlist(sort(unique(list_of_cols_w_missing)))
      new_full_data = matrix(0,nrow=nrow(full_data),ncol = (ncol(full_data)+length(list_of_cols_w_missing)) )
      new_full_data[,1:p] = full_data
      
      count = 1
      for(column in list_of_cols_w_missing){
        # 
        new_col = matrix(as.numeric(is.na(full_data[,column]),ncol=1))
        new_full_data[,p+count] = new_col
        count = count + 1
      }
      full_data = new_full_data
      old_p                                = p
      p                                    = ncol(full_data)
      new_graph_structure                  = matrix(0,p,p)
      new_graph_structure[1:old_p,1:old_p] = graph_structure
      graph_structure                      = new_graph_structure
      
      # 
    }
  } else {
    ## Otherwise, let's attempt to add a noteable structure to the missingness pattern itself.
    ### I'll assume n_cols_missing are the columns for which missingness is present.
    ### Thus, the last n_cols_missing are discrete variables, and the next-to-last n_cols_missing
    ### colums are the ONLY columns in which missingness is present.
    # 
    # Set up the missingness so that the last 5 columns represent the next to last 5 column's missingness.
    for(j in 11:15){
      for(i in 1:nrow(full_data)){
        if(full_data[i,j]){
          full_data[i,j-5] = NA
        }
      }
    }
    
  }
  
  ### Fill the data into a list for use by the T2 Control Chart
  lower_index   = 1
  upper_index   = 2*base_data_size
  data_set_list = list()
  for(day in 1:(length(day_dts)-1)){
    data_set_list[[day]] = full_data[lower_index:upper_index,]
    lower_index   = lower_index + 2*base_data_size
    upper_index   = upper_index + 2*base_data_size
  }
  
  # Uncomment these to remake the example data files.
  # save(full_data, file = "data/example_data.RData")
  # save(day_of_observations, file = "data/day_of_observations.RData")
  # save(day_dts, file = "data/day_dts.RData")
  return(list(full_data=full_data, day_of_observations=day_of_observations, day_dts=day_dts))
}
