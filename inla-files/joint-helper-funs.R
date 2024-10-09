# housekeeping ---------------------------------

`%>%` <- dplyr::`%>%`
`%m+%` <- lubridate::`%m+%`
`%do%` <- foreach::`%do%`

# Source helper functions ----------------
source(here::here("inla-files/helper-functions.R"))
source(here::here("inla-files/inla-create-all-variables.R"))
source(here::here("inla-files/inla-prediction-funs.R"))
source(here::here("usa_adj_mat_files/adj-mat-functions.R"))
source(here::here("utilities/date-functions.R"))


# Helper functions (maybe move later) --------


add_list <- function(old_list,
                     new_list){
  
  no <- length(old_list)
  nl <- length(new_list)
  p <- no + nl
  out <- vector("list", length = p)
  
  for(i in 1:p){
    if(i <= no){
      out[[i]] <- old_list[[i]]
    }else{
      out[[i]] <- new_list[[(i-no)]]
    }
  }
  
 out 
}


append_list <- function(init_list,
                        outcome_list,
                        offset_list,
                        id_list,
                        month_list,
                        month_list_lag = NULL){
  
  if(is.null(init_list) ==FALSE){
    outl <- init_list$outcome_list
    offl <- init_list$offset_list
    idl <- init_list$id_list
    ml <- init_list$month_list
    
    
    
    np_old <- length(outl)
    np_new <- length(outcome_list)
    
    np <- np_old + np_new
    
    outcome_list <- add_list(old_list = outl,
                             new_list = outcome_list)
    
    offset_list <- add_list(old_list = offl,
                             new_list = offset_list)
    
   id_list <- add_list(old_list = idl,
                       new_list = id_list)

    month_list <- add_list(old_list = ml,
                           new_list = month_list)
    
    
   # outcome_list <- append(outl,outcome_list)  
   # offset_list <- append(offl,offset_list) 
    
    
   # id_list <- append(idl,id_list)
  #  month_list <- append(ml, month_list)
    
    if(is.null(month_list_lag) ==F){
      mll <- init_list$month_list_lag
      
      month_list_lag <- add_list(old_list = mll,
                             new_list = month_list_lag)
      
     # month_list_lag <- append(mll,month_list_lag)
    }
    
  }else{
    np <- length(outcome_list)
  }
  
  out <- list(outcome_list = outcome_list,
              offset_list = offset_list,
              id_list = id_list,
              month_list = month_list,
              month_list_lag = month_list_lag,
              np = np)
  
  out
  
  
}  



# Joint Model functions at the county level --------------





case_prep <- function(df,
                      n,
                      jhu_zeros = T,
                      init_list = NULL){
  
  if(jhu_zeros == T){
    np <- 3
    outcome_list <- list(df$cases,
                         df$jhu_zero,
                         df$jhu_cases)
    
    offset_list <- rep_vector(np = np,
                              x = df$population_scale)
    
    id_list <- rep_vector(np = np,
                          x = df$remapped_id)
    
    month_list <- rep_vector(np = np,
                             x = df$month_id)
  }else{
    np <- 2
    outcome_list <- list(df$cases,
                         df$jhu_cases)
    
    offset_list <- rep_vector(np = np,
                              x = df$population_scale)
    
    id_list <- rep_vector(np = np,
                          x = df$remapped_id)
    
    month_list <- rep_vector(np = np,
                             x = df$month_id)
    
  }
  
  
  out <-  append_list(init_list = init_list,
                      outcome_list = outcome_list,
                      offset_list = offset_list,
                      id_list = id_list,
                      month_list = month_list)
  
  out
  
}


####### Hospitalization prep

gen_prep <- function(df,
                     n,
                     out_names,
                     init_list = NULL){
  
  out_cols <- df %>% dplyr::select(dplyr::all_of(out_names))
  
  outcome_list <- as.list(out_cols)
  
  np <- ncol(out_cols)
  
  offset_list <- rep_vector(np = np,
                            x = df$population_scale)
  
  id_list <- rep_vector(np = np,
                        x = df$remapped_id)
  
  month_list <- rep_vector(np = np,
                           x = df$month_id)
  
  
  month_list_lag <- rep_vector(np = np,
                               x = ifelse(df$month_id > 1,df$month_id - 1,NA))
  
  
  
  out <-  append_list(init_list = init_list,
                      outcome_list = outcome_list,
                      offset_list = offset_list,
                      id_list = id_list,
                      month_list = month_list,
                      month_list_lag = month_list_lag)
  
  out
  
}


###### outcome vec append ------------------

append_vec <- function(init_vec = NULL,
                       new_vec){
  
  if(is.null(init_vec) == TRUE){
    out <- new_vec
  }else{
    out <- c(init_vec,new_vec)
  }
  out
}

# append list gen function ------------
append_list_gen <- function(initial_list = NULL,
                            new_list,
                            new_names){
  
  if(is.null(initial_list)==FALSE){
    w <- length(new_list)
    
    out <- foreach::foreach(i =1:w)%do%{
      nl <- list(new_list[[i]])
      names(nl) <- new_names[i]
      
      initial_list <- append(initial_list,nl)
    }
    
  }else{
    out <- new_list
    names(out) <- new_names
  }
  out 
  
}

# append list 2

append_list_2 <- function(initial_list = NULL,
                          new_list,
                          new_names){
  
  
  if(is.null(initial_list)==FALSE){
    l_old <- length(initial_list)
    nms_old <- names(initial_list)
    w <- length(new_list)
    
    J <- l_old + w
    
    out <- foreach::foreach(i =1:J)%do%{
      
      if(i <= l_old){
        x <- initial_list[[i]]
      }else{
        x <-  new_list[[(i-l_old)]]
      }
      
      x
      
    }
    
    names(out) <- c(nms_old,new_names)
    
  }else{
    out <- new_list
    names(out) <- new_names
  }
  out 
  
  
}

# link set function --------------

link_set_fun <- function(n,
                         likelihood_vec){
  
  
  p <- length(likelihood_vec)
  
  out <- unlist(lapply(c(1:p),function(i){
    rep(i,n)
  }))
  
  out 
  
}


create_re_intercept <- function(total_vars,
                                n,
                                which_vars,
                                var_name){
  
  int_count <- 1
  
  x_list <- foreach::foreach(i =1:total_vars,.combine = c)%do%{
    
    if(which_vars[i]== 1){
      out <- rep(int_count,n)
      int_count = int_count + 1
    }else{
      out <- rep(NA,n)
    }
    out 
  }
  
  
  
  # inla_vectors(n = n,
  #             np = total_vars,
  #            x_list = x_list,
  #           base_name = var_name)
  
  x_list
}



outcome_type_intercept <- function(cases_outcomes,
                                   hosp_outcomes,
                                   deaths_outcomes,
                                   cases_likelihoods,
                                   hosp_likelihoods,
                                   deaths_likelihoods){
  
  
  all_outcomes <- c(cases_outcomes,hosp_outcomes,deaths_outcomes)
  
  c_l <- length(cases_outcomes)
  
  wc <-  grepl("zero",cases_outcomes)
  
  c_l_z <- sum(wc)
  c_l_p <- c_l - c_l_z
  
  hc <-  grepl("zero",hosp_outcomes)
  
  h_l <- length(hosp_outcomes)
  h_l_z <- sum(hc)
  h_l_p <- h_l - h_l_z
  
  dc <- grepl("zero", deaths_outcomes)
  
  d_l <- length(deaths_outcomes)
  d_l_z <- sum(dc)
  d_l_p <- d_l - d_l_z
  
  
  n_total <- sum(c_l,h_l,d_l)
  
  cs <- c(1:c_l)
  hs <- c((c_l+1):(c_l+h_l))
  ds <- c((c_l+h_l+1):n_total)
  
  w_cz <- rep(0,n_total)
  w_cz[cs][wc] <- rep(1,sum(wc == T))
  w_cp <- rep(0,n_total)
  w_cp[cs][!wc] <- rep(1,sum(wc == F))
  
  
  w_hz <- rep(0,n_total)
  w_hz[hs][hc] <- rep(1,sum(hc == T))
  w_hp <- rep(0,n_total)
  w_hp[hs][!hc] <- rep(1,sum(hc == F))
  
  w_dz <- rep(0,n_total)
  w_dz[ds][dc] <- rep(1,sum(dc == T))
  w_dp <- rep(0,n_total)
  w_dp[ds][!dc] <- rep(1,sum(dc == F))
  
  
  case_zero_list <-   create_re_intercept(total_vars = n_total,
                                          n = n,
                                          which_vars = w_cz,
                                          var_name = "case_zero")
  
  case_poisson_list <- create_re_intercept(total_vars = n_total,
                                           n = n,
                                           which_vars = w_cp,
                                           var_name = "case_poisson")
  
  
  hosp_zero_list <-   create_re_intercept(total_vars = n_total,
                                          n = n,
                                          which_vars = w_hz,
                                          var_name = "hosp_zero")
  
  hosp_poisson_list <- create_re_intercept(total_vars = n_total,
                                           n = n,
                                           which_vars = w_hp,
                                           var_name = "hosp_poisson")
  
  
  death_zero_list <-   create_re_intercept(total_vars = n_total,
                                           n = n,
                                           which_vars = w_dz,
                                           var_name = "death_zero")
  
  death_poisson_list <- create_re_intercept(total_vars = n_total,
                                            n = n,
                                            which_vars = w_dp,
                                            var_name = "death_poisson")
  
  
  out <- list(case_zero_list = case_zero_list,
              case_poisson_list = case_poisson_list,
              hosp_zero_list = hosp_zero_list,
              hosp_poisson_list = hosp_poisson_list,
              death_zero_list = death_zero_list,
              death_poisson_list = death_poisson_list)
  
  out
  
}

#######


# create copied ids

copy_ids <- function(id_list,
                     new_name = "idc"){
  
  n_l <- length(id_list) 
  
  names(id_list) <- paste0(new_name,c(1:n_l))
  
  id_list
  
}





inla_formula_joint <- function(np,
                               spatial = T,
                               spatial_model = "'bym2'",
                               grouped_effect = T,
                               grouped_type = "'ar2'",
                               temporal_model = "ar",
                               order = 5,
                               mar = T,
                               vac = T,
                               late = T){
  
  return(NULL)
  
}
