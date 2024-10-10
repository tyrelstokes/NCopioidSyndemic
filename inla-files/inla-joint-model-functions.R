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
source(here::here("inla-files/joint-helper-funs.R"))
source(here::here("inla-files/formula-funs.R"))
library(INLA)

# Define some priors

U <- 1
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(U, 0.01)))

beta.prior <- list(prec = list(prior = "normal", param = c(1, 1)),
                   initial = 4, fixed = FALSE)

# Set the prior for the space and time effects

h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.7, 0.5))) # time effect

#hlist <- list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),  # spatial effect
#             prec.spatial=list(prior="loggamma",param=c(1,0.001)))

hlist = list(theta1 = list("PCprior", c(2, 0.01)),  
             #  Pr(sd<1) = 0.01,# unlikely to have rr>3just based on the spatial confounding
             theta2 = list("PCprior", c(0.5, 0.1))) 

####################3


get_fe_vec <- function(ind_intercept = T,
                       mar = T,
                       vac = T,
                       n_outcomes_cases,
                       n_outcomes_hosp,
                       n_outcomes_deaths){
  
  n_t <- sum(n_outcomes_cases,n_outcomes_hosp,n_outcomes_deaths)
  

  
  if(ind_intercept == T){
    out <- paste0("int",c(1:n_t))
  }else{
    out <- NULL
  }
  
  if(mar == T){
    out <- c(out,paste0("mar",c(1:n_t)))
  }
  
  if(vac == T){
    out <- c(out,paste0("vac",c(1:n_t)))
  }
  
  out <- out[is.null(out)==F]
  
  out
  
}




valid_zeros <- function(df,
                        out_names,
                        out_likelihoods,
                        thres = 0.075){
  
  out_cols <- df %>% dplyr::select(dplyr::all_of(out_names))
  
 zvec <-  apply(out_cols,2,function(x){
   x <-  mean(x >0, na.rm =T)
   if(x >= thres){
     out <- F
   }else{
     out <- T
   }
   out
  })
 


if(sum(zvec)>0){
  output_names <- out_names[!zvec]
  output_likelihoods <- out_likelihoods[!zvec]
  
  
}else{
  output_names <- out_names
  output_likelihoods <- out_likelihoods
}
  
 out <- list(output_names = output_names,
             output_likelihoods = output_likelihoods,
             zvec = zvec) 
 
 out
}


which_re_effects <- function(cases_likelihoods,
                             hosp_likelihoods,
                             deaths_likelihoods){
  
  total_res <- c("case_zeros","case_poisson","hosp_zeros","hosp_poisson","death_zeros","death_poisson")
  out_res <- total_res
  
  czeros <- grepl("binomial",cases_likelihoods)
  cpois <- !czeros
  if(sum(czeros)<=1){
    out_res <- out_res[-which(out_res == "case_zeros")]
  }
  if(sum(cpois) <= 1){
    out_res <- out_res[-which(out_res == "case_poisson")] 
  }
  
  
  hzeros <- grepl("binomial",hosp_likelihoods)
  hpois <- !hzeros
  
  if(sum(hzeros)<=1){
    out_res <- out_res[-which(out_res == "hosp_zeros")]
  }
  if(sum(hpois) <= 1){
    out_res <- out_res[-which(out_res == "hosp_poisson")] 
  }
  
  
  
  dzeros <- grepl("binomial",deaths_likelihoods)
  dpois <- !dzeros
  
  if(sum(dzeros)<=1){
    out_res <- out_res[-which(out_res == "death_zeros")]
  }
  if(sum(dpois) <= 1){
    out_res <- out_res[-which(out_res == "death_poisson")] 
  }
  
  out_res
}


# Run the inla model on a single county ---------------------
inla_run_county_joint <- function(state_df,
                            county_fips,
                            ar_order = 5,
                            time_intercepts = TRUE,
                            beta.prior = beta.prior,
                            hyper.prec = hyper.prec,
                            prec_fixed = 1,
                            run = T,
                            verbose = T,
                            cases_outcomes = NULL,
                            cases_likelihoods = NULL,
                            hosp_outcomes = NULL,
                            hosp_likelihoods = NULL,
                            deaths_outcomes = NULL,
                            deaths_likelihoods = NULL,
                            cases_group = NULL,
                            hosp_group = NULL,
                            deaths_group = NULL,
                            lag_hosp = T,
                            lag_deaths = T,
                            ind_intercept = T,
                            mar =T,
                            vac = F){
  
  
  df <- state_df%>% dplyr::filter(county_fips_code == county_fips)
  n <- nrow(df)
  
  
  init_list <- NULL
  likelihood_vec <- NULL
  
  if(is.null(cases_outcomes)==FALSE){
    
case_list <- valid_zeros(df = df,
                         out_names = cases_outcomes,
                         out_likelihoods = cases_likelihoods,
                         thres = 0.075)
    
 
cases_outcomes_mod <- case_list$output_names
cases_likelihoods_mod <- case_list$output_likelihoods
cases_rm <- case_list$zvec

cases_group_mod <- cases_group[!cases_rm]
    
 init_list <-  gen_prep(df = df,
                        n = n,
                        out_names = cases_outcomes_mod,
                        init_list = init_list)
 
likelihood_vec <-  append_vec(init_vec = likelihood_vec,
                              new_vec = cases_likelihoods_mod)
  
 
  }
  
  if(is.null(hosp_outcomes) ==FALSE){
    
    if(sum(hosp_outcomes %in% names(df)) != length(hosp_outcomes)){
   df <-   hosp_outcome_create(df)
  
    
    }
    
    
    hosp_list <- valid_zeros(df = df,
                             out_names = hosp_outcomes,
                             out_likelihoods = hosp_likelihoods,
                             thres = 0.075)
    
    
    hosp_outcomes_mod <- hosp_list$output_names
    hosp_likelihoods_mod <- hosp_list$output_likelihoods
    
    hosp_rm <- hosp_list$zvec
    
    hosp_group_mod <- hosp_group[!hosp_rm]
    
    init_list <-  gen_prep(df = df,
                           n = n,
                           out_names = hosp_outcomes_mod,
                           init_list = init_list)
    
    likelihood_vec <-  append_vec(init_vec = likelihood_vec,
                                  new_vec = hosp_likelihoods_mod)
    
  }
  
  if(is.null(deaths_outcomes)==FALSE){
    
    if(sum(deaths_outcomes %in% names(df)) != length(deaths_outcomes)){
      df <- death_outcome_create(df)
      
    }
    
    deaths_list <- valid_zeros(df = df,
                             out_names = deaths_outcomes,
                             out_likelihoods = deaths_likelihoods,
                             thres = 0.075)
    
    
    deaths_outcomes_mod <- deaths_list$output_names
    deaths_likelihoods_mod <- deaths_list$output_likelihoods
    
    deaths_rm <- deaths_list$zvec
    
    deaths_group_mod <- deaths_group[!deaths_rm]
    
    init_list <-  gen_prep(df = df,
                           n = n,
                           out_names = deaths_outcomes_mod,
                           init_list = init_list)
    
    likelihood_vec <-  append_vec(init_vec = likelihood_vec,
                                  new_vec = deaths_likelihoods_mod) 
    
    
  }
  
  
  outcome_list <- init_list$outcome_list
  offset_list <- init_list$offset_list
  id_list <- init_list$id_list
  month_list <- init_list$month_list
  np <- init_list$np  
  month_list_lag <- init_list$month_list_lag
  

  
  
 outcome_type_res <-   outcome_type_intercept(cases_outcomes = cases_outcomes_mod,
                                              hosp_outcomes = hosp_outcomes_mod,
                                              deaths_outcomes = deaths_outcomes_mod,
                                              cases_likelihoods = cases_likelihoods_mod,
                                              hosp_likelihoods = hosp_likelihoods_mod,
                                              deaths_likelihoods = deaths_likelihoods_mod)
    
    inla_lists <-  inla_data_list(df = df,
                                  n = n,
                                  np = np,
                                  outcome_list = outcome_list,
                                  offset_list = offset_list,
                                  id_list = id_list,
                                  month_list = month_list,
                                  extra_time_int = TRUE,
                                  month_list_lag = month_list_lag)
    
    
    # extract the data 
    y <- inla_lists$y
    int_list <- inla_lists$int_list
    id_list <- inla_lists$id_list
    month_list <- inla_lists$month_list
    offset_list <- inla_lists$offset_list
    time_list <- inla_lists$time_lists
    march_list <-  time_list[[1]]
    vac_list <- time_list[[2]]
    month_list_lag <- inla_lists$month_list_lag
    idc_list <- copy_ids(id_list = id_list,
                         new_name = "idc",
                         group_id_list = month_list_lag)
    
    
  link_set <-   link_set_fun(n = n,
                             likelihood_vec = likelihood_vec)
    
    
  data_list <- NULL
    
  
  nl <- length(likelihood_vec)
    
data_list <- append_list_2(initial_list = data_list,
                             new_list = id_list,
                             new_names = paste0("id",c(1:nl)))

data_list <- append_list_2(initial_list = data_list,
                             new_list = int_list,
                             new_names = paste0("int",c(1:nl)))

data_list <- append_list_2(initial_list = data_list,
                             new_list = month_list,
                             new_names = paste0("m",c(1:nl)))

data_list <- append_list_2(initial_list = data_list,
                             new_list = offset_list,
                             new_names = paste0("E",c(1:nl)))

data_list <- append_list_2(initial_list = data_list,
                             new_list = march_list,
                             new_names = paste0("mar",c(1:nl)))
    
data_list <- append_list_2(initial_list = data_list,
                             new_list = vac_list,
                             new_names = paste0("vac",c(1:nl)))

data_list <- append_list_2(initial_list = data_list,
                           new_list = month_list_lag,
                           new_names = paste0("mlag",c(1:nl)))


data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$case_zero_list),
                             new_names = paste0("case_zeros"))

data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$case_poisson_list),
                             new_names = paste0("case_poisson"))


data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$hosp_zero_list),
                             new_names = paste0("hosp_zeros"))

data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$hosp_poisson_list),
                             new_names = paste0("hosp_poisson"))

data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$death_zero_list),
                             new_names = paste0("death_zeros"))

data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$death_poisson_list),
                             new_names = paste0("death_poisson"))


data_list <- append_list_2(initial_list = data_list,
                           new_list = idc_list,
                           new_names = paste0("idc",c(1:nl)))


n_outcomes_cases <- length(cases_outcomes_mod)
n_outcomes_hosp <- length(hosp_outcomes_mod)
n_outcomes_deaths <- length(deaths_outcomes_mod)

fe_vec <- get_fe_vec(ind_intercept = ind_intercept,
                     mar = mar,
                     vac = vac,
                     n_outcomes_cases = n_outcomes_cases,
                     n_outcomes_hosp = n_outcomes_hosp,
                     n_outcomes_deaths = n_outcomes_deaths)


re_effects <- which_re_effects(cases_likelihoods = cases_likelihoods_mod,
                              hosp_likelihoods = hosp_likelihoods_mod,
                              deaths_likelihoods = deaths_likelihoods_mod)

in_form <- formula_county(outcome = "y",
                          rm_int = T,
                          fixed_effects = fe_vec,
                          re_effects = re_effects,
                          lag_spatial_hosp = lag_hosp,
                          lag_spatial_deaths = lag_deaths,
                          n_outcomes_cases = n_outcomes_cases,
                          n_outcomes_hosp = n_outcomes_hosp,
                          n_outcomes_deaths = n_outcomes_deaths)
    
   
      
      test_mod <- INLA::inla(as.formula(in_form),
                             family=c(cases_likelihoods_mod,hosp_likelihoods_mod,deaths_likelihoods_mod),
                             data=data_list,
                             control.compute=list(waic=TRUE),
                             control.predictor = list(compute = TRUE,
                                                      link = link_set),
                             control.inla=list(cmin=0),
                             control.fixed = list(prec = prec_fixed),
                             safe = TRUE,
                             verbose = verbose)
      
      
  pred_list <- inla_predict_gen(inla_mod = test_mod,
                       n = n,
                       np = np,
                       df = df,
                       cases_outcomes =cases_outcomes_mod,
                       hosp_outcomes = hosp_outcomes_mod,
                       deaths_outcomes = deaths_outcomes_mod,
                       cases_likelihoods = cases_likelihoods_mod,
                       hosp_likelihoods = hosp_likelihoods_mod,
                       deaths_likelihoods = deaths_likelihoods_mod,
                       cases_group_mod = cases_group_mod,
                       hosp_group_mod = hosp_group_mod,
                       deaths_group_mod = deaths_group_mod,
                       y = y,
                       which_pred = "mean",
                       month_id_vec = df$month_id,
                       county_id_vec = df$county_fips_code)
      
      
  avg_pred_df <- pred_list$avg_pred_df
  ind_outcomes_prediction_list <- pred_list$pred_list
  grouped_by_type_prediction_list <- pred_list$grouped_list
  
  
  out <- list(avg_pred_df = avg_pred_df,
              ind_outcomes_prediction_list = ind_outcomes_prediction_list,
              grouped_by_type_prediction_list = grouped_by_type_prediction_list,
              inla_mod = test_mod,
              data_list = data_list)
  
  
  out
      
      
      
    }
    



###########################3
###############################
#####################################
# Run on state -----------------------
#####################################3


inla_run_joint <-        function(df,
                                  spatial_model = 'bym2',
                                  time_intercepts = TRUE,
                                  beta.prior = beta.prior,
                                  hyper.prec = hyper.prec,
                                  prec_fixed = 1,
                                  run = T,
                                  verbose = T,
                                  cases_outcomes = NULL,
                                  cases_likelihoods = NULL,
                                  hosp_outcomes = NULL,
                                  hosp_likelihoods = NULL,
                                  deaths_outcomes = NULL,
                                  deaths_likelihoods = NULL,
                                  cases_group = NULL,
                                  hosp_group = NULL,
                                  deaths_group = NULL,
                                  lag_hosp = F,
                                  lag_deaths = F,
                                  ind_intercept = T,
                                  mar =T,
                                  vac = F,
                                  ar_order = 5){
  
  
  n <- nrow(df)
  
  
  init_list <- NULL
  likelihood_vec <- NULL
  
  if(is.null(cases_outcomes)==FALSE){
    
    case_list <- valid_zeros(df = df,
                             out_names = cases_outcomes,
                             out_likelihoods = cases_likelihoods,
                             thres = 0.075)
    
    
    cases_outcomes_mod <- case_list$output_names
    cases_likelihoods_mod <- case_list$output_likelihoods
    cases_rm <- case_list$zvec
    
    cases_group_mod <- cases_group[!cases_rm]
    
    init_list <-  gen_prep(df = df,
                           n = n,
                           out_names = cases_outcomes_mod,
                           init_list = init_list)
    
    likelihood_vec <-  append_vec(init_vec = likelihood_vec,
                                  new_vec = cases_likelihoods_mod)
    
    
  }
  
  if(is.null(hosp_outcomes) ==FALSE){
    
    if(sum(hosp_outcomes %in% names(df)) != length(hosp_outcomes)){
      df <-   hosp_outcome_create(df)
      
      
    }
    
    
    hosp_list <- valid_zeros(df = df,
                             out_names = hosp_outcomes,
                             out_likelihoods = hosp_likelihoods,
                             thres = 0.075)
    
    
    hosp_outcomes_mod <- hosp_list$output_names
    hosp_likelihoods_mod <- hosp_list$output_likelihoods
    
    hosp_rm <- hosp_list$zvec
    
    hosp_group_mod <- hosp_group[!hosp_rm]
    
    init_list <-  gen_prep(df = df,
                           n = n,
                           out_names = hosp_outcomes_mod,
                           init_list = init_list)
    
    likelihood_vec <-  append_vec(init_vec = likelihood_vec,
                                  new_vec = hosp_likelihoods_mod)
    
  }
  
  if(is.null(deaths_outcomes)==FALSE){
    
    if(sum(deaths_outcomes %in% names(df)) != length(deaths_outcomes)){
      df <- death_outcome_create(df)
      
    }
    
    deaths_list <- valid_zeros(df = df,
                               out_names = deaths_outcomes,
                               out_likelihoods = deaths_likelihoods,
                               thres = 0.075)
    
    
    deaths_outcomes_mod <- deaths_list$output_names
    deaths_likelihoods_mod <- deaths_list$output_likelihoods
    
    deaths_rm <- deaths_list$zvec
    
    deaths_group_mod <- deaths_group[!deaths_rm]
    
    init_list <-  gen_prep(df = df,
                           n = n,
                           out_names = deaths_outcomes_mod,
                           init_list = init_list)
    
    likelihood_vec <-  append_vec(init_vec = likelihood_vec,
                                  new_vec = deaths_likelihoods_mod) 
    
    
  }
  
  
  outcome_list <- init_list$outcome_list
  offset_list <- init_list$offset_list
  id_list <- init_list$id_list
  month_list <- init_list$month_list
  np <- init_list$np  
  month_list_lag <- init_list$month_list_lag
  
  
  
  
  outcome_type_res <-   outcome_type_intercept(cases_outcomes = cases_outcomes_mod,
                                               hosp_outcomes = hosp_outcomes_mod,
                                               deaths_outcomes = deaths_outcomes_mod,
                                               cases_likelihoods = cases_likelihoods_mod,
                                               hosp_likelihoods = hosp_likelihoods_mod,
                                               deaths_likelihoods = deaths_likelihoods_mod)
  
  inla_lists <-  inla_data_list(df = df,
                                n = n,
                                np = np,
                                outcome_list = outcome_list,
                                offset_list = offset_list,
                                id_list = id_list,
                                month_list = month_list,
                                extra_time_int = TRUE,
                                month_list_lag = month_list_lag)
  
  
  # extract the data 
  y <- inla_lists$y
  int_list <- inla_lists$int_list
  id_list <- inla_lists$id_list
  month_list <- inla_lists$month_list
  offset_list <- inla_lists$offset_list
  time_list <- inla_lists$time_lists
  march_list <-  time_list[[1]]
  vac_list <- time_list[[2]]
  month_list_lag <- inla_lists$month_list_lag
  idc_list <- copy_ids(id_list = id_list,
                       new_name = "idc",
                       group_id_list = month_list_lag)
  
  
  link_set <-   link_set_fun(n = n,
                             likelihood_vec = likelihood_vec)
  
  
  data_list <- NULL
  
  
  nl <- length(likelihood_vec)
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = id_list,
                             new_names = paste0("id",c(1:nl)))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = int_list,
                             new_names = paste0("int",c(1:nl)))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = month_list,
                             new_names = paste0("m",c(1:nl)))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = offset_list,
                             new_names = paste0("E",c(1:nl)))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = march_list,
                             new_names = paste0("mar",c(1:nl)))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = vac_list,
                             new_names = paste0("vac",c(1:nl)))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = month_list_lag,
                             new_names = paste0("mlag",c(1:nl)))
  
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$case_zero_list),
                             new_names = paste0("case_zeros"))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$case_poisson_list),
                             new_names = paste0("case_poisson"))
  
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$hosp_zero_list),
                             new_names = paste0("hosp_zeros"))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$hosp_poisson_list),
                             new_names = paste0("hosp_poisson"))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$death_zero_list),
                             new_names = paste0("death_zeros"))
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = list(outcome_type_res$death_poisson_list),
                             new_names = paste0("death_poisson"))
  
  
  data_list <- append_list_2(initial_list = data_list,
                             new_list = idc_list,
                             new_names = paste0("idc",c(1:nl)))
  
  
  n_outcomes_cases <- length(cases_outcomes_mod)
  n_outcomes_hosp <- length(hosp_outcomes_mod)
  n_outcomes_deaths <- length(deaths_outcomes_mod)
  
  fe_vec <- get_fe_vec(ind_intercept = ind_intercept,
                       mar = mar,
                       vac = vac,
                       n_outcomes_cases = n_outcomes_cases,
                       n_outcomes_hosp = n_outcomes_hosp,
                       n_outcomes_deaths = n_outcomes_deaths)
  
  
  re_effects <- which_re_effects(cases_likelihoods = cases_likelihoods_mod,
                                 hosp_likelihoods = hosp_likelihoods_mod,
                                 deaths_likelihoods = deaths_likelihoods_mod)
  
 
  
in_form <-   joint_formula(outcome = "y",
                            rm_int = T,
                            fixed_effects = fe_vec,
                            re_effects = re_effects,
                            full_spatial =T,
                            ar_order = ar_order,
                            lag_spatial_hosp = lag_hosp,
                            lag_spatial_deaths = lag_deaths,
                            n_outcomes_cases = n_outcomes_cases,
                            n_outcomes_hosp = n_outcomes_hosp,
                            n_outcomes_deaths = n_outcomes_deaths) 
  
  
  
  ##### copied
  


test_mod <- INLA::inla(as.formula(in_form),
                       family=c(cases_likelihoods_mod,hosp_likelihoods_mod,deaths_likelihoods_mod),
                       data=data_list,
                       control.compute=list(waic=TRUE),
                       control.predictor = list(compute = TRUE,
                                                link = link_set),
                       control.inla=list(cmin=0),
                       control.fixed = list(prec = prec_fixed),
                       safe = TRUE,
                       verbose = verbose)


pred_list <- inla_predict_gen(inla_mod = test_mod,
                              n = n,
                              np = np,
                              df = df,
                              cases_outcomes =cases_outcomes_mod,
                              hosp_outcomes = hosp_outcomes_mod,
                              deaths_outcomes = deaths_outcomes_mod,
                              cases_likelihoods = cases_likelihoods_mod,
                              hosp_likelihoods = hosp_likelihoods_mod,
                              deaths_likelihoods = deaths_likelihoods_mod,
                              cases_group_mod = cases_group_mod,
                              hosp_group_mod = hosp_group_mod,
                              deaths_group_mod = deaths_group_mod,
                              y = y,
                              which_pred = "mean",
                              month_id_vec = df$month_id,
                              county_id_vec = df$county_fips_code)


avg_pred_df <- pred_list$avg_pred_df
ind_outcomes_prediction_list <- pred_list$pred_list
grouped_by_type_prediction_list <- pred_list$grouped_list


out <- list(avg_pred_df = avg_pred_df,
            ind_outcomes_prediction_list = ind_outcomes_prediction_list,
            grouped_by_type_prediction_list = grouped_by_type_prediction_list,
            inla_mod = test_mod,
            data_list = data_list)


out
  
  
}