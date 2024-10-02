# housekeeping ---------------------------------

`%>%` <- dplyr::`%>%`
`%m+%` <- lubridate::`%m+%`

# Source helper functions ----------------
source(here::here("inla-files/helper-functions.R"))
source(here::here("inla-files/inla-create-all-variables.R"))
source(here::here("inla-files/inla-prediction-funs.R"))
source(here::here("usa_adj_mat_files/adj-mat-functions.R"))
source(here::here("utilities/date-functions.R"))
library(INLA)

# Notes
# Add the confidence intervals to the output!!!
# Think about areas where it is possible to extend the Dave and Staci framework!

# y_cases = alpha_cases + S(X,T) + e_cases
# y_hosp = alpha_hosp + B_hosp*S(X) + e_hosp
# y_deaths = alpha_death + B_death*S(X) + e_death

# [1, B_hosp, B_death] * S(X)

#S(X)

# f(S(X))

# 


# default priors -----------------------------------------
U <- 1
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(U, 0.01)))

beta.prior <- list(prec = list(prior = "normal", param = c(1, 1)),
                   initial = 4, fixed = FALSE)

# inla run county function ------------------------------

inla_run_county <- function(state_df,
                            county_fips,
                            ar_order = 5,
                            time_intercepts = TRUE,
                            jhu_zeros = T,
                            beta.prior = beta.prior,
                            hyper.prec = hyper.prec,
                            prec_fixed = 1,
                            run = T,
                            verbose = T){
  
  
  df <- state_df%>% dplyr::filter(county_fips_code == county_fips)
  n <- nrow(df)
  
  if(run == TRUE){
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

  
 inla_lists <-  inla_data_list(df = df,
                             n = n,
                             np = np,
                             outcome_list = outcome_list,
                             offset_list = offset_list,
                             id_list = id_list,
                             month_list = month_list,
                             extra_time_int = TRUE,
                             month_list_lag = NULL)
  
  
 # extract the data 
 y <- inla_lists$y
 int_list <- inla_lists$int_list
 id_list <- inla_lists$id_list
 month_list <- inla_lists$month_list
 offset_list <- inla_lists$offset_list
 time_list <- inla_lists$time_lists
 month_list_lag <- inla_lists$month_list_lag
  
  if(np == 2){
    
    link_set <- c(rep(1,nrow(df)),
                  rep(2,nrow(df)))
    
    
    # name the appropriate vectors
    
    ####################
    id1 <- id_list[[1]]
    id2 <- id_list[[2]]
    ###################
    
    #####################
    int1 <- int_list[[1]]
    int2 <- int_list[[2]]
    #####################
    
    ########################
    m1 <- month_list[[1]]
    m2 <- month_list[[2]]
    ########################
    
    #####################
    E1 <- offset_list[[1]]
    E2 <- offset_list[[2]]
    #####################
    
    ######################
    march_list <-  time_list[[1]]
    vac_list <- time_list[[2]]
    
    #######################
    mar1 <- march_list[[1]]
    mar2 <- march_list[[2]]
    #######################
    
    #######################
    vac1 <- vac_list[[1]]
    vac2 <- vac_list[[2]]
    #######################
    
    test_mod <- INLA::inla(y ~ -1+ int1 + int2+
                             mar1 + mar2 +
                             vac1 + vac2 +
                             f(m1,
                               model= "ar",
                               order = order,
                               adjust.for.con.comp = TRUE,
                               hyper = hyper.prec)+
                             f(m2,
                               copy = "m1",
                               #group = m2,
                               hyper = list(beta = beta.prior)),
                           family=c("poisson","zeroinflatedpoisson0"),
                           data=list(y =y,
                                     E1 = E1,
                                     E2 = E2,
                                     int1 = int1,
                                     int2 = int2,
                                     id1 = id1,
                                     id2 = id2,
                                     mar1 = mar1,
                                     mar2 = mar2,
                                     vac1 = vac1,
                                     vac2 = vac2),
                           control.compute=list(waic=TRUE),
                           control.predictor = list(compute = TRUE,
                                                    link = link_set),
                           control.inla=list(cmin=0),
                           control.fixed = list(prec = prec_fixed),
                           safe = TRUE,
                           verbose = verbose)
    
    
    predictions <- inla_predict(inla_mod = test_mod,
                                n = n,
                                np = np,
                                df = df)
    
    df <- predictions$df
    preds <- predictions$preds
    
    
    
  }else{
    
    
    link_set <- c(rep(1,nrow(df)),
                  rep(2,nrow(df)),
                  rep(3,nrow(df)))
  
    
    # name the appropriate vectors
    
    ####################
    id1 <- id_list[[1]]
    id2 <- id_list[[2]]
    id3 <- id_list[[3]]
    ###################
    
    #####################
    int1 <- int_list[[1]]
    int2 <- int_list[[2]]
    int3 <- int_list[[3]]
    #####################
    
    ########################
    m1 <- month_list[[1]]
    m2 <- month_list[[2]]
    m3 <- month_list[[3]]
    ########################
    
    #####################
    E1 <- offset_list[[1]]
    E2 <- offset_list[[2]]
    E3 <- offset_list[[3]]
    #####################
    
    ######################
   march_list <-  time_list[[1]]
   vac_list <- time_list[[2]]
   
   #######################
   mar1 <- march_list[[1]]
   mar2 <- march_list[[2]]
   mar3 <- march_list[[3]]
   #######################
   
   #######################
   vac1 <- vac_list[[1]]
   vac2 <- vac_list[[2]]
   vac3 <- vac_list[[3]]
   #######################
   
    test_mod <- INLA::inla(y ~ -1+ int1 + int2 +int3+
                             mar1 + mar2 + mar3 +
                             vac1 + vac2 + vac3+
                             f(m1,
                               model= "ar",
                               order = ar_order,
                               adjust.for.con.comp = TRUE,
                               hyper = hyper.prec)+
                             f(m2,
                               copy = "m1",
                               #group = m2,
                               hyper = list(beta = beta.prior))+
                             f(m3,
                               copy = "m1",
                               # group = m3,
                               hyper = list(beta = beta.prior)),
                           family=c("poisson","binomial","zeroinflatedpoisson0"),
                           data=list(y =y,
                                     E1 = E1,
                                     E2 = E2,
                                     E3 = E3,
                                     int1 = int1,
                                     int2 = int2,
                                     int3 = int3,
                                     id1 = id1,
                                     id2 = id2,
                                     id3 = id3,
                                     mar1 = mar1,
                                     mar2 = mar2,
                                     mar3 = mar3,
                                     vac1 = vac1,
                                     vac2 = vac2,
                                     vac3 = vac3),
                           control.compute=list(waic=TRUE),
                           control.predictor = list(compute = TRUE,
                                                    link = link_set),
                           control.inla=list(cmin=0),
                           control.fixed = list(prec = prec_fixed),
                           safe = TRUE,
                           verbose = verbose)
   
   
   predictions <- inla_predict(inla_mod = test_mod,
                               n = n,
                               np = np,
                               df = df)
   
   df <- predictions$df
   preds <- predictions$preds
   

    
  }
 
  }else{
    out <- list(df = df,
                preds = NULL,
                inla_mod = NULL)
  }
  
 out <- list(df = df,
             preds = preds,
             inla_mod = test_mod)
 
 out
  
}

# Inla run state function -------------------------



inla_run_state <- function(state_df,
                           state_abr = "NY",
                           spatial = T,
                           spatial_model = "'bym2'",
                           grouped_effect = T,
                           grouped_type = "'rw2'",
                           temporal_model = "no",
                           ar_order = 5,
                           time_intercepts = TRUE,
                           jhu_zeros = T,
                           beta.prior = beta.prior,
                           hyper.prec = hyper.prec,
                           prec_fixed = 1,
                           mar = T,
                           vac = T,
                           late = T,
                           run_inla = T,
                           rm_nas = T,
                           data_prep = T,
                           imputed_outcome = F,
                           verbose = T){
  
  
  
  usa_adj <- read.csv(here::here("usa_adj_mat_files/aug_adj_mat.csv"))
  
  
  state_adj_list <-  state_adj_fun(usa_adj = usa_adj,
                                   state_abr = state_abr,
                                   fips = T) 
  
  adj.mat <- state_adj_list$state_mat
  
  
  if(data_prep == T){

  
  # add remapped ids from 1:n_counties to the state data set ---------------
  state_df <- map_adj_ids(adj_list = state_adj_list,
                    state_df = state_df)
  
  
  # join in the population data --------------
  
  source(here::here("census-data/clean-pop-data.R")) # note this step may not work for all states, see source code for explanation
  
  state_df$population <- plyr::mapvalues(state_df$county_fips_code,
                                   from = pop$fips_num,
                                   to = pop$est_july_2020)
  
  state_df$population_scale <- state_df$population/(median(state_df$cases,na.rm = T)*50)
  
  

  
  
  if(rm_nas == T){
    na_df <- state_df %>%
      dplyr::group_by(case_month) %>%
      dplyr::summarise(na_cases = mean(is.na(cases)),
                       na_jhu = mean(is.na(jhu_cases)))
  nas <- which((na_df$na_cases == 1) & (na_df$na_jhu ==1))   
  
  state_df <- state_df %>% 
    dplyr::filter(!(case_month %in% na_df$case_month[nas]))
  
  }
  
  state_df <- map_month_ids(df = state_df,
                      beg = NULL, # default beg and end will be first and last month respectively
                      end = NULL)
  
  
  # set the dimensions of the flattened data for inla transformation

  
  # artificially set the negative counts to 0 (Change this later) 
  
  state_df$jhu_cases_reset <- ifelse(state_df$jhu_cases < 0, 0,state_df$jhu_cases)
  state_df$jhu_zero <- ifelse(state_df$jhu_cases_reset ==0,1,ifelse(is.na(state_df$jhu_cases),NA,0))
  
  state_df$jhu_cases <- ifelse(is.na(state_df$jhu_cases),NA,state_df$jhu_cases_reset)
  
  
  # order the dataframe to avoid NA issues in the predict step
  
  state_df  <- state_df %>% dplyr::mutate(case_na = ifelse(is.na(cases),1,0),
                              jhu_na = ifelse(is.na(jhu_cases),1,0))
  
  
  }
  state_df <- state_df %>% dplyr::arrange(case_na,
                              jhu_na,
                              county_fips_code)
  
  n <- nrow(state_df)
  
  #######################################3
  
  if(jhu_zeros == T){
    np <- 3
    
    if(imputed_outcome == F){
    outcome_list <- list(state_df$cases,
                         state_df$jhu_zero,
                         state_df$jhu_cases)
    }else{
      
      outcome_list <- list(state_df$cdc_impute,
                           state_df$jhu_zero,
                           state_df$jhu_impute)
    }
    
    offset_list <- rep_vector(np = np,
                              x = state_df$population_scale)
    
    id_list <- rep_vector(np = np,
                          x = state_df$remapped_id)
    
    month_list <- rep_vector(np = np,
                             x = state_df$month_id)
  }else{
    np <- 2
    if(imputed_outcome == F){
    outcome_list <- list(state_df$cases,
                         state_df$jhu_cases)
    }else{
      outcome_list <- list(state_df$cdc_impute,
                           state_df$jhu_impute)
    }
    
    offset_list <- rep_vector(np = np,
                              x = state_df$population_scale)
    
    id_list <- rep_vector(np = np,
                          x = state_df$remapped_id)
    
    month_list <- rep_vector(np = np,
                             x = state_df$month_id)
    
  }
  
  
  inla_lists <-  inla_data_list(df = state_df,
                                n = n,
                                np = np,
                                outcome_list = outcome_list,
                                offset_list = offset_list,
                                id_list = id_list,
                                month_list = month_list,
                                extra_time_int = TRUE)
  
  
  # extract the data 
  y <- inla_lists$y
  int_list <- inla_lists$int_list
  id_list <- inla_lists$id_list
  month_list <- inla_lists$month_list
  offset_list <- inla_lists$offset_list
  time_list <- inla_lists$time_lists
  
  
  if(np == 2){
    
    link_set <- c(rep(1,nrow(state_df)),
                  rep(2,nrow(state_df)))
    
    
    # name the appropriate vectors
    
    ####################
    id1 <- id_list[[1]]
    id2 <- id_list[[2]]
    ###################
    
    #####################
    int1 <- int_list[[1]]
    int2 <- int_list[[2]]
    #####################
    
    ########################
    m1 <- month_list[[1]]
    m2 <- month_list[[2]]
    ########################
    
    #####################
    E1 <- offset_list[[1]]
    E2 <- offset_list[[2]]
    #####################
    
    ######################
    march_list <-  time_list[[1]]
    vac_list <- time_list[[2]]
    late_list <- time_list[[3]]
    
    #######################
    mar1 <- march_list[[1]]
    mar2 <- march_list[[2]]
    #######################
    
    #######################
    vac1 <- vac_list[[1]]
    vac2 <- vac_list[[2]]
    #######################
    
    #######################
    late1 <- late_list[[1]]
    late2 <- late_list[[2]]
    #######################
    
   inla_form <- inla_formula(np = np,
                             spatial = spatial,
                             spatial_model = spatial_model,
                             grouped_effect = grouped_effect,
                             grouped_type = grouped_type,
                             temporal_model = temporal_model,
                             order = ar_order,
                             mar = mar,
                             vac = vac,
                             late = late) 
    
    
    if(run_inla == T){
    test_mod <- INLA::inla(as.formula(inla_form),
                           family=c("poisson","zeroinflatedpoisson0"),
                           data=list(y =y,
                                     E1 = E1,
                                     E2 = E2,
                                     int1 = int1,
                                     int2 = int2,
                                     id1 = id1,
                                     id2 = id2,
                                     mar1 = mar1,
                                     mar2 = mar2,
                                     vac1 = vac1,
                                     vac2 = vac2,
                                     late1 = late1,
                                     late2 = late2),
                           control.compute=list(waic=TRUE),
                           control.predictor = list(compute = TRUE,
                                                    link = link_set),
                           control.inla=list(cmin=0),
                           control.fixed = list(prec = prec_fixed),
                           safe = TRUE,
                           verbose = verbose)
    
    
    
    predictions <- inla_predict(inla_mod = test_mod,
                                n = n,
                                np = np,
                                df = df)
    
    
    df <- predictions$df
    preds <- predictions$preds
    
    }else{
      
      
      
      df <- state_df
      preds <- NULL
      test_mod <- NULL
    }
    
 
    
    
    
  }else{
    
    link_set <- c(rep(1,nrow(state_df)),
                  rep(2,nrow(state_df)),
                  rep(3,nrow(state_df)))
    
    
    # name the appropriate vectors
    
    ####################
    id1 <- id_list[[1]]
    id2 <- id_list[[2]]
    id3 <- id_list[[3]]
    ###################
    
    #####################
    int1 <- int_list[[1]]
    int2 <- int_list[[2]]
    int3 <- int_list[[3]]
    #####################
    
    ########################
    m1 <- month_list[[1]]
    m2 <- month_list[[2]]
    m3 <- month_list[[3]]
    ########################
    
    #####################
    E1 <- offset_list[[1]]
    E2 <- offset_list[[2]]
    E3 <- offset_list[[3]]
    #####################
    
    ######################
    march_list <-  time_list[[1]]
    vac_list <- time_list[[2]]
    late_list <- time_list[[3]]
    
    #######################
    mar1 <- march_list[[1]]
    mar2 <- march_list[[2]]
    mar3 <- march_list[[3]]
    #######################
    
    #######################
    vac1 <- vac_list[[1]]
    vac2 <- vac_list[[2]]
    vac3 <- vac_list[[3]]
    #######################
    
    #######################
    late1 <- late_list[[1]]
    late2 <- late_list[[2]]
    late3 <- late_list[[3]]
    #######################
    
    inla_form <- inla_formula(np = np,
                              spatial = spatial,
                              spatial_model = spatial_model,
                              grouped_effect = grouped_effect,
                              grouped_type = grouped_type,
                              temporal_model = temporal_model,
                              order = ar_order,
                              mar = mar,
                              vac = vac,
                              late = late) 
    
    #inla_form
    
    if(run_inla == TRUE){
    test_mod <- INLA::inla(as.formula(inla_form),
                           family=c("poisson","binomial","zeroinflatedpoisson0"),
                           data=list(y =as.matrix(y),
                                     E1 = E1,
                                     E2 = E2,
                                     E3 = E3,
                                     int1 = int1,
                                     int2 = int2,
                                     int3 = int3,
                                     id1 = id1,
                                     id2 = id2,
                                     id3 = id3,
                                     mar1 = mar1,
                                     mar2 = mar2,
                                     mar3 = mar3,
                                     vac1 = vac1,
                                     vac2 = vac2,
                                     vac3 = vac3,
                                     late1 = late1,
                                     late2 = late2,
                                     late3 = late3,
                                     m1 = m1,
                                     m2 = m2,
                                     m3 = m3,
                                     adj.mat = as.matrix(adj.mat)),
                           control.compute=list(waic=TRUE),
                           control.predictor = list(compute = TRUE,
                                                    link = link_set),
                           control.inla=list(cmin=0),
                           control.fixed = list(prec = prec_fixed),
                           safe = TRUE,
                           verbose = verbose)
    
    
    predictions <- inla_predict(inla_mod = test_mod,
                                n = n,
                                np = np,
                                df = state_df)
    
    df <- predictions$df
    preds <- predictions$preds
    
    }else{
      
      
      
      df <- state_df
      preds <- NULL
      test_mod <- NULL
    }
    
  }
  
  
out <- list(df = df,
            preds = preds)  
  
 out 
}