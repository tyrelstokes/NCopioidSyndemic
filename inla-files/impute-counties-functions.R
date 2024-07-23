# Housekeeping -------
`%>%` <- dplyr::`%>%`
`%do%` <- foreach::`%do%`

source(here::here("inla-files/inla_state_run_functions.R"))

# impute county function -----------
county_impute <- function(state_df,
                          fips,
                          ar_order = 5,
                          time_intercepts = TRUE,
                          jhu_zeros = T,
                          beta.prior = beta.prior,
                          hyper.prec = hyper.prec,
                          prec_fixed = 1,
                          impute_month_vec = c(38:48),
                          outcomes = c("jhu_cases"),
                          run = TRUE){
  
 
  
  county_model <-   inla_run_county(state_df = state_df,
                                county_fips = fips,
                                ar_order = ar_order,
                                time_intercepts = time_intercepts,
                                jhu_zeros = jhu_zeros,
                                beta.prior = beta.prior,
                                hyper.prec = hyper.prec,
                                prec_fixed = prec_fixed,
                                run = run) 
  
  
  county_df <- county_model$df
  
 
  month_nas <- county_df$month_id[which(is.na(county_df$jhu_cases))]
  
  imputed_months <- month_nas[month_nas %in% impute_month_vec]
  
  cdc <- county_df$cases
  jhu <- county_df$jhu_cases
  
  if("cases" %in% outcomes){
    county_df$cdc_impute <- ifelse(county_df$month_id %in% imputed_months,county_df$pred_cdc,county_df$cases)
  }else{
    county_df$cdc_impute <- county_df$cases
  }
  
 
  if("jhu_cases" %in% outcomes){
    if(jhu_zeros == T){
    county_df$jhu_impute <- ifelse(county_df$month_id %in% imputed_months,county_df$jhu_comb,county_df$jhu_cases)
    }else{
      county_df$jhu_impute <- ifelse(county_df$month_id %in% imputed_months,county_df$pred_jhu,county_df$jhu_cases)
      
    }
  }else{
    county_df$jhu_impute <- county_df$cases
  }
  
 county_df 
   
}


# impute many counties --------------


county_impute_many <- function(state_df,
                          fips_vec,
                          ar_order = 5,
                          time_intercepts = TRUE,
                          jhu_zeros = T,
                          beta.prior = beta.prior,
                          hyper.prec = hyper.prec,
                          prec_fixed = 1,
                          impute_month_vec = c(38:48),
                          outcomes = c("jhu_cases")){
  
 
  n_counties <- length(fips_vec)
  
  counties <- foreach::foreach(i = 1:n_counties,.combine = rbind)%do%{
    
    county_impute(state_df = state_df,
                  fips = fips_vec[i],
                  ar_order = ar_order,
                  time_intercepts = time_intercepts,
                  jhu_zeros = jhu_zeros,
                  beta.prior = beta.prior,
                  hyper.prec = hyper.prec,
                  prec_fixed = prec_fixed,
                  impute_month_vec = impute_month_vec,
                  outcomes = outcomes,
                  run = TRUE) 
    
    
  }
  
  counties <- counties %>% dplyr::arrange(county_fips_code,month_id)
  state_df <- state_df %>% dplyr::arrange(county_fips_code,month_id)
  
  ind <- which(state_df$county_fips_code %in% fips_vec)
  
  state_df$cdc_impute <- state_df$cases
  state_df$jhu_impute <- state_df$jhu_cases
  
  state_df$cdc_impute[ind] <- counties$cdc_impute
  state_df$jhu_impute[ind] <- counties$jhu_impute
 
 
  state_df
  
}


# Run state with imputation --------------------

run_state_with_impute <- function(state_df,
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
                                  impute_fips_vec,
                                  ar_order = 5,
                                  time_intercepts = TRUE,
                                  jhu_zeros = T,
                                  impute_month_vec = c(38:48),
                                  outcomes = c("jhu_cases","cases")){
  
  
  
  state_list <- inla_run_state(state_df = state_df,
                               state_abr = state_abr,
                               spatial = spatial,
                               spatial_model = spatial_model,
                               grouped_effect = grouped_effect,
                               grouped_type = grouped_type,
                               temporal_model = temporal_model,
                               ar_order = ar_order,
                               time_intercepts = time_intercepts,
                               jhu_zeros = jhu_zeros,
                               beta.prior = beta.prior,
                               hyper.prec = hyper.prec,
                               prec_fixed = prec_fixed,
                               mar = mar,
                               vac = vac,
                               late = late,
                               run_inla = F)
  
  state_df <- state_list$df
  
state_df_impute <-   county_impute_many(state_df = state_df,
                     fips_vec = fips_vec,
                     ar_order = ar_order,
                     time_intercepts = time_intercepts,
                     jhu_zeros = jhu_zeros,
                     beta.prior = beta.prior,
                     hyper.prec = hyper.prec,
                     prec_fixed = prec_fixed,
                     impute_month_vec = impute_month_vec,
                     outcomes = outcomes) 



out_list <- inla_run_state(state_df = state_df_impute,
                             state_abr = state_abr,
                             spatial = spatial,
                             spatial_model = spatial_model,
                             grouped_effect = grouped_effect,
                             grouped_type = grouped_type,
                             temporal_model = temporal_model,
                             ar_order = ar_order,
                             time_intercepts = time_intercepts,
                             jhu_zeros = jhu_zeros,
                             beta.prior = beta.prior,
                             hyper.prec = hyper.prec,
                             prec_fixed = prec_fixed,
                             mar = mar,
                             vac = vac,
                             late = late,
                             run_inla = T,
                             data_prep = F) 
  
}

  
