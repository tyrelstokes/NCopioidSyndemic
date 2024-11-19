# test state functions ---------------

source(here::here("inla-files/inla_state_run_functions.R"))
source(here::here("inla-files/impute-counties-functions.R"))
source(here::here("utilities/plot-functions.R"))
ny <- read.csv(here::here("state-data/NY_COVID_data.csv"))


# find out which counties have lot's of missingness

ms_df <- ny %>% dplyr::group_by(county_fips_code) %>%
  dplyr::summarise(miss_cases = mean(is.na(cases)),
                   miss_jhu = mean(is.na(jhu_cases)))



impute_vec <- unique(ny$county_fips_code)
# removing that one weird county 36097?


out <- run_state_with_impute(state_df = ny,
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
                             impute_fips_vec = impute_vec,
                             impute_month_vec = c(38:48),
                             outcomes = c("jhu_cases","cases"),
                             verbose = F)


pred_df <- out$df


# fips: 36097 has tonnes of missing data and tends to behave poorly, it's not clear if imputation
# alone is helpful because one of the data sources is often missing. might have to think of a 
# strategy for this county

plot_ny <- plot_state(state_df = pred_df,
                      random_counties = T,
                      n_sample = 3,
                      which_counties = NULL,
                      ipsum = T,
                      imputed = T)

plot_ny + ggplot2::ggtitle("Random Sample of Counties - Predictions vs Actual")



#####################################################################
#####################################################################

# NEVADA

#####################3
######################


ne <- read.csv(here::here("state-data/NE_COVID_data.csv"))

ne <- ne %>% dplyr::filter(is.na(county_fips_code) == FALSE)

# find out which counties have lot's of missingness

ms_df <- ne %>% dplyr::group_by(county_fips_code) %>%
  dplyr::summarise(miss_cases = mean(is.na(cases)),
                   miss_jhu = mean(is.na(jhu_cases)))


ms_df2 <-  ne %>% dplyr::group_by(case_month) %>%
  dplyr::summarise(miss_cases = mean(is.na(cases)),
                   miss_jhu = mean(is.na(jhu_cases)))



impute_vec <- unique(ne$county_fips_code)
#impute_vec <- impute_vec[-length(impute_vec)]



out_ne <- run_state_with_impute(state_df = ne,
                             state_abr = "NE",
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
                             impute_fips_vec = impute_vec,
                             impute_month_vec = c(3:39),
                             outcomes = c("jhu_cases","cases"),
                             verbose = F)


pred_df_ne <- out_ne$df


plot_ne <- plot_state(state_df = pred_df_ne,
                      random_counties = T,
                      n_sample = 3,
                      which_counties = NULL,
                      ipsum = T,
                      imputed = T)

plot_ne




write.csv(pred_df_ne, 
          here::here("inla-files/model-results/ne_results.csv"))


write.csv(pred_df, 
          here::here("inla-files/model-results/ny_results.csv"))




