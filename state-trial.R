# test state functions ---------------

source(here::here("inla-files/inla_state_run_functions.R"))
source(here::here("inla-files/impute-counties-functions.R"))
ny <- read.csv(here::here("state-data/NY_COVID_data.csv"))


# find out which counties have lot's of missingness

ms_df <- ny %>% dplyr::group_by(county_fips_code) %>%
  dplyr::summarise(miss_cases = mean(is.na(cases)),
                   miss_jhu = mean(is.na(jhu_cases)))



impute_vec <- unique(ny$county_fips_code)
impute_vec <- impute_vec[-62]


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


out$preds
