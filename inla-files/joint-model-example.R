# Here is an example using the new joint modelling functions on the NY data set
# By joint, I mean modelling of cases, hospitalizations, and counts simultaneously

# Housekeeping -----------------------

# These are special functions which must be loaded (unless you prefer to load dplyr and lubridate full libraries)
`%>%` <- dplyr::`%>%`
`%m+%` <- lubridate::`%m+%`

# Source important files with functions

source(here::here("inla-files/inla-joint-model-functions.R"))

# Load case data ---------

ny <- read.csv(here::here("state-data/NY_COVID_data.csv"))

# create adj matrix ----------


source(here::here("usa_adj_mat_files/adj-mat-functions.R"))


usa_adj <- read.csv(here::here("usa_adj_mat_files/aug_adj_mat.csv"))


ny_adj_list <-  state_adj_fun(usa_adj = usa_adj,
                              state_abr = "NY",
                              fips = T) 

adj.mat <- ny_adj_list$state_mat

# add remapped ids from 1:n_counties to the state data set ---------------
ny <- map_adj_ids(adj_list = ny_adj_list,
                  state_df = ny)


# join in the population data --------------

source(here::here("census-data/clean-pop-data.R")) # note this step may not work for all states, see source code for explanation

ny$population <- plyr::mapvalues(ny$county_fips_code,
                                 from = pop$fips_num,
                                 to = pop$est_july_2020)

scale_factor <- 500 # NOTE: This may have to be changed in some cases, especially state to state due to potential overflow issues

ny$population_scale <- ny$population/(median(ny$cases,na.rm = T)*scale_factor)


# Convert the dates to numeric 1:n_months --------

source(here::here("utilities/date-functions.R")) # load mapping function


#ny <- ny %>% dplyr::filter(!(case_month %in% c("2020-01-01",
# "2020-02-01")))

na_df <- ny %>% dplyr::group_by(case_month) %>% dplyr::summarise(na_cases = mean(is.na(cases)),
                                                                 na_jhu = mean(is.na(jhu_cases)))


nas <- which((na_df$na_cases == 1) & (na_df$na_jhu ==1) )   
ny <- ny %>% dplyr::filter(!(case_month %in% na_df$case_month[nas]))
ny <- map_month_ids(df = ny,
                    beg = NULL, # default beg and end will be first and last month respectively
                    end = NULL)


# artificially set the negative counts to 0 (Change this later) 

ny$jhu_cases_reset <- ifelse(ny$jhu_cases < 0, 0,ny$jhu_cases)
ny$jhu_zero <- ifelse(ny$jhu_cases_reset ==0,1,ifelse(is.na(ny$jhu_cases),NA,0))

ny$jhu_cases <- ifelse(is.na(ny$jhu_cases),NA,ny$jhu_cases_reset)




# order the dataframe to avoid NA issues in the predict step

ny  <- ny %>% dplyr::mutate(case_na = ifelse(is.na(cases),1,0),
                            jhu_na = ifelse(is.na(jhu_cases),1,0))

ny <- ny %>% dplyr::arrange(case_na,
                            jhu_na,
                            county_fips_code)


ny <- hosp_outcome_create(ny) # This creates additional variables
ny <- death_outcome_create(ny)

# Create the state dataframe -----------

state_df <- ny

#county_fips <- 36005

#df <- state_df%>% dplyr::filter(county_fips_code == county_fips)

# Decide the names of the columns you'd like to model and the corresponding likelihoods for each of the outcome types
# Note: This should be set up such that the modelled columns can change state to state or include new columns not listed here

## Cases
cases_outcomes = c("cases","jhu_zero", "jhu_cases")
cases_likelihoods = c("poisson","binomial","zeroinflatedpoisson0")
cases_group <- c(1,2,2) # Notice that jhu_zero and jhu_cases are really the same column but we are modelling the zeros and counts separately, this constitutes a group

# hopsitalizations
hosp_outcomes = c("hosp_yes","hosp_zero","new_hosp_total","new_hosp_zero")
hosp_likelihoods = c("zeroinflatedpoisson0","binomial","zeroinflatedpoisson0","binomial")
hosp_group <- c(1,1,2,2) # note there are two full groups here

# Deaths 
deaths_outcomes = c("death_yes","death_zero","CovidDeathCount","covid_death_zero","jhu_deaths","jhu_deaths_zero")
deaths_likelihoods = c("zeroinflatedpoisson0","binomial","zeroinflatedpoisson0","binomial","zeroinflatedpoisson0","binomial")
deaths_group <- c(1,1,2,2,3,3)


# Run the model -------------------------------------

## Set model parameters

spatial_model <- 'bym2' # which model type
time_intercepts <- T # should we calculate intercepts for important time points in the pandemic
verbose <- T # should the inla model return a verbose output
hospitalizations <- T

lag_hosp <- F #  Should we included an additional lagged spatial effect for hospitalizations?
lag_deaths <- F # Should we included an additional lagged spatial effect for deaths?

ind_intercept <- T # should each column get it's own intercept (as opposed to a within type-likelihood re only)?
mar <- T # include intercept for march 2020 and onwards
vac <- F # include intercept for rough demarcation of the availability of the vacine
ar_order <- 5 # what order should the ar process be?


model_output <- inla_run_joint(df,
                               spatial_model = spatial_model,
                               time_intercepts = time_intercepts,
                               beta.prior = beta.prior,
                               hyper.prec = hyper.prec,
                               prec_fixed = 1,
                               run = T,
                               verbose = verbose,
                               cases_outcomes = cases_outcomes,
                               cases_likelihoods = cases_likelihoods,
                               hosp_outcomes = hosp_outcomes,
                               hosp_likelihoods = hosp_likelihoods,
                               deaths_outcomes = deaths_outcomes,
                               deaths_likelihoods = deaths_likelihoods,
                               cases_group = cases_group,
                               hosp_group = hosp_group,
                               deaths_group = deaths_group,
                               lag_hosp = lag_hosp,
                               lag_deaths = lag_deaths,
                               ind_intercept = ind_intercept,
                               mar = mar,
                               vac = vac,
                               ar_order = ar_order)

