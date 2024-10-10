# Housekeeping -----------

`%>%` <- dplyr::`%>%`
`%m+%` <- lubridate::`%m+%`

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

ny$population_scale <- ny$population/(median(ny$cases,na.rm = T)*500)


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


ny <- hosp_outcome_create(ny)
ny <- death_outcome_create(ny)

# death outcomes create -----------






state_df <- ny

county_fips <- 36005

df <- state_df%>% dplyr::filter(county_fips_code == county_fips)


cases_outcomes = c("cases","jhu_zero", "jhu_cases")
cases_likelihoods = c("poisson","binomial","zeroinflatedpoisson0")
hosp_outcomes = c("hosp_yes","hosp_zero","new_hosp_total","new_hosp_zero")
hosp_likelihoods = c("zeroinflatedpoisson0","binomial","zeroinflatedpoisson0","binomial")
deaths_outcomes = c("death_yes","death_zero","CovidDeathCount","covid_death_zero","jhu_deaths","jhu_deaths_zero")
deaths_likelihoods = c("zeroinflatedpoisson0","binomial","zeroinflatedpoisson0","binomial","zeroinflatedpoisson0","binomial")

cases_group <- c(1,2,2)
hosp_group <- c(1,1,2,2)
deaths_group <- c(1,1,2,2,3,3)






