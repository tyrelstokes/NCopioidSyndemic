# housekeeping -------

`%>%` <- dplyr::`%>%`

if(!require("urbnmapr")){
devtools::install_github("UrbanInstitute/urbnmapr")
}

# get data for plotting --------

get_county_data <- function(state_name){
  dt <- urbnmapr::counties %>% 
    dplyr::filter(state_name == state_name)
}

# plot inla model ------------

plot_inla_model <- function()

