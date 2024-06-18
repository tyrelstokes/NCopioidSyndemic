source(here::here("usa_adj_mat_files/adj-mat-functions.R"))


usa_adj <- read.csv(here::here("usa_adj_mat_files/countyadj.csv"))


ny_adj_matrix <- state_adj_fun(usa_adj = usa_adj,
                               state_abr = "NY")


neb_adj_matrix <-  state_adj_fun(usa_adj = usa_adj,
                                 state_abr = "NE") 
