
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(cmdstanr)
set_cmdstan_path(path = "C:/Users/mckayc04/AppData/Local/R/win-library/4.4/cmdstanr/")
set_cmdstan_path(path = "C:/Users/mckayc04/.cmdstan/cmdstan-2.35.0")
set_cmdstan_path(path = NULL)
library(loo)

#set_cmdstan_path()

#path_dir <- "C:/Users/mckayc04/Documents/NCopioidSyndemic-explore/"
path_dir <- "C:/Users/mckayc04/Documents/growth-files/"
setwd(path_dir)

#source(paste0(path_dir, "growth-files/stan-fit-functions.R"))
source(paste0(path_dir, "stan-fit-functions.R"))


############


nflis <- read.csv(paste0(path_dir, "nflis_data_RBeast.csv"))

nflis_region1 <- nflis %>% dplyr::select(Connecticut, Maine, Massachusetts, New.Hampshire, Rhode.Island,
                                         Vermont)
nflis_region2 <- nflis %>% dplyr::select(New.Jersey, New.York)
nflis_region3 <- nflis %>% dplyr::select(Delaware, District.of.Columbia, Maryland, Pennsylvania, Virginia,
                                         West.Virginia)
nflis_region4 <- nflis %>% dplyr::select(Alabama, Florida, Georgia, Kentucky, Mississippi, North.Carolina,
                                         South.Carolina, Tennessee)
nflis_region5 <- nflis %>% dplyr::select(Illinois, Indiana, Michigan, Minnesota, Ohio, Wisconsin)
nflis_region6 <- nflis %>% dplyr::select(Arkansas, Louisiana, New.Mexico, Oklahoma, Texas)
nflis_region7 <- nflis %>% dplyr::select(Iowa, Kansas, Missouri, Nebraska)
nflis_region8 <- nflis %>% dplyr::select(Colorado, Montana, North.Dakota, South.Dakota, Utah, Wyoming)
nflis_region9 <- nflis %>% dplyr::select(Arizona, California, Nevada)
nflis_region10 <- nflis %>% dplyr::select(Idaho, Oregon, Washington)


df <-   do.call(rbind,lapply(c(1:ncol(nflis_region9)),
                             function(i){
                               y <- nflis_region9[,i]
                               t <- c(1:nrow(nflis_region9))
                               state_id <- rep(i,nrow(nflis_region9))
                               
                               out <- data.frame(y = y,t = t, state_id = state_id)
                               out
                             }))



n <- nrow(nflis_region9)
n_s <- ncol(nflis_region9)

###################


c_ints <- c(1:n_s)


# Single hurdle -----------------------------------------------------------

start_time <- Sys.time()


p_delta_h_sin <- c()
pred_list_h_sin <- c()

for (i in 1:n) {
  
  h_sin <- run_many(n = n, # number of observations per state
                    df = df, # data frame will all data
                    t_cutoff = i, # year number we consider the 'covid effect' >=
                    c_ints = c_ints, # state ids we would like to run on
                    eps = 0.0001, # threshold below which we consider growth off!
                    model = "hurdle", # hurdle logistic model
                    fill = 0.0001, # what number to fill in zeros with
                    draws = F, # do want the posterior draws to be stored and outputted 
                    snames = colnames(nflis_region2)[c_ints], # vector of names for the states (for plots only right now)
                    prob_fill = T, # in the default plot should the histogram be filled with the p(delta >0) as opposed to the mean delta value
                    waic = T) # report WAIC for each run
  
  h_sin$pdelta$t_cutoff <- i
  
  p_delta_h_sin <- rbind(p_delta_h_sin, h_sin$pdelta)
  pred_list_h_sin <- c(pred_list_h_sin, h_sin$pred_list)
  
  
  
}

write.csv(p_delta_h_sin,here::here("hurdle_single_region2.csv"))


end_time <- Sys.time()

# Calculate the time difference
execution_time <- end_time - start_time

# Single logistic ---------------------------------------------------------

start_time <- Sys.time()

p_delta_g_sin <- c()
pred_list_g_sin <- c()

for (i in 1:n) {
  
  g_sin <- run_many(n = n,
                    df = df,
                    t_cutoff = i,
                    c_ints = c_ints,
                    eps = 0.001,
                    model = "logistic",
                    fill = 0.0001,
                    draws = F,
                    snames = colnames(nflis_region2)[c_ints],
                    prob_fill = T,
                    waic = T)
  
  g_sin$pdelta$t_cutoff <- i
  
  p_delta_g_sin <- rbind(p_delta_g_sin, g_sin$pdelta)
  pred_list_g_sin <- c(pred_list_g_sin, g_sin$pred_list)
  
  
  
}

write.csv(p_delta_g_sin,here::here("logistic_single_region2.csv"))


end_time <- Sys.time()

# Calculate the time difference
execution_time <- end_time - start_time

# Hier hurdle -------------------------------------------------------------

start_time <- Sys.time()

p_delta_h_heir <- c()
#pred_list_h_hier <- c()

for (i in 1:n) {
  
  h_hier <- run_heir_alt(model = "hurdle",
                         fill = 0.001,
                         eps = 0.0001,
                         iter_warmup = 1250, # how many posterior warmup samples
                         iter_sampling = 1250, # how many posterior actual stored samples
                         predictions = 0, # 1 = output predictions, otherwise predictions not calculated
                         t_cutoff = i,
                         sname_list = c("Arizona", "California", "Nevada"),
                         waic = T)
  
  # h_hier <- run_heir(model = "hurdle",
  #                    fill = 0.001,
  #                    eps = 0.0001,
  #                    iter_warmup = 1250, # how many posterior warmup samples
  #                    iter_sampling = 1250, # how many posterior actual stored samples
  #                    predictions = 1) # 1 = output predictions, otherwise predictions not calculated
  
  
  p_delta_h_heir <- rbind(p_delta_h_heir, h_hier[[1]])
  #pred_list_h_hier <- c(pred_list_h_hier, h_hier$pred_list)
  
}

write.csv(p_delta_h_heir,here::here("hurdle_hier_region9.csv"))


end_time <- Sys.time()

# Calculate the time difference
execution_time <- end_time - start_time

# Hier logistic -----------------------------------------------------------

start_time <- Sys.time()

p_delta_l_heir <- c()
#pred_list_l_hier <- c()

for (i in 1:n) {
  
  l_hier <- run_heir_alt(model = "logistic",
                         fill = 0.001,
                         eps = 0.0001,
                         iter_warmup = 1250,
                         iter_sampling = 1250,
                         predictions = 0, # 1 = output predictions, otherwise predictions not calculated
                         t_cutoff = i,
                         sname_list = c("Arizona", "California", "Nevada"),
                         waic = T)
  
  # l_hier <- run_heir(model = "logistic",
  #                    fill = 0.001,
  #                    eps = 0.0001,
  #                    iter_warmup = 1250,
  #                    iter_sampling = 1250,
  #                    predictions = 1)
  
  p_delta_l_heir <- rbind(p_delta_l_heir, l_hier[[1]])
  #pred_list_l_hier <- c(pred_list_l_hier, l_hier$pred_list)
  
}

write.csv(p_delta_l_heir,here::here("logistic_hier_region9.csv"))

end_time <- Sys.time()

# Calculate the time difference
execution_time <- end_time - start_time



# Hier exp hurdle ---------------------------------------------------------

start_time <- Sys.time()

p_delta_e_hier_1 <- c()
#pred_list_e_hier_1 <- c()

for (i in 1:n) {
  
  e_hier_1 <- run_heir_alt(model = "exp-hurdle",
                           fill = 0.001,
                           eps = 0.0001,
                           iter_warmup = 1250,
                           iter_sampling = 1250,
                           predictions = 0, # 1 = output predictions, otherwise predictions not calculated
                           t_cutoff = i,
                           sname_list = c("Arizona", "California", "Nevada"),
                           waic = T)
  
  # e_hier_1 <- run_heir(model = "exp-hurdle",
  #                      fill = 0.001,
  #                      eps = 0.0001,
  #                      iter_warmup = 1250,
  #                      iter_sampling = 1250,
  #                      predictions = 1)
  
  p_delta_e_hier_1 <- rbind(p_delta_e_hier_1, e_hier_1[[1]])
  #pred_list_e_hier_1 <- c(pred_list_e_hier_1, e_hier_1$pred_list)
  
}

write.csv(p_delta_e_hier_1,here::here("exp_hurdle_hier_region9.csv"))

end_time <- Sys.time()

# Calculate the time difference
execution_time <- end_time - start_time



# Hier exp ----------------------------------------------------------------

start_time <- Sys.time()

p_delta_e_hier_2 <- c()
#pred_list_e_hier_2 <- c()

for (i in 1:n) {
  
  e_hier_2 <- run_heir_alt(model = "exp",
                           fill = 0.001,
                           eps = 0.0001,
                           iter_warmup = 1250,
                           iter_sampling = 1250,
                           predictions = 0, # 1 = output predictions, otherwise predictions not calculated
                           t_cutoff = i,
                           sname_list = c("Arizona", "California", "Nevada"),
                           waic = T)
  
  # e_hier_2 <- run_heir(model = "exp",
  #                      fill = 0.001,
  #                      eps = 0.0001,
  #                      iter_warmup = 1250,
  #                      iter_sampling = 1250,
  #                      predictions = 1)
  
  
  p_delta_e_hier_2 <- rbind(p_delta_e_hier_2, e_hier_2[[1]])
  #pred_list_e_hier_2 <- c(pred_list_e_hier_2, e_hier_2$pred_list)
  
}

write.csv(p_delta_e_hier_2,here::here("exp_hier_region9.csv"))

end_time <- Sys.time()

# Calculate the time difference
execution_time <- end_time - start_time




# Compare log likelihood --------------------------------------------------


loglik_df <- rbind(p_delta_h_heir,
                   p_delta_l_heir,
                   p_delta_e_hier_1,
                   p_delta_e_hier_2)


summary_df9 <- loglik_df %>%
  dplyr::group_by(model) %>%
  dplyr::summarize(
    mean_elpd_loo = mean(elpd_loo, na.rm = TRUE),
    sd_elpd_loo = sd(elpd_loo, na.rm = TRUE),
    mean_p_loo = mean(p_loo, na.rm = TRUE),
    sd_p_loo = sd(p_loo, na.rm = TRUE),
    mean_looic = mean(looic, na.rm = TRUE),
    sd_looic = sd(looic, na.rm = TRUE),
    .groups = "drop"
  )





# Extension ideas
# It might make more sense to have the information sharing operate with respect to an adjacency matrix
# or perhaps just regionally which would be easier to implement.
# look into the regional classification!!! That would probably work.

