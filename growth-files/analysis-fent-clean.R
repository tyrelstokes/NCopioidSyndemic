source(here::here("growth-files/stan-fit-functions.R"))


############


nflis <- read.csv(here::here("growth-files/nflis_data_RBeast.csv"))


df <-   do.call(rbind,lapply(c(1:ncol(nflis)),
                             function(i){
                               y <- nflis[,i]
                               t <- c(1:nrow(nflis))
                               state_id <- rep(i,nrow(nflis))
                               
                               out <- data.frame(y = y,t = t, state_id = state_id)
                               out
                             }))



n <- nrow(nflis)
n_s <- ncol(nflis)

###################


c_ints <- c(1:51)


h_sin <- run_many(n = n, # number of observations per state
                  df = df, # data frame will all data
                  t_cutoff = 14, # year number we consider the 'covid effect' >=
                  c_ints = c_ints, # state ids we would like to run on
                  eps = 0.0001, # threshold below which we consider growth off!
                  model = "hurdle", # hurdle logistic model
                  fill = 0.0001, # what number to fill in zeros with
                  draws = T, # do want the posterior draws to be stored and outputted 
                  snames = colnames(nflis)[c_ints], # vector of names for the states (for plots only right now)
                  prob_fill = T) # in the default plot should the histogram be filled with the p(delta >0) as opposed to the mean delta value

 write.csv(h_sin,here::here("growth-files/hurdle_single.csv"))

g_sin <- run_many(n = n,
                  df = df,
                  t_cutoff = 14,
                  c_ints = c_ints,
                  eps = 0.001,
                  model = "logistic",
                  fill = 0.0001,
                  draws = T,
                  snames = colnames(nflis)[c_ints],
                  prob_fill = T)


write.csv(g_sin,here::here("growth-files/logistic_single.csv"))


h_hier <- run_heir(model = "hurdle",
                   fill = 0.001,
                   eps = 0.0001,
                   iter_warmup = 1250, # how many posterior warmup samples
                   iter_sampling = 1250, # how many posterior actual stored samples
                   predictions = 1) # 1 = output predictions, otherwise predictions not calculated

write.csv(h_hier,here::here("growth-files/hurdle_hier.csv"))


l_hier <- run_heir(model = "logistic",
                   fill = 0.001,
                   eps = 0.0001,
                   iter_warmup = 1250,
                   iter_sampling = 1250,
                   predictions = 1)

write.csv(l_hier,here::here("growth-files/logistic_hier.csv"))

e_hier_1 <- run_heir(model = "exp-hurdle",
                   fill = 0.001,
                   eps = 0.0001,
                   iter_warmup = 1250,
                   iter_sampling = 1250,
                   predictions = 1)

write.csv(e_hier_1,here::here("growth-files/exp_hurdle_hier.csv"))


e_hier_2 <- run_heir(model = "exp",
                     fill = 0.001,
                     eps = 0.0001,
                     iter_warmup = 1250,
                     iter_sampling = 1250,
                     predictions = 1)

write.csv(e_hier_2,here::here("growth-files/exp_hier.csv"))


# Extension ideas
# It might make more sense to have the information sharing operate with respect to an adjacency matrix
# or perhaps just regionally which would be easier to implement.
# look into the regional classification!!! That would probably work.

