
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

stan_list <-  list(N = n*n_s,
                   y = ifelse(df$y ==0, 0.001,df$y),
                   t = df$t,
                   t_cutoff = 14,
                   state = df$state_id,
                   n_s = n_s,
                   eps = 0.001)


smod <- cmdstanr::cmdstan_model(here::here("growth-files/growth-model-states-gamma-simple.stan"))

smod <- cmdstanr::cmdstan_model(here::here("growth-files/growth-hurdle-many.stan"))


fit <- smod$sample(data = stan_list,
                   iter_warmup = 2000,
                   iter_sampling = 2000)

sumfit <-  fit$summary()



run_single <- function(n,
                       j,
                       df,
                       t_cutoff,
                       eps = 0.001,
                       fill = 0.0001,
                       model = "logistic"){
  
 df_s <- df %>% dplyr::filter(state_id ==j) 
 
 if(model %in% c("hurdle","exp-hurdle")){
   hurdle <- T
 }else{
   hurdle <- F
 }
 
 if(hurdle == F){
   y <- ifelse(df_s$y ==0,fill,df_s$y)
 }else{
   y <- df_s$y
 }
 

 
 sl <- list(y = y,
            N = n,
            t = df_s$t,
            t_cutoff = 14,
            eps = eps)
 
 if(model == "hurdle"){
 smods <- cmdstanr::cmdstan_model(here::here("growth-files/growth-hurdle-single.stan"))
 }
 
if(model == "logistic"){
   smods <- cmdstanr::cmdstan_model(here::here("growth-files/growth-gamma-single.stan")) 
}
 
 if(model == "exp"){
   smods <- cmdstanr::cmdstan_model(here::here("growth-files/growth-exp-single.stan")) 
 }
 
 if(model == "exp-hurdle"){
   smods <- cmdstanr::cmdstan_model(here::here("growth-files/growth-exp-single-hurdle.stan")) 
 }
  
 fit <- smods$sample(data = sl,
                    iter_warmup = 2000,
                    iter_sampling = 2000)
 
 sumfit <-  fit$summary()
 
 if(model == "hurdle"){
 param_df <- sumfit[1:10,]
 pred_df <- sumfit[11:nrow(sumfit),]
 }
 
 if(model == "logistic"){
   param_df <- sumfit[1:8,]
   pred_df <- sumfit[9:nrow(sumfit),]  
 }
 
 if(model == "exp"){
   param_df <- sumfit[1:4,]
   pred_df <- sumfit[5:nrow(sumfit),]   
 }
 
 if(model =="exp-hurdle"){
   param_df <- sumfit[1:6,]
   pred_df <- sumfit[7:nrow(sumfit),]   
 }
 
 pdelta <- param_df %>% dplyr::filter(variable =="p_delta")
 if(hurdle == T){
 pdelta$model <- "hurdle-single"
 }else{
   pdelta$model <- "gamma-single"
 }
 
 df_pred <- data.frame(y = sl$y,yhat = pred_df$mean)
 
 r2 = sqrt(mean((df_pred$y - df_pred$yhat)^2))
 
 out <- list(r2 = r2,
             param_df = param_df,
             df_pred = df_pred,
             pdelta = pdelta)
 
 out
}


run_many <- function(n,
                     df,
                     t_cutoff,
                     c_ints,
                     eps = 0.001,
                     hurdle = F){


out <- foreach::foreach(i = c_ints,.combine = rbind)%do%{
inter <- run_single(n = n,
                    j = i,
                    d = df,
                    t_cutoff = t_cutoff,
                    eps = eps,
                    hurdle = hurdle)

pdelta <- inter$pdelta
pdelta$state_id <- i
#inter$

#inter$state_id <- i

#inter

pdelta

}

out
}




run_heir <- function(model = "logistic",
                     fill = 0.001,
                     eps = 0.001,
                     iter_warmup = 1250,
                     iter_sampling = 1250){

  
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


if(model %in% c("logistic","exp")){
  y <- ifelse(df$y ==0, fill,df$y)
}else{
  y <- df$y
}


y_0 <- y[which(df$t == 1)]


stan_list <-  list(N = n*n_s,
                   y =y,
                   t = df$t,
                   t_cutoff = 14,
                   state = df$state_id,
                   n_s = n_s,
                   eps = eps,
                   y_0 = y_0)

if(model == "logistic"){
smod <- cmdstanr::cmdstan_model(here::here("growth-files/growth-model-states-gamma-simple.stan"))
}

if(model == "hurdle"){

smod <- cmdstanr::cmdstan_model(here::here("growth-files/growth-hurdle-many.stan"))

}

if(model == "exp"){
  smod <- cmdstanr::cmdstan_model(here::here("growth-files/growth-exp-many.stan"))
}

if(model =="exp-hurdle"){
  smod <- cmdstanr::cmdstan_model(here::here("growth-files/growth-exp-hurdle-many.stan"))
}

fit <- smod$sample(data = stan_list,
                   iter_warmup = iter_warmup,
                   iter_sampling = iter_sampling)

sumfit <-  fit$summary()


pdelta <-  sumfit %>% dplyr::filter(grepl("p_delta",variable))
pdelta$state_id <- c(1:51)

pdelta$model <- paste(model,"heir")

pdelta


}







h_sin <- run_many(n = n,
                  df = df,
                  t_cutoff = 14,
                  c_ints = c_ints,
                  eps = 0.001,
                  hurdle = T)

g_sin <- run_many(n = n,
                  df = df,
                  t_cutoff = 14,
                  c_ints = c_ints,
                  eps = 0.001,
                  hurdle = F)






