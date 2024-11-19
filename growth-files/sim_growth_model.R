`%do%` <- foreach::`%do%`



growth_sim <- function(n,
                       n_s,
                       tau_k,
                       tau_delta,
                       tau_L,
                       tau_x,
                       alpha_k,
                       alpha_delta,
                       alpha_l,
                       alpha_x,
                       sigma,
                       t_0){
  
  
  k <- alpha_k + rnorm(n_s)*tau_k
  delta <- alpha_delta + rnorm(n_s)*tau_delta
  L <-  boot::inv.logit(alpha_l + rnorm(n_s)*tau_L)
  x_0 <- boot::inv.logit(alpha_x + rnorm(n_s)*tau_x)
  
  t <- c(1:n)/n
  
  mu <-matrix(nrow = n,unlist(lapply(c(1:n_s),function(i){
    x <- vector(length = n)
    for(j in 1:n){
      if(t[j] >= t_0){
   x[j] <-  L[i]*boot::inv.logit((k[i]+ delta[i])*(t[j]-t_0))
   

      }else{
        
   x[j] <-  L[i]*boot::inv.logit((k[i])*(t[j]-t_0))
        
      }
      
  }
    
   x} 
    )))
  
 obs_df <- foreach::foreach(i=1:n_s,.combine = rbind)%do%{
  
   y <- truncnorm::rtruncnorm(n = n, a = 0, b =1, mean = mu[,i],sd = sigma)
   df <- data.frame(y = y, state_int = i, t= t)
    df 
  }
    
 
 
 out <- list(k = k,
             delta = delta,
             L = L,
             x_0 = x_0,
             obs_df = obs_df) 
  
out  
}
  

n <- 17
n_s <- 10
tau_k <- 0.3
tau_delta <- 0.2
tau_x <- 0.1

alpha_k <- 1
alpha_delta <- 0
alpha_l <- 0
alpha_x <- 0
sigma <- 0.01
t_0 <- 0.7

growth_sim(n = n,
           n_s = n_s,
           tau_k = tau_k,
           tau_delta = tau_delta,
           tau_L = tau_L,
           tau_x = tau_x,
           alpha_k = alpha_k,
           alpha_delta = alpha_delta,
           alpha_l = alpha_l,
           alpha_x = alpha_x,
           sigma = sigma,
           t_0 = t_0)




growth_stan <- function(n = n,
                        n_s = n_s,
                        tau_k = tau_k,
                        tau_delta = tau_delta,
                        tau_L = tau_L,
                        tau_x = tau_x,
                        alpha_k = alpha_k,
                        alpha_delta = alpha_delta,
                        alpha_l = alpha_l,
                        alpha_x = alpha_x,
                        sigma = sigma,
                        t_0 = t_0){
  
  total_sim <- growth_sim(n = n,
                          n_s = n_s,
                          tau_k = tau_k,
                          tau_delta = tau_delta,
                          tau_L = tau_L,
                          tau_x = tau_x,
                          alpha_k = alpha_k,
                          alpha_delta = alpha_delta,
                          alpha_l = alpha_l,
                          alpha_x = alpha_x,
                          sigma = sigma,
                          t_0 = t_0)
  
  obs_df <- total_sim$obs_df
  
  
  stan_list <- list(N = n*n_s,
                    y = obs_df$y,
                    t = obs_df$t,
                    t_cutoff = t_0,
                    state = obs_df$state_int,
                    n_s = n_s)
  
  
  
  smod <- cmdstanr::cmdstan_model(here::here("growth-files/growth-model-states-gamma.stan"))
  

  fit <- smod$sample(data = stan_list,
                     iter_warmup = 2000,
                     iter_sampling = 2000)
  

  
}
