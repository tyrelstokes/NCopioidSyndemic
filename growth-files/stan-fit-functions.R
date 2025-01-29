# housekeeping ------------------------------------
`%>%` <- dplyr::`%>%`
`%do%` <- foreach::`%do%`



# Compute mu (needed for log lik and WAIC) --------------------------------

compute_mu_for_draw <- function(draw, L, k, delta, x_0, t) {
  # Extract parameters for this draw
  L_draw <- L[draw]
  k_draw <- k[draw]
  delta_draw <- delta[draw]
  x_0_draw <- x_0[draw]
  
  # Compute t_scale (normalize time vector t)
  t_min <- min(t)
  t_max <- max(t)
  t_scale <- (t - t_min) / (t_max - t_min)
  
  # Set t_0 (cutoff for time)
  t_cutoff <- 14  # Adjust this value as per your model
  t_0 <- (t_cutoff - t_min) / (t_max - t_min)
  
  # Compute k_eff based on t_scale
  k_eff <- ifelse(t_scale >= t_0, k_draw + delta_draw, k_draw)
  
  # Compute mu
  mu <- L_draw * plogis(k_eff * (t_scale - x_0_draw))
  return(mu)
}

# single function --------------

run_single <- function(n,
                       j,
                       df,
                       t_cutoff,
                       eps = 0.001,
                       fill = 0.0001,
                       model = "logistic",
                       draws = F,
                       waic = T){
  
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
             t_cutoff = t_cutoff, # changed from 14
             eps = eps)
  
  if(model == "hurdle"){
    smods <- cmdstanr::cmdstan_model(here::here("growth-hurdle-single.stan"))
  }
  
  if(model == "logistic"){
    smods <- cmdstanr::cmdstan_model(here::here("growth-gamma-single.stan")) 
  }
  
  if(model == "exp"){
    smods <- cmdstanr::cmdstan_model(here::here("growth-exp-single.stan")) 
  }
  
  if(model == "exp-hurdle"){
    smods <- cmdstanr::cmdstan_model(here::here("growth-exp-single-hurdle.stan")) 
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
  
  ## calculate WAIC
  posterior_samples <- fit$draws(format = "df")
  t <- df_s$t
  L <- posterior_samples$L  # Example parameter, adapt as needed
  k <- posterior_samples$k
  delta <- posterior_samples$delta
  x_0 <- posterior_samples$x_0
  
  # Compute log-likelihood matrix for all posterior draws
  log_lik_matrix <- sapply(1:nrow(posterior_samples), function(draw) {
    mu <- compute_mu_for_draw(draw, L, k, delta, x_0, t)  # Recreate mu per draw
    log_lik <- dgamma(y, shape = mu, scale = 1 / mu, log = TRUE)  # Example likelihood
    return(log_lik)
  })
  
  # Transpose the matrix if needed (rows = posterior draws, columns = observations)
  log_lik_matrix <- t(log_lik_matrix)
  
  # Convert Inf to -Inf
  log_lik_df <- as.data.frame(log_lik_matrix)
  
  log_lik_cols <- colnames(log_lik_df)
  
  log_lik_df <- log_lik_df %>%
    dplyr::mutate(across(all_of(log_lik_cols), ~ replace(., . == Inf, -Inf)))
  
  # Compute WAIC
  waic_result <- loo::waic(as.matrix(log_lik_df))
  
  waic_df <- data.frame(elpd_waic = waic_result[["elpd_waic"]],
                        se_elpd_waic = waic_result[["se_elpd_waic"]],
                        p_waic = waic_result[["p_waic"]],
                        se_p_waic = waic_result[["se_p_waic"]],
                        waic = waic_result[["waic"]],
                        se_waic = waic_result[["se_waic"]])
  
  if(draws == F){
    
    if(waic == F) {
      out <- list(r2 = r2,
                  param_df = param_df,
                  df_pred = df_pred,
                  pdelta = pdelta)
    } else {
      
      # include WAIC
      out <- list(r2 = r2,
                  param_df = param_df,
                  df_pred = df_pred,
                  pdelta = pdelta,
                  waic = waic_df)
    }
    
    
  }else{
    
    draw_df <- fit$draws(variables = "delta",format = "df")
    draw_df$state_id <- j
    draw_df$model <- model
    
    if(waic == F) {
    
      out <- list(r2 = r2,
                  param_df = param_df,
                  df_pred = df_pred,
                  pdelta = pdelta,
                  draw_df = draw_df)  
    } else {
      
      # include WAIC
      out <- list(r2 = r2,
                  param_df = param_df,
                  df_pred = df_pred,
                  pdelta = pdelta,
                  draw_df = draw_df,
                  waic = waic_df)  
    }
    
  }
  
  out
}

# single many times function -------------------------
run_many <- function(n,
                     df,
                     t_cutoff,
                     c_ints,
                     eps = 0.001,
                     model = "logistic",
                     fill = 0.0001,
                     draws = F,
                     snames = colnames(nflis)[c_ints],
                     prob_fill = F,
                     waic = T){
  
  if(draws == T){
    draw_list <- vector("list", length = length(c_ints))
  }
  
  
  pred_list <- vector("list",length = length(c_ints))
  
  out <- foreach::foreach(i = c_ints,.combine = rbind)%do%{
    inter <- run_single(n = n,
                        j = i,
                        d = df,
                        t_cutoff = t_cutoff,
                        eps = eps,
                        fill = fill,
                        model = model,
                        draws = draws,
                        waic = waic)
    
    pdelta <- inter$pdelta
    pdelta$state_id <- i
    pdelta <- cbind(pdelta, inter$waic)
    
    #inter
    if(draws == T){
      ddf <- inter$draw_df
      ddf$state <- snames[i]
      draw_list[[i]] <- ddf
    }
    
    
    preds <- inter$df_pred
    preds$state <- i
    
    pred_list[[i]] <- preds
    
    return(pdelta)
    
    
    
  }
  
  if(draws == T){
    draw_cmb <- do.call(rbind,draw_list)
    
    mn_df <- draw_cmb %>% dplyr::group_by(state) %>% dplyr::summarise(avg = mean(delta))
    mn_df$prob <- plyr::mapvalues(x = mn_df$state,
                                  from = mn_df$state,
                                  to = out$mean)
    
    if(prob_fill == T){
      mn_df <- mn_df %>% dplyr::arrange(prob)
      
    }else{
      
      mn_df <- mn_df %>% dplyr::arrange(avg) 
    }
    
    
    draw_cmb$state <- factor(x = draw_cmb$state,levels= mn_df$state)
    draw_cmb$mn <- plyr::mapvalues(x = draw_cmb$state,from = mn_df$state,to = as.numeric(mn_df$avg))
    draw_cmb$mn <- as.numeric(as.character(draw_cmb$mn))
    
    draw_cmb$prob <- plyr::mapvalues(x = draw_cmb$state,from = mn_df$state,to = as.numeric(mn_df$prob))
    
    
    if(prob_fill == F){
      dplot <- draw_cmb %>%
        ggplot2::ggplot(ggplot2::aes(x = delta,
                                     y = state)) + ggridges::geom_density_ridges2(ggplot2::aes(fill = mn))+
        ggplot2::scale_fill_continuous(low = "firebrick", high = "dodgerblue")
    }else{
      dplot <- draw_cmb %>%
        ggplot2::ggplot(ggplot2::aes(x = delta,
                                     y = state)) + ggridges::geom_density_ridges2(ggplot2::aes(fill = prob))+
        ggplot2::scale_fill_continuous(low = "firebrick", high = "dodgerblue")  
      
      
    }
    
    out <- list(pdelta = out,
                draw_cmb = draw_cmb,
                dplot = dplot,
                pred_list = pred_list)
  }else{
    
    out <- list(pdelta = out,
                pred_list = pred_list)
    
  }
  
  out
}


# run hierarchical models --------------------------------

run_heir <- function(model = "logistic",
                     fill = 0.001,
                     eps = 0.001,
                     iter_warmup = 1250,
                     iter_sampling = 1250,
                     predictions = 1){
  
  
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
                     y_0 = y_0,
                     pred_ind = predictions,
                     snames = colnames(nflis))
  
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
  
  draw_df <- fit$draws(variables = "delta",format = "df")
  draw_cmb <- do.call(rbind,lapply(c(1:51),function(i){
    delta <-  unlist(as.vector(draw_df[,i]))
    state <- snames[i]
    ot <- data.frame(delta = delta,state = state)
    ot
  }))
  mn_df <- draw_cmb %>% dplyr::group_by(state) %>% dplyr::summarise(avg = mean(delta)) %>% dplyr::arrange(avg)
  
  draw_cmb$state <- factor(x = draw_cmb$state,levels= mn_df$state)
  draw_cmb$mn <- plyr::mapvalues(x = draw_cmb$state,from = mn_df$state,to = as.numeric(mn_df$avg))
  draw_cmb$mn <- as.numeric(as.character(draw_cmb$mn))
  
  dplot <- draw_cmb %>%
    ggplot2::ggplot(ggplot2::aes(x = delta,
                                 y = state)) + ggridges::geom_density_ridges2(ggplot2::aes(fill = mn))+
    ggplot2::scale_fill_continuous(low = "firebrick", high = "dodgerblue") + 
    ggplot2::geom_vline(xintercept = 0)
  
  #pdelta
  
  if(predictions == 1){
    pred_df <- sumfit %>% dplyr::filter(grepl("pred",variable))
    pred_df$t <- df$t
    pred_df$state_id <- df$state_id
    
    
    draws_df <- fit$draws(format = "df")
    
    out <- list(pdelta = pdelta,
                pred_df = pred_df,
                dplot = dplot)
    
    
  }else{
    out = pdelta
  }
  
  out
}


# run hierarchical models (alt with any time point) -----------------------


run_heir_alt <- function(model = "logistic",
                         fill = 0.001,
                         eps = 0.001,
                         iter_warmup = 1250,
                         iter_sampling = 1250,
                         predictions = 1,
                         t_cutoff,
                         sname_list = c(),
                         waic = T){
  
  
  nflis <- read.csv(here::here("nflis_data_RBeast.csv"))
  nflis_test <- nflis %>% dplyr::select(all_of(sname_list))
  
  
  df <-   do.call(rbind,lapply(c(1:ncol(nflis_test)),
                               function(i){
                                 y <- nflis_test[,i]
                                 t <- c(1:nrow(nflis_test))
                                 state_id <- rep(i,nrow(nflis_test))
                                 
                                 out <- data.frame(y = y,t = t, state_id = state_id)
                                 out
                               }))
  
  
  
  n <- nrow(nflis_test)
  n_s <- ncol(nflis_test)
  
  
  if(model %in% c("logistic","exp")){
    y <- ifelse(df$y ==0, fill,df$y)
  }else{
    y <- df$y
  }
  
  
  y_0 <- y[which(df$t == 1)]
  
  
  stan_list <-  list(N = n*n_s,
                     y =y,
                     t = df$t,
                     t_cutoff = t_cutoff,
                     state = df$state_id,
                     n_s = n_s,
                     eps = eps,
                     y_0 = y_0,
                     pred_ind = predictions)
  
  snames = colnames(nflis_test)[c_ints]
  
  if(model == "logistic"){
    #smod <- cmdstanr::cmdstan_model(here::here("growth-model-states-gamma-simple.stan"))
    smod <- cmdstanr::cmdstan_model(here::here("growth-model-states-gamma-simple_alt.stan"))
  }
  
  if(model == "hurdle"){
    
    #smod <- cmdstanr::cmdstan_model(here::here("growth-hurdle-many.stan"))
    smod <- cmdstanr::cmdstan_model(here::here("growth-hurdle-many_alt.stan"))
    
    
  }
  
  if(model == "exp"){
    #smod <- cmdstanr::cmdstan_model(here::here("growth-exp-many.stan"))
    smod <- cmdstanr::cmdstan_model(here::here("growth-exp-many_alt.stan"))
  }
  
  if(model =="exp-hurdle"){
    #smod <- cmdstanr::cmdstan_model(here::here("growth-exp-hurdle-many.stan"))
    smod <- cmdstanr::cmdstan_model(here::here("growth-exp-hurdle-many_alt.stan"))
  }
  
  fit <- smod$sample(data = stan_list,
                     iter_warmup = iter_warmup,
                     iter_sampling = iter_sampling)
  
  sumfit <-  fit$summary()
  
  
  pdelta <-  sumfit %>% dplyr::filter(grepl("p_delta",variable))
  pdelta <- pdelta %>% dplyr::filter(!is.na(mean))
  pdelta$model <- paste(model,"heir")
  pdelta$state_id <- c(1:n_s)
  
  log_lik <- fit$draws("log_lik")
  log_lik_matrix <- posterior::as_draws_matrix(log_lik)
  loo <- loo::loo(log_lik_matrix)
  
  pdelta$elpd_loo <- loo$elpd_loo
  pdelta$p_loo <- loo$p_loo
  pdelta$looic <- loo$looic
  
  pdelta$se_elpd_loo <- loo$se_elpd_loo
  pdelta$se_p_loo <- loo$se_p_loo
  pdelta$se_looic <- loo$se_looic
  
  pdelta$t_cutoff <- t_cutoff
  
  draw_df <- fit$draws(variables = "delta",format = "df")
  draw_cmb <- do.call(rbind,lapply(c(1:n_s),function(i){
    delta <-  unlist(as.vector(draw_df[,i]))
    state <- snames[i]
    ot <- data.frame(delta = delta,state = state)
    ot
  }))
  mn_df <- draw_cmb %>% dplyr::group_by(state) %>% dplyr::summarise(avg = mean(delta)) %>% dplyr::arrange(avg)
  
  draw_cmb$state <- factor(x = draw_cmb$state,levels= mn_df$state)
  draw_cmb$mn <- plyr::mapvalues(x = draw_cmb$state,from = mn_df$state,to = as.numeric(mn_df$avg))
  draw_cmb$mn <- as.numeric(as.character(draw_cmb$mn))
  
  dplot <- draw_cmb %>%
    ggplot2::ggplot(ggplot2::aes(x = delta,
                                 y = state)) + ggridges::geom_density_ridges2(ggplot2::aes(fill = mn))+
    ggplot2::scale_fill_continuous(low = "firebrick", high = "dodgerblue") + 
    ggplot2::geom_vline(xintercept = 0)
  
  #pdelta
  
  # if(predictions == 1){
  #   pred_df <- sumfit %>% dplyr::filter(grepl("pred",variable))
  #   pred_df$t <- df$t
  #   pred_df$state_id <- df$state_id
  #   
  #   
  #   draws_df <- fit$draws(format = "df")
  #   
  #   out <- list(pdelta = pdelta,
  #               pred_df = pred_df,
  #               dplot = dplot)
  #   
  #   
  # }else{
  #   out = pdelta
  # }
  
  if (waic == T) {
    
    out <- list(pdelta,
                loo)
    
  } else {
    out = pdelta
  }
  
  out
}



# OLD CODE ----------------------------------------------------------------
# 
# # rbind WAIC
# ## calculate WAIC
# # posterior_samples <- fit$draws(format = "df")
# # t <- df_s$t
# # L <- posterior_samples$L  # Example parameter, adapt as needed
# # k <- posterior_samples$k
# # delta <- posterior_samples$delta
# # x_0 <- posterior_samples$x_0
# 
# posterior_samples <- fit$draws(format = "df")
# t <- df$t
# L <- posterior_samples %>% dplyr::select(starts_with("L")) %>% dplyr::select(!starts_with("L_0")) %>% dplyr::select(!starts_with("lp"))
# k <- posterior_samples %>% dplyr::select(starts_with("k["))
# delta <- posterior_samples %>% dplyr::select(starts_with("delta["))
# x_0 <- posterior_samples %>% dplyr::select(starts_with("x_0"))
# 
# 
# model_mu <- posterior_samples %>% dplyr::select(starts_with("mu"))
# fit_mu <- sumfit %>% dplyr::filter(grepl("mu",variable))
# 
# alpha <- posterior_samples %>% dplyr::select(starts_with("alpha["))
# b <- posterior_samples %>% dplyr::select(starts_with("b["))
# 
# invphi <- posterior_samples %>% dplyr::select(invphi)
# 
# waic_df_list <- c()
# 
# for (iState in 1:n_s) {
#   
#   if (dim(L)[2] > 1) {
#     this_L <- as.matrix(L[,iState])
#   } else {
#     this_L <- as.matrix(L)
#   }
#   
#   this_k <- as.matrix(k[,iState])
#   this_delta <- as.matrix(delta[,iState])
#   
#   if (dim(x_0)[2] > 1) {
#     this_x_0 <- as.matrix(x_0[,iState])
#   } else {
#     this_x_0 <- as.matrix(x_0)
#   }
#   
#   # this_alpha <- as.matrix(alpha[,iState])
#   # this_b <- as.matrix(b[,iState])
#   # 
#   # calc_mu <- exp(-this_k * t)
#   # beta <- invphi / calc_mu
#   # 
#   # log_lik <- numeric(length(y))
#   # for (iT in seq_along(y)) {
#   #   # Hurdle probability
#   #   hurdle_prob <- plogis(this_alpha + this_b * t[iT])
#   #   
#   #   if (y[iT] > 0) {
#   #     # Log-likelihood for crossing the hurdle and Gamma density
#   #     log_lik[iT] <- log(hurdle_prob) + dgamma(y[iT], shape = invphi, rate = beta, log = TRUE)
#   #   } else {
#   #     # Log-likelihood for not crossing the hurdle
#   #     log_lik[iT] <- log(1 - hurdle_prob)
#   #   }
#   # }
#   # 
#   # get_mu <- sapply(1:nrow(posterior_samples), function(draw) {
#   #   mu <- compute_mu_for_draw(draw, this_L, this_k, this_delta, this_x_0, df_s$t)  # Recreate mu per draw
#   #   return(mu)
#   # })
#   
#   # Compute log-likelihood matrix for all posterior draws
#   log_lik_matrix <- sapply(1:nrow(posterior_samples), function(draw) {
#     mu <- compute_mu_for_draw(draw, this_L, this_k, this_delta, this_x_0, df$t)  # Recreate mu per draw
#     log_lik <- dgamma(y, shape = mu, scale = 1 / mu, log = TRUE)  # Example likelihood
#     return(log_lik)
#     #return(posterior_samples[draw,])
#   })
#   
#   log_lik_matrix <- t(log_lik_matrix)
#   
#   was_debugged <- 0
#   
#   # Replace Inf values with a large negative number only if there is an error
#   waic_result <- tryCatch(
#     {
#       # Attempt to compute WAIC without modifications
#       loo::waic(log_lik_matrix)
#       
#     },
#     error = function(e) {
#       message("Error in loo::waic: ", e$message)
#       message("Replacing Inf/-Inf values with large finite numbers...")
#       
#       was_debugged <<- 1
#       
#       log_lik_df <- as.data.frame(log_lik_matrix)
#       
#       # Replace Inf and -Inf with large numbers for debugging
#       log_lik_df_clean <- log_lik_df %>%
#         dplyr::mutate(across(everything(), ~ ifelse(is.infinite(.), -1e10, .)))
#       
#       # Confirm there are no remaining issues
#       if (any(!is.finite(as.matrix(log_lik_df_clean)))) {
#         stop("Still contains non-finite values after cleaning!")
#       }
#       
#       # Reattempt WAIC with modified log-likelihood matrix
#       loo::waic(as.matrix(log_lik_df_clean))
#       
#     }
#   )
#   
#   waic_df <- data.frame(elpd_waic = waic_result[["elpd_waic"]],
#                         se_elpd_waic = waic_result[["se_elpd_waic"]],
#                         p_waic = waic_result[["p_waic"]],
#                         se_p_waic = waic_result[["se_p_waic"]],
#                         waic = waic_result[["waic"]],
#                         se_waic = waic_result[["se_waic"]],
#                         was_debugged = was_debugged)
#   
#   waic_df_list <- rbind(waic_df_list, waic_df)
#   
# }
# 
# pdelta <- cbind(pdelta, waic_df_list)
# 
