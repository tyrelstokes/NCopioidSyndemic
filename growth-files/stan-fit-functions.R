# housekeeping ------------------------------------
`%>%` <- dplyr::`%>%`
`%do%` <- foreach::`%do%`



# single function --------------

run_single <- function(n,
                       j,
                       df,
                       t_cutoff,
                       eps = 0.001,
                       fill = 0.0001,
                       model = "logistic",
                       draws = F){
  
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
  
  if(draws == F){
  
  out <- list(r2 = r2,
              param_df = param_df,
              df_pred = df_pred,
              pdelta = pdelta)
  
  }else{
    draw_df <- fit$draws(variables = "delta",format = "df")
    draw_df$state_id <- j
    draw_df$model <- model
    
    out <- list(r2 = r2,
                param_df = param_df,
                df_pred = df_pred,
                pdelta = pdelta,
                draw_df = draw_df)  
    
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
                     prob_fill = F){
  
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
                        draws = draws)
    
    pdelta <- inter$pdelta
    pdelta$state_id <- i
    #inter$
    
    #inter$state_id <- i
    
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

