inla_predict <- function(inla_mod,
                         n,
                         np,
                         df){
  
  preds <- inla_mod$summary.fitted.values
  if(np == 3){
  preds$cases <- c(df$cases,df$cases,df$cases)
  preds$jhu_cases <- c(df$jhu_cases,df$jhu_cases,df$jhu_cases)
  preds$county <- c(df$county_fips_code,df$county_fips_code,df$county_fips_code)
  preds$month <- c(df$month_id,df$month_id,df$month_id)
  
  preds_cdc <- preds[1:nrow(df),]
  preds_jhu <- preds[(2*nrow(df)+1):(3*nrow(df)),]
  preds_zero <- preds[(nrow(df)+1):(2*nrow(df)),]
  
  
  df$pred_cdc <- preds_cdc$`0.5quant`
  df$pred_jhu <- preds_jhu$`0.5quant`
  df$pred_zero <- preds_zero$`0.5quant`
  df$jhu_comb <- (1-df$pred_zero)*df$pred_jhu
  df$avg_pred <- (df$jhu_comb+df$pred_cdc)/2
  df$avg <- (df$jhu_cases+df$cases)/2
  
  }else{
    
    preds$cases <- c(df$cases,df$cases)
    preds$jhu_cases <- c(df$jhu_cases,df$jhu_cases)
    preds$county <- c(df$county_fips_code,df$county_fips_code)
    preds$month <- c(df$month_id,df$month_id)
    
    preds_cdc <- preds[1:nrow(df),]
    preds_jhu <- preds[(nrow(df)+1):(2*nrow(df)),]

    
    
    df$pred_cdc <- preds_cdc$`0.5quant`
    df$pred_jhu <- preds_jhu$`0.5quant`

    df$jhu_comb <-  df$pred_jhu
    df$avg_pred <- (df$jhu_comb+df$pred_cdc)/2
    df$avg <- (df$jhu_cases+df$cases)/2
    
    
    
  }
  
  
  
  df$diff1 <- abs(df$avg - df$avg_pred)
  df$diff2 <- ifelse(is.na(df$avg),abs(df$pred_cdc - df$cases),df$diff1 ) 
  
  
  out <- list(df = df,
              preds = preds)
  
  out
  
  
  
}

#### joint/general prediction function -----------


inla_pred_list <- function(preds,
                             n,
                             np,
                             out_names,
                             y){
  out <- foreach::foreach(i = 1:np)%do%{
    
    beg <- 1 + (i-1)*n
    end <- beg -1 + n
    
    actual <- y[beg:end,i]
    
    inter <- preds[beg:end,]
    int_nms <- names(inter)
    names(inter) <- paste0(out_names[i],"_",int_nms)
    
    inter$actual <- actual
    inter$abs_diff <- abs(inter[,1] - inter$actual)
    
    inter
    
  }
  
  out
}


## outcome type function ----------

outcome_type_fun <- function(cases_outcomes,
                             hosp_outcomes,
                             deaths_outcomes){
  
  n_c <- length(cases_outcomes)
  n_h <- length(hosp_outcomes)
  n_d <- length(deaths_outcomes)
  
  nt <- n_c + n_h + n_d
  
  outcome_type <- foreach::foreach(i =1:nt,.combine = c)%do%{
    if(i <= n_c){
      x <- cases_outcomes[i]
      if(is.null(x) ==F){
        out <- "cases"
      }else{
        out <- NULL
      }
    }else{
      if(i <= n_c + n_h){
        x <- hosp_outcomes[i - n_c]
        if(is.null(x) ==F){
          out <- "hosp"
        }else{
          out <- NULL
        }
        
      }else{
        x <- deaths_outcomes[i - n_c - n_h]
        if(is.null(x) ==F){
          out <- "deaths"
        }else{
          out <- NULL
        }
      }
    }
    out
  }
  
  outcome_type
}



single_group_prediction <- function(type_list,
                                    pred_inds,
                                    group_likelihoods,
                                    name = "not_named"){
  
  nliks <- length(group_likelihoods)
  
  if(nliks == 1){
    out <- data.frame(total_pred = type_list[[pred_inds]])

  }else{
    
    which_binom <- which(group_likelihoods == "binomial")
    ind <- pred_inds[which_binom]
    binom_pred <- type_list[[which_binom]]
    
    
    pois_ind <- pred_inds[-which_binom]
    
    poisson_pred <- type_list[[pois_ind]]
    
    out <- data.frame(is_zero_pred = binom_pred, count_pred = poisson_pred)
    out$total_pred <- (1.0 - out$is_zero_pred)*out$count_pred
    
  }
  
  out$name <- "not_named"
  out
  
}


group_prediction <- function(type_list,
                             type_groups,
                             type_likelihoods,
                             type_names = NULL){
  
  
  unique_groups <- unique(type_groups)
  n_groups <- length(unique_groups)
  
  
 out <-  foreach::foreach(i = 1:n_groups)%do%{
   
   active_group <- unique_groups[i]
   group_inds <- which(type_groups == active_group)
  
   if(is.null(type_names)==F){
   nm <- type_names[group_inds[1]]
     
   }else{
     nm = "not_named"
   }
   
  single_group_prediction(type_list = type_list,
                          pred_inds = group_inds,
                          group_likelihoods = type_likelihoods[group_inds],
                          name = nm)
    
    
 }
  
  out
}


group_pred_df <- function(grouped_predictions,
                          type,
                          month_id_vec = NULL,
                          county_id_vec = NULL){
  
  n_g <- length(grouped_predictions)
  
out <-   foreach::foreach(i=1:n_g,.combine = cbind)%do%{
  inter <-   grouped_predictions[[i]]
  names(inter) <- paste0(names(inter),"_g_",i)
  inter
}

out$type <- type

if(is.null(month_id_vec)==F){
  out$month_id <- month_id_vec
}

if(is.null(county_id_vec)==F){
  out$county_id <- county_id_vec
}

out
  
}


average_type_pred <- function(grouped_predictions){
  
 n_groups <- length(grouped_predictions) 
  
total_preds <- foreach::foreach(i=1:n_groups,.combine = cbind)%do%{
  grouped_predictions[[i]]$total_pred
  
}

avg_pred <- apply(total_preds,1,mean)

 avg_pred 
}


pred_by_type <- function(pred_list,
                         outcome_type,
                         likelihoods,
                         groups,
                         which_pred = "mean",
                         type_names = NULL,
                         month_id_vec = NULL,
                         county_id_vec = NULL){
  
  
  unique_types <- unique(outcome_type)
  nt <- length(unique_types)
  
#  grouped_list <- vector("list",length)
  
 out <- foreach::foreach(i = 1:nt)%do%{
    
    tp <- unique_types[i]
    
  inds <-    which(outcome_type == tp)
    
   liks <- likelihoods[inds]
   
   unique_liks <-unique(likelihoods)
   
   type_groups <- groups[inds]
   unique_type_groups <- unique(type_groups)
   type_likelihoods <- likelihoods[inds]
   
  type_list <-  foreach::foreach(j = inds)%do%{
    
   preds <-  pred_list[[j]]
   
   x <- preds[,grepl(which_pred,names(preds))]
   x
   
  }
  
  
grouped_predictions <-   group_prediction(type_list = type_list,
                               type_groups = type_groups,
                               type_likelihoods = type_likelihoods,
                               type_names = type_names)
  
group_pred_df <- group_pred_df(grouped_predictions = grouped_predictions,
                               type = type_names,
                               month_id_vec = month_id_vec,
                               county_id_vec = county_id_vec) 
avg_pred <-  average_type_pred(grouped_predictions = grouped_predictions)

inter <- list(avg_pred = avg_pred,
              group_pred_df = group_pred_df,
              grouped_predictions = grouped_predictions)

inter
    
  }
 
 out 
}


avg_pred_combine <- function(outcome_type,
                             inter_list,
                             month_id_vec = NULL,
                             county_id_vec = NULL){
  unique_cases <- unique(outcome_type) 
  
  n_types <- length(unique_cases)
  
 out<- foreach::foreach(i=1:n_types,.combine = 'cbind')%do%{
   x <-  inter_list[[i]]$avg_pred
   nm <- paste0(unique_cases[i],"_avg_pred")
  inter <-  data.frame(x = x)
  names(inter) <- nm
  
  if(is.null(month_id_vec)==F){
    inter$month_id <- month_id_vec
  }
  
  if(is.null(county_id_vec)==F){
    inter$county_id <- county_id_vec
  }
  
  inter
 }
 out
}

get_grouped_list <- function(outcome_type,
                             inter_list){
  
  unique_cases <- unique(outcome_type) 
  
  n_types <- length(unique_cases)
  
 out <- foreach::foreach(i = 1:n_types)%do%{
   inter_list[[i]]$group_pred_df
 }
 
 out
}

inla_predict_gen <- function(inla_mod,
                             n,
                             np,
                             df,
                             cases_outcomes = NULL,
                             hosp_outcomes = NULL,
                             deaths_outcomes = NULL,
                             cases_likelihoods = NULL,
                             hosp_likelihoods = NULL,
                             deaths_likelihoods = NULL,
                             cases_group_mod = NULL,
                             hosp_group_mod = NULL,
                             deaths_group_mod = NULL,
                             y = y,
                             which_pred = "mean",
                             month_id_vec = NULL,
                             county_id_vec = NULL){
  
  
  preds <- inla_mod$summary.fitted.values
  
  
  outcome_names <- c(cases_outcomes,hosp_outcomes, deaths_outcomes)
  
  n_c <- length(cases_outcomes)
  n_h <- length(hosp_outcomes)
  n_d <- length(deaths_outcomes)
  
  
  likelihoods <- c(cases_likelihoods,hosp_likelihoods,deaths_likelihoods)
  
  
 outcome_type <-  outcome_type_fun(cases_outcomes = cases_outcomes,
                                   hosp_outcomes = hosp_outcomes,
                                   deaths_outcomes = deaths_outcomes)

 
 ind_null <- which(is.null(outcome_names))
 
 if(length(ind_null)>0){
   outcome_names <- outcome_names[!ind_null]
   outcome_type <- outcome_type[!ind_null]
   likelihoods <- likelihoods[!ind_null]
 }
  
  pred_list <- inla_pred_list(preds = preds,
                              n = n,
                              np = np,
                              out_names = outcome_names,
                              y =y)
  
  
  groups <- c(cases_group_mod,hosp_group_mod,deaths_group_mod)
  
 inter_list <-  pred_by_type(pred_list = pred_list,
                             outcome_type = outcome_type,
                             likelihoods = likelihoods,
                             groups = groups,
                             which_pred = which_pred,
                             type_names = NULL,
                             month_id_vec = month_id_vec,
                             county_id_vec = month_id_vec) 
  
  

 avg_pred_df <- avg_pred_combine(outcome_type = outcome_type,
                                 inter_list = inter_list,
                                 month_id_vec = month_id_vec,
                                 county_id_vec = county_id_vec) 
 
  
 
grouped_list <- get_grouped_list(outcome_type = outcome_type,
                                 inter_list = inter_list)


out <- list(pred_list = pred_list,
            avg_pred_df = avg_pred_df,
            grouped_list = grouped_list)
  
  
out
}


