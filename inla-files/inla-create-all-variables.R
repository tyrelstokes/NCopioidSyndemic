rep_vector <- function(np,
                       x){
  
  out <- vector("list", length = np)
  
  for(i in 1:np){
    out[[i]] <- x
  }
  
out  
  
}


time_vector <- function(np,
                        x,
                        base_name,
                        n){
  
  
 xlist <- rep_vector(np = np,
                     x = x)
  
 out <-  inla_vectors(n = n,
                            np = np,
                            x_list = xlist,
                            base_name = base_name)
  
 out 
}







# create additional intercepts ------------------


inla_time_intercepts <- function(df,
                                 np,
                                 march = T,
                                 aft_2021 = T,
                                 late = T,
                                 n){
  
  
 effs <- c(march,aft_2021,late)
  
  n_ef <- sum(effs)
  
  out <- vector("list",length = n_ef)
  j<- 1
  
  if(march == T){
  
  df$after_march_2020 <- ifelse(df$case_month > "2020-02-01",1,0)
  
 march_list <-  time_vector(np = np,
                          x = df$after_march_2020,
                          base_name = "march",
                          n = n)
  
 out[[j]] <- march_list
 
 names(out[[j]]) <- "march_list"
 
 j <- j+1
 
  }
  
if(aft_2021 == T){
  
  df$after_2021 <- ifelse(df$case_month > "2020-12-01",1,0)
  
  vac_list <-  time_vector(np = np,
                             x = df$after_2021,
                             base_name = "vac",
                             n = n)
  
  out[[j]] <- vac_list
  
  names(out[[j]]) <- "vac_list"
  
  j <- j+1
  
}
  
  
if(late == T){
  
  df$after_2023 <- ifelse(df$case_month > "2022-12-01",1,0)
  
  late_list <- time_vector(np = np,
                           x = df$after_2023,
                           base_name = "late",
                           n = n)
  
  out[[j]] <- late_list
  
  names(out[[j]]) <- "late_list"
  
  j <- j+1
  
}
  

  
out  
  
  
  
  
}















# Functions to put all variables in a list ----------


inla_data_list <- function(df,
                           n,
                           np,
                           outcome_list,
                           offset_list,
                           id_list,
                           month_list,
                           extra_time_int = TRUE,
                           outcome_inters = NULL,
                           month_list_lag = NULL){
  
  
  # Create INLA outcome matrix for multidimensional modelling
  y <- inla_matrix(n = n,
                   np = np,
                   y_list = outcome_list)
  
  # Create proper intercepts 
  
  int_xlist <- rep_vector(np = np,
                           x = rep(1,n))
  int_list <- inla_vectors(n =n,
                           np = np,
                           x_list = int_xlist,
                           base_name = "int")
  
  
  # Create the month indices
  
  month_list <- inla_vectors(n = n,
                             np = np,
                             x_list = month_list,
                             base_name = "month")
  
  
  # Offset List
  
  offset_list <- inla_vectors(n = n,
                              np = np,
                              x_list = offset_list,
                              base_name = "E")
  
  
  # Id list
  
  id_list <- inla_vectors(n = n,
                          np = np,
                          x_list = id_list,
                          base_name = "id")
  
  # create time intercepts if desired
  
  if(is.null(month_list_lag)==FALSE){
    month_list_lag <-  inla_vectors(n = n,
                                    np = np,
                                    x_list = month_list_lag,
                                    base_name = "month_lag")
  }
  
  
  
  if(extra_time_int == T){
  time_lists <-   inla_time_intercepts(df = df,
                                     np = np,
                                     march = T,
                                     aft_2021 = T,
                                     late = T,
                                     n = n) 
  
  march_list <- time_lists[[1]]
  vac_list <- time_lists[[2]]
  late_list <- time_lists[[3]]
  
  out <- list(y = y,
              int_list = int_list,
              month_list = month_list,
              offset_list = offset_list,
              id_list = id_list,
              time_lists = time_lists,
              month_list_lag = month_list_lag)
    
  }else{
    
    out <- list(y = y,
                int_list = int_list,
                month_list = month_list,
                offset_list = offset_list,
                id_list = id_list,
                month_list_lag = month_list_lag)  
    
    
    
  }
  
out  
  
}


# hospital outcomes create ---------

hosp_outcome_create <- function(df){
  df$hosp_zero <- ifelse(df$hosp_yes ==0,1,df$hosp_yes)
  df$new_hosp_zero <- ifelse(df$new_hosp_total ==0, 1, df$new_hosp_total)
  
  df
}

# death outcomes create -----------

death_outcome_create <- function(df){
  
  df$death_zero <- ifelse(df$death_yes ==0,1,0)
  df$covid_death_zero <- ifelse(df$CovidDeathCount ==0,1,0)
  df$jhu_deaths_zero <- ifelse(df$jhu_deaths ==0,1,0)
  
  df
}

