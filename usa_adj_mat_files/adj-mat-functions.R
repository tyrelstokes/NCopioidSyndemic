# Create state abreviation matrices ----------

state_adj_fun <- function(usa_adj,
                          state_abr = "NY",
                          fips = F){
  
  select_state_cols <- which(grepl(paste0(".",state_abr),colnames(usa_adj)) == TRUE)
  select_state_rows <- which(grepl(state_abr,usa_adj$X) == TRUE)
  
  state_mat <- usa_adj[select_state_rows, select_state_cols]
  rownames(state_mat) <- usa_adj$X[select_state_rows]
  
  if(fips ==T){
   county_ids <- usa_adj$county_id
   county_ids <- county_ids[select_state_rows]
   county_df <- data.frame(county_ids = county_ids,
                           county_name = usa_adj$X[select_state_rows])
   out <- list(state_mat = state_mat,
               county_df = county_df)
  }else{
  out <-  state_mat
  }
  
 
  out
}


######

map_adj_ids <- function(adj_list,
                        state_df){

  
  adj <- adj_list$state_mat
  fips <- adj_list$county_df  
  
  fips$remapped_id <- c(1:nrow(fips))
  
  state_df$remapped_id <- plyr::mapvalues(state_df$county_fips_code,
                                    from = fips$county_ids,
                                    to = fips$remapped_id)
  
  state_df
}



