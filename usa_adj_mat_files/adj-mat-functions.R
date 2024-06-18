# Create state abreviation matrices ----------

state_adj_fun <- function(usa_adj,
                          state_abr = "NY"){
  
  select_state_cols <- which(grepl(paste0(".",state_abr),colnames(usa_adj)) == TRUE)
  select_state_rows <- which(grepl(state_abr,usa_adj$X) == TRUE)
  
  state_mat <- usa_adj[select_state_rows, select_state_cols]
  rownames(state_mat) <- usa_adj$X[select_state_rows]
  
  state_mat
  
}
