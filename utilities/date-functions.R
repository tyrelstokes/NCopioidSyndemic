# map the dates to a numeric 1:n_months -----------

map_month_ids <- function(df,
                          beg = NULL,
                          end = NULL){
  
  df$month <- as.Date(df$case_month) %>% lubridate::ymd()
  
  month_data <- as.Date(df$case_month)
  
  months <- unique(month_data) %>% lubridate::ymd()
  
  if(is.null(beg)==TRUE){
  beg <- min(months)
  end <- max(months)
  }
  
  full_seq <- beg %m+% months(seq(0, round(lubridate::interval(beg, end) / months(1)), 1))
  
  month_id <- c(1:length(full_seq))
  m_df <- data.frame(month = full_seq,
                     month_id = month_id)
  
  df <- dplyr::left_join(df,m_df, by = "month") 
  
 df 
  
}

