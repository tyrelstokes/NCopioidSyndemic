
plot_state <- function(state_df,
                       random_counties = F,
                       n_sample = 5,
                       which_counties = NULL,
                       ipsum = F){

pdf <- data.table::melt(data.table::setDT(state_df),id.vars = c("county_fips_code",
                                                          "month_id"),
                                             measure.vars = c("cases",
                                                              "jhu_cases",
                                                              "pred_cdc",
                                                              "avg_pred",
                                                              "avg"))


if(random_counties == F){
  w_count <- which_counties
}else{
  
  unique_counties <- unique(state_df$county_fips_code)
  
  n_count <- length(unique_counties)
  
  w_count <- sample(unique_counties,size = n_sample)
  
  
}

plot_out <- pdf %>% 
              dplyr::filter(county_fips_code %in% w_count) %>%
              ggplot2::ggplot(ggplot2::aes(x = month_id, y = value))+
              ggplot2::geom_line(ggplot2::aes(group = variable,color = variable),alpha = .3)+
              ggplot2::facet_wrap(~county_fips_code)
            
if(ipsum == T){
  plot_out <- plot_out + hrbrthemes::theme_ipsum()
}

plot_out

}