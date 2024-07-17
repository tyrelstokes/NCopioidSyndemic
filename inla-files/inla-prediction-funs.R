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
