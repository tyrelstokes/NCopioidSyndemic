er <- ny %>% dplyr::filter(county_fips_code =="36029")




# Create INLA outcome matrix for multidimensional modelling

n <- nrow(er)
y <- inla_matrix(n = n,
                 np = np,
                 y_list = list(er$cases,
                               er$jhu_zero,
                               er$jhu_cases))

# Create proper intercepts 
int_list <- inla_vectors(n =n,
                         np = np,
                         x_list = list(rep(1,n),
                                       rep(1,n),
                                       rep(1,n)),
                         base_name = "int")

int1 <- int_list$int_1
int2 <- int_list$int_2
int3 <- int_list$int_3

# Create the E vectors for rate modelling
E_list <- inla_vectors(n = n,
                       np = np,
                       x_list = list(er$population_scale,
                                     er$population_scale,
                                     er$population_scale),
                       base_name = "E")


E1 <- E_list$E_1
E2 <- E_list$E_2
E3 <- E_list$E_3



# Create Ids for spatial effects

id_list <- inla_vectors(n = n,
                        np = np,
                        x_list = list(er$remapped_id,
                                      er$remapped_id,
                                      er$remapped_id),
                        base_name = "id")


id1 <- id_list$id_1
id2 <- id_list$id_2
id3 <- id_list$id_3

id_slm1 <- id1
id_slm2 <- id2
id_slm3 <- id3
# Create the month indices

month_list <- inla_vectors(n = n,
                           np = np,
                           x_list = list(er$month_id,
                                         er$month_id,
                                         er$month_id),
                           base_name = "month")

m1 <- month_list$month_1
m2 <- month_list$month_2
m3 <- month_list$month_3

# lag cases ----------


# intercepts for time events

er$after_march_2020 <- ifelse(er$case_month > "2020-02-01",1,0)

m20_list <-  inla_vectors(n = n,
                          np = np,
                          x_list = list(er$after_march_2020,
                                        er$after_march_2020,
                                        er$after_march_2020),
                          base_name = "march")

mar1 <- m20_list$march_1
mar2 <- m20_list$march_2
mar3 <- m20_list$march_3


er$after_2021 <- ifelse(er$case_month > "2020-12-01",1,0)

vac_list <-  inla_vectors(n = n,
                          np = np,
                          x_list = list(er$after_2021,
                                        er$after_2021,
                                        er$after_2021),
                          base_name = "vac")

vac1 <- vac_list$vac_1
vac2 <- vac_list$vac_2
vac3 <- vac_list$vac_3



















link_set <- c(rep(1,nrow(er)),
              rep(2,nrow(er)),
              rep(3,nrow(er)))




test_err <- INLA::inla(y ~ -1+ int1 + int2 +int3+
                           mar1 + mar2 + mar3 +
                           vac1 + vac2 + vac3+
                           f(m1,
                             model= "ar",
                             order = 2,
                             adjust.for.con.comp = TRUE)+
                           f(m2,
                             copy = "m1",
                             #group = m2,
                             hyper = list(beta = beta.prior))+
                           f(m3,
                             copy = "m1",
                            # group = m3,
                             hyper = list(beta = beta.prior)),
                         family=c("poisson","binomial","zeroinflatedpoisson0"),
                         data=list(y =y,
                                   E1 = E1,
                                   E2 = E2,
                                   E3 = E3,
                                   int1 = int1,
                                   int2 = int2,
                                   int3 = int3,
                                   id1 = id1,
                                   id2 = id2,
                                   id3 = id3,
                                   id_slm1 = id_slm1,
                                   id_slm2 = id_slm2,
                                   id_slm3 = id_slm3,
                                   l1 = l1,
                                   l2 = l2,
                                   l3 = l3,
                                   mar1 = mar1,
                                   mar2 = mar2,
                                   mar3 = mar3,
                                   vac1 = vac1,
                                   vac2 = vac2,
                                   vac3 = vac3),
                         control.compute=list(waic=TRUE),
                         control.predictor = list(compute = TRUE,
                                                  link = link_set),
                         control.inla=list(cmin=0),
                         control.fixed = list(prec = 1),
                         safe = TRUE,
                         verbose = TRUE)


summary(test_err)
######################

preds <- test_err$summary.fitted.values
preds$cases <- c(er$cases,er$cases,er$cases)
preds$jhu_cases <- c(er$jhu_cases,er$jhu_cases,er$jhu_cases)
preds$county <- c(er$county_fips_code,er$county_fips_code,er$county_fips_code)
preds$month <- c(er$month_id,er$month_id,er$month_id)

preds_cdc <- preds[1:nrow(er),]
preds_jhu <- preds[(2*nrow(er)+1):(3*nrow(er)),]
preds_zero <- preds[(nrow(er)+1):(2*nrow(er)),]

w_count <-  head(er$county_fips_code)

er$pred_cdc <- preds_cdc$`0.5quant`
er$pred_jhu <- preds_jhu$`0.5quant`
er$pred_zero <- preds_zero$`0.5quant`
er$jhu_comb <- (1-er$pred_zero)*er$pred_jhu
er$avg_pred <- (er$jhu_comb+er$pred_cdc)/2
er$avg <- (er$jhu_cases+er$cases)/2

er$diff1 <- abs(er$avg - er$avg_pred)
er$diff2 <- ifelse(is.na(er$avg),abs(er$pred_cdc - er$cases),er$diff1 )

pdf <- data.table::melt(data.table::setDT(er),id.vars = c("county_fips_code",
                                                          "month_id"),
                        measure.vars = c("cases",
                                         "jhu_cases",
                                         "pred_cdc",
                                         "avg_pred",
                                         "avg"))



pdf %>% 
  dplyr::filter(county_fips_code %in% w_count) %>%
  ggplot2::ggplot(ggplot2::aes(x = month_id, y = value))+
  ggplot2::geom_line(ggplot2::aes(group = variable,color = variable),alpha = .3)+
  ggplot2::facet_wrap(~county_fips_code)
