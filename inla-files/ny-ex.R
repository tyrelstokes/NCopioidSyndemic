# Housekeeping -----------

`%>%` <- dplyr::`%>%`
`%m+%` <- lubridate::`%m+%`

# Load case data ---------

ny <- read.csv(here::here("state-data/NY_COVID_data.csv"))

# create adj matrix ----------


source(here::here("usa_adj_mat_files/adj-mat-functions.R"))


usa_adj <- read.csv(here::here("usa_adj_mat_files/aug_adj_mat.csv"))


ny_adj_list <-  state_adj_fun(usa_adj = usa_adj,
                               state_abr = "NY",
                               fips = T) 

adj.mat <- ny_adj_list$state_mat

# add remapped ids from 1:n_counties to the state data set ---------------
ny <- map_adj_ids(adj_list = ny_adj_list,
                  state_df = ny)


# join in the population data --------------

source(here::here("census-data/clean-pop-data.R")) # note this step may not work for all states, see source code for explanation

ny$population <- plyr::mapvalues(ny$county_fips_code,
                                 from = pop$fips_num,
                                 to = pop$est_july_2020)

ny$population_scale <- ny$population/(median(ny$cases,na.rm = T)*500)


# Convert the dates to numeric 1:n_months --------

source(here::here("utilities/date-functions.R")) # load mapping function


#ny <- ny %>% dplyr::filter(!(case_month %in% c("2020-01-01",
                                            # "2020-02-01")))

na_df <- ny %>% dplyr::group_by(case_month) %>% dplyr::summarise(na_cases = mean(is.na(cases)),
                                                               na_jhu = mean(is.na(jhu_cases)))


nas <- which((na_df$na_cases == 1) & (na_df$na_jhu ==1) )   
ny <- ny %>% dplyr::filter(!(case_month %in% na_df$case_month[nas]))
ny <- map_month_ids(df = ny,
                    beg = NULL, # default beg and end will be first and last month respectively
                    end = NULL)


#ny <- ny %>% dplyr::filter(month_id != 47)

# Begin the INLA model -----------------

# load libraries and source the helper functions
library(INLA)
source(here::here("inla-files/helper-functions.R"))

# set the dimensions of the flattened data for inla transformation
n <- nrow(ny)
np <- 3

# artificially set the negative counts to 0 (Change this later) 

ny$jhu_cases_reset <- ifelse(ny$jhu_cases < 0, 0,ny$jhu_cases)
ny$jhu_zero <- ifelse(ny$jhu_cases_reset ==0,1,ifelse(is.na(ny$jhu_cases),NA,0))

ny$jhu_cases <- ifelse(is.na(ny$jhu_cases),NA,ny$jhu_cases_reset)


  

# order the dataframe to avoid NA issues in the predict step
  
  ny  <- ny %>% dplyr::mutate(case_na = ifelse(is.na(cases),1,0),
                              jhu_na = ifelse(is.na(jhu_cases),1,0))
  
  ny <- ny %>% dplyr::arrange(case_na,
                              jhu_na,
                              county_fips_code)
    

  

# Create INLA outcome matrix for multidimensional modelling
y <- inla_matrix(n = n,
                 np = np,
                 y_list = list(ny$cases,
                               ny$jhu_zero,
                               ny$jhu_cases))

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
                       x_list = list(ny$population_scale,
                                     ny$population_scale,
                                     ny$population_scale),
                       base_name = "E")


E1 <- E_list$E_1
E2 <- E_list$E_2
E3 <- E_list$E_3



# Create Ids for spatial effects

id_list <- inla_vectors(n = n,
                        np = np,
                        x_list = list(ny$remapped_id,
                                      ny$remapped_id,
                                      ny$remapped_id),
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
                           x_list = list(ny$month_id,
                                         ny$month_id,
                                         ny$month_id),
                           base_name = "month")

m1 <- month_list$month_1
m2 <- month_list$month_2
m3 <- month_list$month_3

# lag cases ----------


# intercepts for time events

ny$after_march_2020 <- ifelse(ny$case_month > "2020-02-01",1,0)

m20_list <-  inla_vectors(n = n,
                          np = np,
                          x_list = list(ny$after_march_2020,
                                        ny$after_march_2020,
                                        ny$after_march_2020),
                          base_name = "march")

mar1 <- m20_list$march_1
mar2 <- m20_list$march_2
mar3 <- m20_list$march_3


ny$after_2021 <- ifelse(ny$case_month > "2020-12-01",1,0)

vac_list <-  inla_vectors(n = n,
                          np = np,
                          x_list = list(ny$after_2021,
                                        ny$after_2021,
                                        ny$after_2021),
                          base_name = "vac")

vac1 <- vac_list$vac_1
vac2 <- vac_list$vac_2
vac3 <- vac_list$vac_3


ny$after_2023 <- ifelse(ny$case_month > "2022-12-01",1,0)

late_list <-  inla_vectors(n = n,
                          np = np,
                          x_list = list(ny$after_2023,
                                        ny$after_2023,
                                        ny$after_2023),
                          base_name = "late")

late1 <- late_list$late_1
late2 <- late_list$late_2
late3 <- late_list$late_3


#lag_list <- inla_vectors(n = n,
                        #   np = np,
                         #  x_list = list(ny$lag_cases,
                          #               ny$lag_cases,
                           #              ny$lag_cases),
                          # base_name = "lag")

#l1 <- lag_list$lag_1
#l2 <- lag_list$lag_2
#l3 <- lag_list$lag_3


# Set the prior for the space and time effects

h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.7, 0.5))) # time effect

#hlist <- list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),  # spatial effect
#             prec.spatial=list(prior="loggamma",param=c(1,0.001)))

hlist = list(theta1 = list("PCprior", c(2, 0.01)),  
             #  Pr(sd<1) = 0.01,# unlikely to have rr>3just based on the spatial confounding
             theta2 = list("PCprior", c(0.5, 0.1))) 
# Model Code

#inla.setOption(mkl=TRUE)
U <- 1
hyper.prec <- list(theta = list(prior = "pc.prec", param = c(U, 0.01)))

beta.prior <- list(prec = list(prior = "normal", param = c(1, 1)),
                   initial = 4, fixed = FALSE)

#betaprec <- .01
#Q.beta = Diagonal(n=1, betaprec)

#X_mat <- 

#l1 <- data.frame(l1) %>% as.matrix()

link_set <- c(rep(1,nrow(ny)),
              rep(2,nrow(ny)),
              rep(3,nrow(ny)))


test_model <- INLA::inla(y ~ -1+ int1 + int2 +int3+
                           mar1 + mar2 + mar3 +
                           vac1 + vac2 + vac3+
                           late1 + late2 + late3+
                           f(id1,
                             E1,
                             model="bym2",
                             graph=as.matrix(adj.mat),
                            # hyper = hlist,
                             group = m1,
                             scale.model=TRUE,
                             adjust.for.con.comp = TRUE,
                             control.group = list(model = 'rw2',
                                                  hyper = hyper.prec,
                                                  scale.model = TRUE),
                             initial = 3)+
                           f(id2,
                             E2,
                             copy = "id1",
                             group = m2,
                             hyper = list(beta = beta.prior))+
                           f(id3,
                             E3,
                             copy = "id1",
                             group = m3,
                             hyper = list(beta = beta.prior)),#+
                          # f(id_slm1,
                            # model = "slm",
                           #  args.slm = list(W = as.matrix(adj.mat),
                                         #    X = l1,
                                         #    rho.min = 0,
                                           #  rho.max = 1,
                                          #   Q.beta = Q.beta)),
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



summary(test_model)

srr <- test_model$summary.random

sr1 <- srr$id1

preds <- test_model$summary.fitted.values
preds$cases <- c(ny$cases,ny$cases,ny$cases)
preds$jhu_cases <- c(ny$jhu_cases,ny$jhu_cases,ny$jhu_cases)
preds$county <- c(ny$county_fips_code,ny$county_fips_code,ny$county_fips_code)
preds$month <- c(ny$month_id,ny$month_id,ny$month_id)

preds_cdc <- preds[1:nrow(ny),]
preds_jhu <- preds[(2*nrow(ny)+1):(3*nrow(ny)),]
preds_zero <- preds[(nrow(ny)+1):(2*nrow(ny)),]

w_count <-  head(ny$county_fips_code)

ny$pred_cdc <- preds_cdc$`0.5quant`
ny$pred_jhu <- preds_jhu$`0.5quant`
ny$pred_zero <- preds_zero$`0.5quant`
ny$jhu_comb <- (1-ny$pred_zero)*ny$pred_jhu
ny$avg_pred <- (ny$jhu_comb+ny$pred_cdc)/2
ny$avg <- (ny$jhu_cases+ny$cases)/2

ny$diff1 <- abs(ny$avg - ny$avg_pred)
ny$diff2 <- ifelse(is.na(ny$avg),abs(ny$pred_cdc - ny$cases),ny$diff1 )

pdf <- data.table::melt(data.table::setDT(ny),id.vars = c("county_fips_code",
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

latent <- test_model$summary.random
