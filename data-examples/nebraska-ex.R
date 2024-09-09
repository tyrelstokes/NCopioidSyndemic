# Housekeeping -----------

`%>%` <- dplyr::`%>%`
`%m+%` <- lubridate::`%m+%`

# Load case data ---------

ne <- read.csv(here::here("state-data/NE_COVID_data.csv"))

# create adj matrix ----------


source(here::here("usa_adj_mat_files/adj-mat-functions.R"))


usa_adj <- read.csv(here::here("usa_adj_mat_files/aug_adj_mat.csv"))


neb_adj_list <-  state_adj_fun(usa_adj = usa_adj,
                                 state_abr = "NE",
                                 fips = T) 

adj.mat <- neb_adj_list$state_mat

# add remapped ids from 1:n_counties to the state data set ---------------
ne <- map_adj_ids(adj_list = neb_adj_list,
                  state_df = ne)


# join in the population data --------------

source(here::here("census-data/clean-pop-data.R")) # note this step may not work for all states, see source code for explanation

ne$population <- plyr::mapvalues(ne$county_fips_code,
                                 from = pop$fips_num,
                                 to = pop$est_july_2020)


# Convert the dates to numeric 1:n_months --------

source(here::here("utilities/date-functions.R")) # load mapping function


ne <- map_month_ids(df = ne,
                    beg = NULL, # default beg and end will be first and last month respectively
                    end = NULL)

# Begin the INLA model -----------------

# load libraries and source the helper functions
library(INLA)
source(here::here("inla-files/helper-functions.R"))

# set the dimensions of the flattened data for inla transformation
n <- nrow(ne)
np <- 2

# artificially set the negative counts to NA (Change this later) 

ne$jhu_cases <- ifelse(ne$jhu_cases < 0, NA,ne$jhu_cases)

# Create INLA outcome matrix for multidimensional modelling
y <- inla_matrix(n = n,
                 np = np,
                 y_list = list(ne$cases,ne$jhu_cases))

# Create proper intercepts 
int_list <- inla_vectors(n =n,
                         np = np,
                         x_list = list(rep(1,n),rep(1,n)),
                         base_name = "int")

int1 <- int_list$int_1
int2 <- int_list$int_2

# Create the E vectors for rate modelling
E_list <- inla_vectors(n = n,
                       np = np,
                       x_list = list(ne$population,ne$population),
                       base_name = "E")


E1 <- E_list$E_1
E2 <- E_list$E_2

# Create Ids for spatial effects

id_list <- inla_vectors(n = n,
                        np = np,
                        x_list = list(ne$remapped_id,ne$remapped_id),
                        base_name = "id")


id1 <- id_list$id_1
id2 <- id_list$id_2

# Create the month indices

month_list <- inla_vectors(n = n,
                          np = np,
                          x_list = list(ne$month_id,ne$month_id),
                          base_name = "month")

m1 <- month_list$month_1
m2 <- month_list$month_2


# Set the prior for the space and time effects

h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.7, 0.5))) # time effect

hlist <- list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)),  # spatial effect
              prec.spatial=list(prior="loggamma",param=c(1,0.001)))

hlist = list(theta1 = list("PCprior", c(2, 0.01)),  
             # Pr(sd<1) = 0.01, unlikely to have rr>3just based on the spatial confounding
             theta2 = list("PCprior", c(0.5, 0.75))) 
# Model Code

#inla.setOption(mkl=TRUE)

test_model <- INLA::inla(y ~ -1+ int1 + int2 + f(id1,E1,model="bym",
                                         graph=as.matrix(adj.mat),
                                         hyper = hlist,
                                         group = m1,
                                         control.group = list(model = 'ar1', hyper = h.spec),
                                         initial = 3)+
                           f(id2,E2,copy = "id1",
                             group = m2,
                             hyper = list(beta = list(fixed = FALSE))),
                 family=c("poisson","zeroinflatedpoisson0"),
                 data=list(y =y,
                           E1 = E1,
                           E2 = E2,
                           int1 = int1,
                           int2 = int2,
                           id1 = id1,
                           id2 = id2),
                 control.compute=list(waic=TRUE),
                 control.predictor = list(compute = TRUE,
                                          link =1),
                 verbose = TRUE)



summary(test_model)


preds <- test_model$summary.fitted.values

latent <- test_model$summary.random
