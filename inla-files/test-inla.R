# Source Necessary Files

source(here::here("inla-files/load_nc_data.R"))
source(here::here("inla-files/helper-functions.R"))



ldf <- data.frame(y.T = y.T, ET = ET,
                  id = adj_id,
                  year = year_id)



## spefiying the priors for the unstri and str 

h.spec <- list(prec = list(prior = 'pc.cor0', param = c(0.03, 0.03)))

hlist <- list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)), 
              prec.spatial=list(prior="loggamma",param=c(1,0.001)))

m1 <- INLA::inla(y.T ~ f(id,
                         ET,
                         model="bym",
                         graph=adj.mat,
                         scale.model=TRUE,
                         hyper = hlist,
                         group = year,
                         control.group = list(model = 'iid', hyper = h.spec,
                                              scale.model = TRUE)),
                 family="poisson",
                 data=ldf,
                 control.compute=list(waic=TRUE),
                 control.predictor = list(compute = TRUE,
                                          link =1))


summary(m1)
sfit <- m1$summary.fitted.values


preds <- sfit$mean


pred_df <- data.frame(preds = preds, y = y.T)


# Try with multiple likelihoods -------

# Set some of the parameters
np <- 2 # dimension of the outcome
n <- nrow(outcomedataNC)


# Create INLA outcome
y <- inla_matrix(n = n,
                 np = np,
                 y_list = list(y.T,y.E))



# Create proper intercepts
int_list <- inla_vectors(n =n,
                         np = np,
                         x_list = list(rep(1,n),rep(1,n)),
                         base_name = "int")

int1 <- int_list$int_1
int2 <- int_list$int_2

# Create the E vectors
E_list <- inla_vectors(n = n,
                      np = np,
                      x_list = list(ET,EE),
                      base_name = "E")


ET1 <- E_list$E_1
EE1 <- E_list$E_2

# Create Ids for spatial effects

id_list <- inla_vectors(n = n,
                        np = np,
                        x_list = list(adj_id,adj_id),
                        base_name = "id")


id1 <- id_list$id_1
id2 <- id_list$id_2
id3 <- id2

# Create the year indices

year_list <- inla_vectors(n = n,
                          np = np,
                          x_list = list(year_id,year_id),
                          base_name = "year")

y1 <- year_list$year_1
y2 <- year_list$year_2




# Set the prior for the time dimension

h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.5)))


# Model Code



m2 <- INLA::inla(y ~ -1+ int1 + int2 + f(id1,ET1,model="bym",
                                         graph=adj.mat, scale.model=TRUE,
                                         hyper = hlist,
                                         group = y1,
                                         control.group = list(model = 'ar1', hyper = h.spec,
                                                              scale.model = TRUE))+
                   f(id2,EE1,model="bym",
                     graph=adj.mat, scale.model=TRUE,
                     hyper = hlist,
                     group = y2,
                     control.group = list(model = 'ar1', hyper = h.spec,
                                          scale.model = TRUE))+
                   f(id3,EE1,copy = "id1",
                     group = y2,
                     hyper = list(beta = list(fixed = FALSE))),
                 family=c("poisson","poisson"),
                 data=list(y =y,
                           EE1 = EE1,
                           ET1 = ET1,
                           int1 = int1,
                           int2 = int2,
                           id1 = id1,
                           id2 = id2,
                           id3 = id3),
                 control.compute=list(waic=TRUE),
                 control.predictor = list(compute = TRUE,
                                          link =1))



summary(m2)


# Extend to 6 outcomes ----------------

# Set some of the parameters
np <- 6
n <- nrow(outcomedataNC)


# Create INLA outcome
y <- inla_matrix(n = n,
                 np = np,
                 y_list = list(y.T,y.E,y.I,y.D,y.C,y.B))



# Create proper intercepts
int_list <- inla_vectors(n =n,
                         np = np,
                         x_list = list(rep(1,n),rep(1,n),rep(1,n),rep(1,n),rep(1,n),rep(1,n)),
                         base_name = "int")

int1 <- int_list$int_1
int2 <- int_list$int_2
int3 <- int_list$int_3
int4 <- int_list$int_4
int5 <- int_list$int_5
int6 <- int_list$int_6

# Create the E vectors
E_list <- inla_vectors(n = n,
                       np = np,
                       x_list = list(ET,EE,EI,ED,EC,EB),
                       base_name = "E")


ET1 <- E_list$E_1
EE1 <- E_list$E_2
EI1 <- E_list$E_3
ED1  <- E_list$E_4
EC1 <- E_list$E_5
EB1 <- E_list$E_6

# Create Ids for spatial effects

id_list <- inla_vectors(n = n,
                        np = np,
                        x_list = list(adj_id,
                                      adj_id,
                                      adj_id,
                                      adj_id,
                                      adj_id,
                                      adj_id),
                        base_name = "id")


idt <- id_list$id_1
ide <- id_list$id_2
idi <- id_list$id_3
idd <- id_list$id_4
idc <- id_list$id_5
idb <- id_list$id_6


id_cs_t <- idt
id_cs_e <- ide
id_cs_i <- idi
id_cs_d <- idd
id_cs_c <- idc
id_cs_b <- idb

id_ct_b <- idb

id_cd_e <- ide

id_cc_i <- idi


# Create the year indices

year_list <- inla_vectors(n = n,
                          np = np,
                          x_list = list(year_id,
                                        year_id,
                                        year_id,
                                        year_id,
                                        year_id,
                                        year_id),
                          base_name = "year")

yt <- year_list$year_1
ye <- year_list$year_2
yi <- year_list$year_3

yd <- year_list$year_4
yc <- year_list$year_5
yb <- year_list$year_6


# Run the model

hyper_copy <- list(beta = list(prior = 'normal', param = c(0.5,.1 )))
h.spec <- list(rho = list(prior = 'pc.cor0', param = c(0.5, 0.8)))

hlist3 <- list(theta1=list(prior="pc.prec",param=c(0.5,0.01)), 
              theta2=list(prior="pc.prec",param=c(0.5,0.5)))

hlist2 <- list(prec=list(prior="pc.prec",param=c(0.5,0.5)))

m3 <- INLA::inla(y ~ -1+ int1 + int2 + int3 + int4 + int5 + int6+
                   f(id_cs_t,ET1,model="bym2",
                    graph=adj.mat, 
                   # scale.model=TRUE,
                    hyper = hlist3,
                    group = yt,
                    constr = TRUE,
                    control.group = list(model = 'ar1', hyper = h.spec,
                                         scale.model = TRUE))+
                   
                   f(id_cs_e,EE1,copy = "id_cs_t",group = ye, hyper = hyper_copy) +
                   f(id_cs_i,EI1,copy = "id_cs_t",group = yi, hyper = hyper_copy)+
                   f(id_cs_d,ED1,copy = "id_cs_t",group = yd, hyper = hyper_copy)+
                   f(id_cs_c,EC1,copy = "id_cs_t",group = yc, hyper = hyper_copy)+
                   f(id_cs_b,EB1,copy = "id_cs_t",group = yb, hyper = hyper_copy)+
                   
                   f(id_ct_b,EB1,model="bym2",
                    graph=adj.mat,
                    hyper = hlist3,
                   # scale.model=TRUE,
                    group = yb,
                    constr = TRUE,
                    control.group = list(model = 'ar1', hyper = h.spec,
                                         scale.model = TRUE))+
                   
                    f(idt,ET1,copy = "id_ct_b",group = yt, hyper = hyper_copy)+
                   f(idd,ED1,model="bym2",
                    graph=adj.mat,
                    hyper = hlist3,
                    group = yd,
                   # scale.model=TRUE,
                    constr = TRUE,
                    control.group = list(model = 'ar1', hyper = h.spec,
                         scale.model = TRUE))+
                     f(id_cd_e,EE1,copy = "idd",group = ye, hyper = hyper_copy)+
                   f(idc,EC1,model="bym2",
                     graph=adj.mat,
                     hyper = hlist3,
                     group = yc,
                    # scale.model=TRUE,
                     constr = TRUE,
                     control.group = list(model = 'ar1', hyper = h.spec,
                                          scale.model = TRUE))+
                  f(id_cc_i,EI1,copy = "idc",group = yi, hyper = hyper_copy),
                 family=c("poisson",
                          "poisson",
                          "poisson",
                          "poisson",
                          "poisson",
                          "poisson"),
                 data=list(y =y,
                           EE1 = EE1,
                           ET1 = ET1,
                           EB1 = EB1,
                           EC1 = EC1,
                           ED1 = ED1,
                           EI1 = EI1,
                           int1 = int1,
                           int2 = int2,
                           int3 = int3,
                           int4 = int4,
                           int5 = int5,
                           int6 = int6,
                           idt = idt,
                           ide = ide,
                           idi =idi,
                           idd =idd,
                           idc =idc,
                           idb =idb,
                           id_cs_t =id_cs_t,
                           id_cs_e = id_cs_e,
                           id_cs_i = id_cs_i,
                           id_cs_d = id_cs_d,
                           id_cs_c = id_cs_c,
                           id_cs_b = id_cs_b,
                           id_ct_b = id_ct_b,
                           id_cd_e = id_cd_e,
                           id_cc_i = id_cc_i),
                 control.compute=list(waic=TRUE),
                 control.predictor = list(compute = TRUE,
                                          link =1),
                 verbose = FALSE)


summary(m3)


## Besag version -------------

m4 <- INLA::inla(y ~ -1+ int1 + int2 + int3 + int4 + int5 + int6+
                   f(id_cs_t,ET1,model="besag",
                     graph=adj.mat,
                     scale.model=TRUE,
                    # hyper = hlist3,
                     group = yt,
                     control.group = list(model = 'ar1', hyper = h.spec,
                                          scale.model = TRUE))+
                   
                   f(id_cs_e,EE1,copy = "id_cs_t",group = ye, hyper = hyper_copy) +
                   f(id_cs_i,EI1,copy = "id_cs_t",group = yi, hyper = hyper_copy)+
                   f(id_cs_d,ED1,copy = "id_cs_t",group = yd, hyper = hyper_copy)+
                   f(id_cs_c,EC1,copy = "id_cs_t",group = yc, hyper = hyper_copy)+
                   f(id_cs_b,EB1,copy = "id_cs_t",group = yb, hyper = hyper_copy)+
                   
                   f(id_ct_b,EB1,model="besag",
                     graph=adj.mat,
                     scale.model=TRUE,
                    # hyper = hlist3,
                     group = yb,
                     control.group = list(model = 'ar1', hyper = h.spec,
                                          scale.model = TRUE))+
                   
                   f(idt,ET1,copy = "id_ct_b",group = yt, hyper = hyper_copy)+
                   f(idd,ED1,model="besag",
                     scale.model=TRUE,
                     graph=adj.mat,
                     #hyper = hlist3,
                     group = yd,
                     control.group = list(model = 'ar1', hyper = h.spec,
                                          scale.model = TRUE))+
                   f(id_cd_e,EE1,copy = "idd",group = ye, hyper = hyper_copy)+
                   f(idc,EC1,model="besag",
                     graph=adj.mat,
                     scale.model=TRUE,
                    # hyper = hlist3,
                     group = yc,
                     control.group = list(model = 'ar1', hyper = h.spec,
                                          scale.model = TRUE))+
                   f(id_cc_i,EI1,copy = "idc",group = yi, hyper = hyper_copy),
                 family=c("poisson",
                          "poisson",
                          "poisson",
                          "poisson",
                          "poisson",
                          "poisson"),
                 data=list(y =y,
                           EE1 = EE1,
                           ET1 = ET1,
                           EB1 = EB1,
                           EC1 = EC1,
                           ED1 = ED1,
                           EI1 = EI1,
                           int1 = int1,
                           int2 = int2,
                           int3 = int3,
                           int4 = int4,
                           int5 = int5,
                           int6 = int6,
                           idt = idt,
                           ide = ide,
                           idi =idi,
                           idd =idd,
                           idc =idc,
                           idb =idb,
                           id_cs_t =id_cs_t,
                           id_cs_e = id_cs_e,
                           id_cs_i = id_cs_i,
                           id_cs_d = id_cs_d,
                           id_cs_c = id_cs_c,
                           id_cs_b = id_cs_b,
                           id_ct_b = id_ct_b,
                           id_cd_e = id_cd_e,
                           id_cc_i = id_cc_i),
                 control.compute=list(waic=TRUE),
                 control.predictor = list(compute = TRUE,
                                          link =1),
                 verbose = FALSE)


summary(m4)

