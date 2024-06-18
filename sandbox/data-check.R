
adj_mat <- load(here::here("Adj_Matrix_NC.RData"))

nc <- readxl::read_excel(here::here("NCdata.xlsx"))


### Use only 2017 or later data
outcomedataNC = nc[which(nc$year>=2017),]

############################################################
###########        Extract the outcomes         ############
############################################################
y.D = outcomedataNC$`Illicit opioid overdose deaths`
y.E = outcomedataNC$`Drug overdose ED visits`
Pop = outcomedataNC$pop
Year = outcomedataNC$year
y.C = outcomedataNC$AcuteHepCpdf + outcomedataNC$ChronicHepC
y.I = outcomedataNC$HIV
y.T = outcomedataNC$`People served by treatment programs`
y.B = outcomedataNC$`Patients receiving buprenorphine`

############################################################
###########           Address census            ############
############################################################
cens = 1*(is.na(outcomedataNC$`People served by treatment programs`))
cenlb = rep(NA,length(y.T))
cenub = rep(NA,length(y.T))
cenlb[which(cens==1)]= 1
cenub[which(cens==1)]= 5

cenlb[which(is.na(cenlb))]=y.T[which(is.na(cenlb))]
cenub[which(is.na(cenub))]=y.T[which(is.na(cenub))]

cen=cens
bd=cbind(cenlb,cenub)

### plug in values for censored 
y.T0 = y.T
y.T0[which(cens == 1)] = 5

############################################################
###########    Load the adjacency matrix        ############
############################################################

#load("Adj_Matrix_NC.RData")
adj.mat = adj.mat
n<-nrow(adj.mat)
TT = length(unique(outcomedataNC$year))
num<-colSums(adj.mat)


##############################################################################
####### Standardize the expected count to the state level rate in 2017 ####### 
##############################################################################

ET = rep(0,n*TT)
EB = rep(0,n*TT)
ED = rep(0,n*TT)
EE = rep(0,n*TT)
EC = rep(0,n*TT)
EI = rep(0,n*TT)
for(t in 1:TT){
  for(i in 1:n){
    ET[((t-1)*n+i)] =  sum(y.T0[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    ED[((t-1)*n+i)] =  sum(y.D[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EB[((t-1)*n+i)] =  sum(y.B[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EE[((t-1)*n+i)] =  sum(y.E[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EC[((t-1)*n+i)] =  sum(y.C[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
    EI[((t-1)*n+i)] =  sum(y.I[which(outcomedataNC$year==2017)])/sum(Pop[which(outcomedataNC$year==2017)])*Pop[((t-1)*n+i)]
  }
}


# Simple INLA model ----------------

library(INLA)

counties <- unique(nc$county)

adj_id <- as.numeric(plyr::mapvalues(outcomedataNC$county, from = counties,
                          to = c(1:100)))


yrs <- unique(outcomedataNC$year)
year_id <- as.numeric(plyr::mapvalues(outcomedataNC$year,
                                     from = yrs, to = c(1:5)))



ldf <- data.frame(y.T = y.T, ET = ET,
                  id = adj_id,
                  year = year_id)


  
## spefiying the priors for the unstri and str 

h.spec <- list(prec = list(prior = 'pc.prec', param = c(0.03, 0.03)))

hlist <- list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)), 
              prec.spatial=list(prior="loggamma",param=c(1,0.001)))

m1 <- INLA::inla(y.T ~ f(id,ET,model="bym",
                         graph=adj.mat, scale.model=TRUE,
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


## Stacked approach

np <- 2
n <- nrow(outcomedataNC)

y <- matrix(ncol = np, nrow = n*np)
ET1 <- vector(length = n*np)
EB1 <- vector(length = n*np)
#int1 <- matrix(ncol = np, nrow = n*np)
#int2 <- matrix(ncol = np, nrow = n*np)

int1 <- vector(length = n*np)
int2 <- vector(length = n*np)

#id1 <- matrix(ncol = np, nrow = n*np)
#id2 <- matrix(ncol = np, nrow = n*np)
id1 <- vector(length = n*np)
id2 <-  vector(length = n*np)


y[1:n,1] <- y.T
y[(n+1):(2*n),2] <- y.E


ET1[1:n] <- ET
ET1[(n+1):(2*n)] <- rep(NA,n)

EB1[(n+1):(2*n)] <- EB
EB1[1:n] <- rep(NA,n)


int1[1:n] <- rep(1,n)
int1[(n+1):(2*n)] <- rep(NA,n)

int2[1:n] <- rep(NA,n)
int2[(n+1):(2*n)] <- rep(1,n)

#id1[1:n,1] <- adj_id
#id2[(n+1):(2*n),2] <- adj_id

id1[1:n] <- adj_id
id1[(n+1):(2*n)] <- rep(NA,n)
id2[(n+1):(2*n)] <- adj_id
id2[1:n] <- rep(NA,n)



y1 <- vector(length = n*np)
y2 <- vector(length = n*np)




####







m1 <- INLA::inla(y ~ -1+ int1 + int2 + f(id1,ET1,model="bym",
                         graph=adj.mat, scale.model=TRUE,
                         hyper = hlist)+
                   f(id2,EB1,model="bym",
                     graph=adj.mat, scale.model=TRUE,
                     hyper = hlist),
                 family=c("poisson","poisson"),
                 data=list(y =y,
                           EB1 = EB1,
                           ET1 = ET1,
                           int1 = int1,
                           int2 = int2,
                           id1 = id1,
                           id2 = id2),
                 control.compute=list(waic=TRUE),
                 control.predictor = list(compute = TRUE,
                                          link =1))



summary(m1)












stack1 <- inla.stack(
  data = list(y.T = y.T),
  A = list(adj.mat),
  effects = list(ET = ET,
                 id.t = adj_id),
  tag = "T"
)


