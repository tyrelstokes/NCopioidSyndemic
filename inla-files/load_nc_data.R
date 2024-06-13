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

