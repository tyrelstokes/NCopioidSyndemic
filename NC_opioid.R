library(nimble)
library(coda)
############################################################
###########             Load data               ############
############################################################
library(readxl)
NCdata <- read_excel("NCdata.xlsx")
### Use only 2017 or later data
outcomedataNC = NCdata[which(NCdata$year>=2017),]

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

load("Adj_Matrix_NC.RData")
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


############################################################
###########         Set up for NIMBLE          ############
############################################################


adj<-NULL
for(j in 1:n){
  adj<-c(adj,which(adj.mat[j,]==1)) ##takes only the neighbors
}
adj<-as.vector(adj)
num<-as.vector(num)
weights<-1+0*adj

## Initialize the loadings matrix
load.mat.transp = matrix(c(1,0,0,0,0,0, 0,1,0,0,0,0, 0,0,1,0,0,0, 0,0,0,1,0,0),nrow= 4, ncol = 6, byrow = T)

# Function to compute the sign of a number
signVec.nim <- nimbleFunction(
  run = function(x = double(1)) {
    returnType(double(1))
    # Initialize an empty vector to store the results
    result <- numeric(length(x))
    
    # Loop through each element of x
    for (i in seq_along(x)) {
      # Check the sign of the current element
      if (x[i] > 0) {
        result[i] <- 1
      } else if (x[i] < 0) {
        result[i] <- -1
      } else {
        result[i] <- 0
      }
    }
    
    # Return the resulting vector
    return(result)
  })

## Function to perform the QR factorization on the transpose of load.mat
qr.decomp.transp<-nimbleFunction(
  run = function(X = double(2)) {
    returnType(double(2))
    #Xt = t(X)
    m<-dim(X)[1]
    n<-dim(X)[2]
    QR<- matrix(0, m, m+n)
    R.full = matrix(0, m,n)
    X1 = X[1:m, 1:m]
    X2 = X[1:m, (m+1):n]
    m1 = dim(X1)[1]
    n1 = dim(X1)[2]
    # initialize Q and R
    Q<-matrix(0,m1,n1)
    R<-matrix(0,n1,n1)
    u1<-X1[1:m1,1]
    Q[1:m1,1]<- u1[1:m1]/sqrt(inprod(u1[1:m1], u1[1:m1]))
    R[1,1]<- inprod(X1[1:m1,1], Q[1:m1, 1])
    R[1,2]<- inprod(X1[1:m1,2], Q[1:m1, 1])
    R[1,3]<- inprod(X1[1:m1,3], Q[1:m1, 1])
    R[1,4]<- inprod(X1[1:m1,4], Q[1:m1, 1])
    u2<- X1[1:m1,2] - R[1,2] * Q[1:m1, 1]
    Q[1:m1, 2]<- u2[1:m1]/sqrt(inprod(u2[1:m1],u2[1:m1]))
    R[2,2]<- inprod(X1[1:m1, 2], Q[1:m1, 2])
    R[2,3]<- inprod(X1[1:m1, 3], Q[1:m1, 2])
    R[2,4]<- inprod(X1[1:m1, 4], Q[1:m1, 2])
    u3<- X1[1:m1,3] - R[1,3] * Q[1:m1, 1] - R[2,3] * Q[1:m1, 2]
    Q[1:m1, 3]<- u3[1:m1]/sqrt(inprod(u3[1:m1],u3[1:m1]))
    R[3,3]<- inprod(X1[1:m1, 3], Q[1:m1, 3])
    R[3,4]<- inprod(X1[1:m1, 4], Q[1:m1, 3])
    u4<- X1[1:m1, 4] - R[1,4] * Q[1:m1, 1] - R[2,4] * Q[1:m1, 2] - R[3,4]*Q[1:m1, 3]
    Q[1:m1, 4]<- u4[1:m1]/sqrt(inprod(u4[1:m1], u4[1:m1]))
    R[4,4]<- inprod(X1[1:m1, 4], Q[1:m1, 4])
    sign.vect <- signVec.nim(diag(R[1:n1, 1:n1]))
    Q1 <- Q[1:m1, 1:n1] %*% diag(sign.vect[1:n1])
    R1 <- diag(sign.vect[1:n1]) %*% R[1:n1, 1:n1]
    R.full[1:m, 1:m] <- R1[1:n1, 1:n1]
    R2 <- t(Q1[1:m1, 1:n1]) %*% X2[1:m, 1:(n-m)]
    R.full[1:m, (m+1):n] <- R2[1:m, 1:(n-m)]
    #R.full[1:m, 1:n] = R.full[1:m, 1:n]%*%P[1:6, 1:6]
    
    ## Bind Q and R together in a big matrix
    QR[1:m, 1:m]<- Q1[1:m, 1:m]
    QR[1:m, (m+1):(n+m)]<- R.full[1:m, 1:n]
    return(QR)
  })

model_code=nimbleCode({
  for (t in T0:T0){
    for (i in 1:n){
      cen[n*(t-T0)+i] ~ dinterval(y.T[n*(t-T0)+i], bd[n*(t-T0)+i,1:2])
      y.D[n*(t-T0)+i] ~ dpois(ED[(n*(t-T0)+i)]*lambdaD[n*(t-T0)+i])
      y.E[n*(t-T0)+i] ~ dpois(EE[(n*(t-T0)+i)]*lambdaE[n*(t-T0)+i])
      y.T[n*(t-T0)+i] ~ dpois(ET[(n*(t-T0)+i)]*lambdaT[n*(t-T0)+i])
      y.C[n*(t-T0)+i] ~ dpois(EC[(n*(t-T0)+i)]*lambdaC[n*(t-T0)+i])
      y.B[n*(t-T0)+i] ~ dpois(EB[(n*(t-T0)+i)]*lambdaB[n*(t-T0)+i])
      y.I[n*(t-T0)+i] ~ dpois(EI[(n*(t-T0)+i)]*lambdaI[n*(t-T0)+i])
      
      U1[n*(t-T0)+i, 1:F]<- (U[n*(t-T0)+i,1:F]+mu[n*(t-T0)+i,1:F]) %*% t(Q.A[1:4, 1:4])
      
      log(lambdaD[n*(t-T0)+i]) <- inprod(L[1,1:F],U1[n*(t-T0)+i, 1:F]) + VD[(n*(t-T0)+i)]
      log(lambdaE[n*(t-T0)+i]) <- inprod(L[2,1:F],U1[n*(t-T0)+i, 1:F]) + VE[(n*(t-T0)+i)]
      log(lambdaT[n*(t-T0)+i]) <- inprod(L[3,1:F],U1[n*(t-T0)+i, 1:F]) + VT[(n*(t-T0)+i)]
      log(lambdaC[n*(t-T0)+i]) <- inprod(L[4,1:F],U1[n*(t-T0)+i, 1:F]) + VC[(n*(t-T0)+i)]
      log(lambdaB[n*(t-T0)+i]) <- inprod(L[5,1:F],U1[n*(t-T0)+i, 1:F]) + VB[(n*(t-T0)+i)]
      log(lambdaI[n*(t-T0)+i]) <- inprod(L[6,1:F],U1[n*(t-T0)+i, 1:F]) + VI[(n*(t-T0)+i)]
      
      
      VD[n*(t-T0)+i] ~ dnorm(0,tau.VD)
      VE[n*(t-T0)+i] ~ dnorm(0,tau.VE)
      VT[n*(t-T0)+i] ~ dnorm(0,tau.VT)
      VC[n*(t-T0)+i] ~ dnorm(0,tau.VC)
      VB[n*(t-T0)+i] ~ dnorm(0,tau.VB) 
      VI[n*(t-T0)+i] ~ dnorm(0,tau.VI) 
      mu[n*(t-T0)+i,1:F] <- mu0[t-T0+1,1:F]
    }
    
    
    
    # ICAR prior for the spatial factors
    U[(n*(t-T0)+1):(n*(t-T0)+n),1] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),1], zero_mean=1)
    U[(n*(t-T0)+1):(n*(t-T0)+n),2] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),2], zero_mean=1)
    U[(n*(t-T0)+1):(n*(t-T0)+n),3] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),3], zero_mean=1)
    U[(n*(t-T0)+1):(n*(t-T0)+n),4] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),4], zero_mean=1)
    
    
    
    tau.U[t-(T0-1),1] ~ dgamma(0.5,0.5)
    tau.U[t-(T0-1),2] ~ dgamma(0.5,0.5)
    tau.U[t-(T0-1),3] ~ dgamma(0.5,0.5)
    tau.U[t-(T0-1),4] ~ dgamma(0.5,0.5)
    
  }  
  
  for (t in (T0+1):(T0+TT-1)){
    for (i in 1:n){
      cen[n*(t-T0)+i] ~ dinterval(y.T[n*(t-T0)+i],bd[n*(t-T0)+i,1:2])
      y.D[n*(t-T0)+i] ~ dpois(ED[(n*(t-T0)+i)]*lambdaD[n*(t-T0)+i])
      y.E[n*(t-T0)+i] ~ dpois(EE[(n*(t-T0)+i)]*lambdaE[n*(t-T0)+i])
      y.T[n*(t-T0)+i] ~ dpois(ET[(n*(t-T0)+i)]*lambdaT[n*(t-T0)+i])
      y.C[n*(t-T0)+i] ~ dpois(EC[(n*(t-T0)+i)]*lambdaC[n*(t-T0)+i])
      y.B[n*(t-T0)+i] ~ dpois(EB[(n*(t-T0)+i)]*lambdaB[n*(t-T0)+i])
      y.I[n*(t-T0)+i] ~ dpois(EI[(n*(t-T0)+i)]*lambdaI[n*(t-T0)+i])
      
      U1[n*(t-T0)+i, 1:F]<- (U[n*(t-T0)+i,1:F]+mu[n*(t-T0)+i,1:F]) %*% t(Q.A[1:4, 1:4])
      #U2[n*(t-T0)+i, 1:F]<- t(U1[1:F, n*(t-T0)+i])
      
      log(lambdaD[n*(t-T0)+i]) <- inprod(L[1,1:F],U1[n*(t-T0)+i, 1:F]) + VD[(n*(t-T0)+i)]
      log(lambdaE[n*(t-T0)+i]) <- inprod(L[2,1:F],U1[n*(t-T0)+i, 1:F]) + VE[(n*(t-T0)+i)]
      log(lambdaT[n*(t-T0)+i]) <- inprod(L[3,1:F],U1[n*(t-T0)+i, 1:F]) + VT[(n*(t-T0)+i)]
      log(lambdaC[n*(t-T0)+i]) <- inprod(L[4,1:F],U1[n*(t-T0)+i, 1:F]) + VC[(n*(t-T0)+i)]
      log(lambdaB[n*(t-T0)+i]) <- inprod(L[5,1:F],U1[n*(t-T0)+i, 1:F]) + VB[(n*(t-T0)+i)]
      log(lambdaI[n*(t-T0)+i]) <- inprod(L[6,1:F],U1[n*(t-T0)+i, 1:F]) + VI[(n*(t-T0)+i)]
      
      VD[n*(t-T0)+i] ~ dnorm(0,tau.VD)
      VE[n*(t-T0)+i] ~ dnorm(0,tau.VE)
      VT[n*(t-T0)+i] ~ dnorm(0,tau.VT)
      VC[n*(t-T0)+i] ~ dnorm(0,tau.VC)
      VB[n*(t-T0)+i] ~ dnorm(0,tau.VB) 
      VI[n*(t-T0)+i] ~ dnorm(0,tau.VI) 
      mu[n*(t-T0)+i,1:F] <- mu0[t-T0+1,1:F]
    }
    
    
    # ICAR prior for the spatial factors
    U[(n*(t-T0)+1):(n*(t-T0)+n),1] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),1], zero_mean=1)
    U[(n*(t-T0)+1):(n*(t-T0)+n),2] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),2], zero_mean=1)
    U[(n*(t-T0)+1):(n*(t-T0)+n),3] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),3], zero_mean=1)
    U[(n*(t-T0)+1):(n*(t-T0)+n),4] ~ dcar_normal(adj[], weights[], num[], tau.U[t-(T0-1),4], zero_mean=1)
    
    tau.U[t-(T0-1),1] ~ dgamma(0.5,0.5)
    tau.U[t-(T0-1),2] ~ dgamma(0.5,0.5)
    tau.U[t-(T0-1),3] ~ dgamma(0.5,0.5)
    tau.U[t-(T0-1),4] ~ dgamma(0.5,0.5)
    
  }
  
  #Priors
  
  ### loadings matrix entries
  load.mat.transp[1,2] ~ dnorm(0,sd=10)
  load.mat.transp[1,3] ~ dnorm(0,sd=10)
  load.mat.transp[1,4] ~ dnorm(0,sd=10)
  load.mat.transp[1,5] ~ dnorm(0,sd=10)
  load.mat.transp[1,6] ~ dnorm(0,sd=10)
  load.mat.transp[2,1] ~ dnorm(0,sd=10)
  load.mat.transp[3,5] ~ dnorm(0,sd=10)
  load.mat.transp[4,6] ~ dnorm(0,sd=10)
  
  #Compute L and Q from the LQ factorization
  QR.A[1:4, 1:10] <- qr.decomp.transp(X = load.mat.transp[1:4,1:6])
  Q.A[1:4, 1:4] <- t(QR.A[1:4, 1:4])
  L[1:6, 1:4] <- t(QR.A[1:4, 5:10])
  
  
  for(t in T0:(T0+TT-1)){
    mu0[t-T0+1,1] ~ dflat()
    mu0[t-T0+1,2] ~ dflat()
    mu0[t-T0+1,3] ~ dflat()
    mu0[t-T0+1,4] ~ dflat()
  }
  
  tau.VD ~ dgamma(0.5,0.5)
  tau.VE ~ dgamma(0.5,0.5)
  tau.VT ~ dgamma(0.5,0.5)
  tau.VB ~ dgamma(0.5,0.5)
  tau.VC ~ dgamma(0.5,0.5)
  tau.VI ~ dgamma(0.5,0.5)
  eta[1] ~ dunif(0,1)
  eta[2] ~ dunif(0,1)
  eta[3] ~ dunif(0,1)
  eta[4] ~ dunif(0,1)
  
})

############################################################
#########              Call NIMBLE                 ###########
############################################################

mod_constants=list(bd=bd,n=n,TT=TT,T0=min(Year),num=num,adj=adj,weights=weights,F=4)

mod_data=list(cen=cen,y.D=y.D,y.E=y.E,y.T=y.T,y.B =y.B,y.C=y.C,y.I=y.I,ED=ED,EE=EE,ET=ET,EB=EB,EC=EC,EI=EI)

mod_inits=list(U=matrix(0,n*TT,4),U1=matrix(0,n*TT,4),load.mat.transp=load.mat.transp,L=matrix(0,6,4),Q.A=matrix(0,4,4),QR.A=matrix(0,4,10),mu0=matrix(0,TT,4),VD=rep(0,n*TT),VE=rep(0,n*TT),VT=rep(0,n*TT),VB=rep(0,n*TT),VC=rep(0,n*TT),VI=rep(0,n*TT),y.T=floor(.5*(bd[,1]+bd[,2])))

# Build the model.
nim_model <- nimbleModel(model_code,mod_constants,mod_data,mod_inits)
compiled_model <- compileNimble(nim_model,resetFunctions = T)



# Set up samplers.
mcmc_conf <- configureMCMC(nim_model,monitors=c("mu0","mu","L","Q.A","load.mat.transp",
                                                "U","U1","eta","VD","VE","VT","VB","VC","VI","lambdaD","lambdaE","lambdaT","lambdaB","lambdaC","lambdaI"),useConjugacy = TRUE)

mod_mcmc<-buildMCMC(mcmc_conf)
compiled_mcmc<-compileNimble(mod_mcmc, project = nim_model,resetFunctions = T)

# Run the model 
MCS=500000
st<-Sys.time()

samples_randomLQ=runMCMC(compiled_mcmc,inits=mod_inits,
                         nchains = 1, nburnin=floor(MCS/2),niter = MCS,samplesAsCodaMCMC = TRUE,thin=50,
                         summary = FALSE, WAIC = FALSE,progressBar=TRUE) 
Sys.time()-st
save(samples_randomLQ,file="ModOutput.Rda")

quit()


nc$adj <- adj


