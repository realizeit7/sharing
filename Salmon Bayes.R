library(rjags)

pinkdata <- read.csv("~/Downloads/pinkdata.csv")
pinkdata$stock2<-as.numeric(pinkdata$Stock)
pinkdata$stock2<-ifelse(pinkdata$stock2==1,100,pinkdata$stock2)
pinkdata$stock2<-ifelse(pinkdata$stock2<5,pinkdata$stock2-1,pinkdata$stock2)
pinkdata$stock2<-ifelse(pinkdata$stock2==100,4,pinkdata$stock2)
pinkdata$Region2<-as.numeric(pinkdata$Region)


pinkdataodds<-subset(pinkdata,pinkdata$BY%%2==1)
Odds<-c()
Ostock<-c()
ORegion<-c()
Oddr<-c()
OddT<-c()
for ( i in 2: 648){
  Odds[i-1]<-ifelse(pinkdataodds$BY[i-1]+2==pinkdataodds$BY[i],pinkdataodds$S[i],0)
  OddT[i-1]<-ifelse(pinkdataodds$BY[i-1]+2==pinkdataodds$BY[i],pinkdataodds$t[i],0)
}
Ostock<-c()
ORegion<-c()
Oddr<-c()
for ( i in 1: 647){
  Oddr[i]<-ifelse(pinkdataodds$BY[i+1]==pinkdataodds$BY[i]+2,pinkdataodds$R[i],0)
  Ostock[i]<-ifelse(pinkdataodds$BY[i+1]==pinkdataodds$BY[i]+2,pinkdataodds$stock2[i],0)
  ORegion[i]<-ifelse(pinkdataodds$BY[i+1]==pinkdataodds$BY[i]+2,pinkdataodds$Region2[i],0)
}
Odds<-subset(Odds,Odds>0)
Oddr<-subset(Oddr,Oddr>0)
Ostock<-subset(Ostock,Ostock>0)
ORegion<-subset(ORegion,ORegion>0)
OddT<-subset(OddT,OddT>0)
pinkdataeven<-subset(pinkdata,pinkdata$BY%%2==0)
Evens<-c()
Estock<-c()
ERegion<-c()
Evenr<-c()
EvenT<-c()
for ( i in 2: 648){
  Evens[i-1]<-ifelse(pinkdataeven$BY[i-1]+2==pinkdataeven$BY[i],pinkdataeven$S[i],0)
  EvenT[i-1]<-ifelse(pinkdataeven$BY[i-1]+2==pinkdataeven$BY[i],pinkdataeven$t[i],0)
  }

for ( i in 1: 647){
  Evenr[i]<-ifelse(pinkdataeven$BY[i+1]==pinkdataeven$BY[i]+2,pinkdataeven$R[i],0)
  Estock[i]<-ifelse(pinkdataeven$BY[i+1]==pinkdataeven$BY[i]+2,pinkdataeven$stock2[i],0)
  ERegion[i]<-ifelse(pinkdataeven$BY[i+1]==pinkdataeven$BY[i]+2,pinkdataeven$Region2[i],0)
  
  }
Evens<-subset(Evens,Evens>0)
Evenr<-subset(Evenr,Evenr>0)
Estock<-subset(Estock,Estock>0)
ERegion<-subset(ERegion,ERegion>0)
EvenT<-subset(EvenT,EvenT>0)

# t is based on 1952(Syear) is 1
# t is based on 1953(Syear) 
OddT<-OddT-1
EvenT<-EvenT-1
###
S<-c(Evens,Odds)*1000
R<-c(Evenr,Oddr)*1000
stock<-c(Estock,Ostock)
region<-c(ERegion,ORegion)
#t<-c(EvenT,OddT2)
t<-c(ceiling(EvenT/2),OddT/2)




# Create a data list
dataList <- list(
  "R" =R,
  "S" =S,
  "stock"=stock,
  "t"=t)





# List of parameters to be monitored  
parameters <- c("Epink","Opink","gamma","tau","sigmaO","sigmaE","region1"
                 ,"region2","region3","region4","region5","region6","region7"
                 ,"region8","region9","region10","region11","region12","region13","region14","MUO","MUE") 


# JAGS Set-up
adaptSteps <- 5000              #number of steps to "tune" the samplers
burnInSteps <- 5000             #number of steps to "burn-in" the samplers
nChains <- 4                    #number of chains to run
numSavedSteps <- 12000          #total number of steps in chains to save
thinSteps <- 1                 	#number of steps to "thin" (1=keep every step)
nIter <- ceiling((numSavedSteps*thinSteps )/nChains) 	#steps per chain

# Create, initialize, and adapt the model
jagsModel1 <- jags.model("Poiss2.txt", 
                         data=dataList, 
                         n.chains=nChains, 
                         n.adapt=adaptSteps)

# Burn-in the algorithm
cat( "Burning in the MCMC chain...\n")
update(jagsModel1, n.iter=burnInSteps)

# Run MCMC algorithm
cat("Sampling final MCMC chain...\n" )
codaSamples <- coda.samples(jagsModel1, 
                            variable.names=parameters, 
                            n.iter=nIter, 
                            thin=thinSteps)

# Make trace plots and density plots


plot(codaSamples)

# Calculate numerical summaries for the posterior samples
summary(codaSamples)

dic1 <- dic.samples(jagsModel1, nIter)
dic1

gelman.diag(codaSamples)
