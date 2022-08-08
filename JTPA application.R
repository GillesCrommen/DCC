fulldata<-read.csv("clean_dataset_JTPA.csv")

library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
source("Functions.R")

father = fulldata[ (fulldata$children == 1 &      # select a subset of fathers
                       fulldata$male == 1 ),]
dataset = subset(father, select = -c(1,5,8))

# Creating the histogram

par(mfrow=c(1,1))
nbreaks=round(sqrt(nrow(dataset)))
col=coloredhist(dataset$days,dataset$delta,nbreaks)
hist(dataset$days, breaks = nbreaks, freq = T, col = col , main= "Distribution of time until employment or censoring", xlab = "Days" )

# Data formatting and application results

Y = as.matrix(log(dataset$days))
Delta = as.matrix(dataset$delta)
dataset$intercept=rep(1,nrow(dataset))
X = as.matrix(subset(dataset, select = c(ncol(dataset),1,2,3,5)))
parl=ncol(X)+2
totparl=2*parl
parlgamma=parl-1
Z = as.matrix(dataset$jtpa)
W = as.matrix(dataset$treatment)
XandW = as.matrix(cbind(X,W))
n=nrow(dataset)
data=as.matrix(cbind(Y,Delta,X,Z,W))
namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,1}","beta_{T,1}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{T,1}","beta_{T,1}","beta_{T,1}","beta_{T,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho")

Application22(data)


