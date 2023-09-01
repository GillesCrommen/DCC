setwd("C:/Users/u0147039/Desktop/PhD/Papers/An instrumental variable approach under dependent censoring/R files/JTPA data application/Data")

fulldata<-read.csv("clean_dataset_JTPA.csv")

setwd("C:/Users/u0147039/Desktop/PhD/Papers/An instrumental variable approach under dependent censoring/R files/R scripts")

library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
source("Functions.R")

father = fulldata[ (fulldata$children == 1 &
                       fulldata$male == 1 ),]
dataset = subset(father, select = -c(1,5,8))

# Summary of the data

summary(dataset)
summary(dataset[dataset$delta==1,])

dataset1=dataset[dataset$treatment==1,]
mean(dataset1$jtpa)
dataset0=dataset[dataset$treatment==0,]
mean(dataset0$jtpa)

summary(dataset1)
summary(dataset1[dataset1$delta==1,])
summary(dataset0)
summary(dataset0[dataset0$delta==1,])

head(dataset)
table(dataset$treatment,dataset$jtpa)

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
namescoef =  c("beta_{T,0}","beta_{T,1}","beta_{T,2}","beta_{T,3}","beta_{T,4}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","beta_{C,2}","beta_{C,3}","beta_{C,4}","alpha_C","lambda_C","sigma_T","sigma_C","rho")

Application22(data)


