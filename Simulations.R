library(MASS)
library(nloptr)
library(numDeriv)
library(xtable)
source("Functions.R")

parN = list(beta=c(2.5,2.6,1.8,2),eta=c(2.8,1.9,1.5,1.2),sd=c(1.1,1.4,0.75),gamma=c(-1,0.6,2.3))    #45-50% censoring
parl = length(parN[[1]])
totparl = 2*parl
parlgamma = (parl-1)
namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho")

samsize= c(250)

# the numbers after SimulationCI indicate whether Z and W are binary/continuous
# 1 = continuous and 2 = binary, the first number indicates Z and the second number indicates W
# the order below also matches the design order in section 4 (simulation study) of the paper

for(l in samsize)
{
  nsim = 1000
  myseed = 876661
  message("sample size = ",l)
  SimulationCI11(l,nsim,myseed)
  SimulationCI12(l,nsim,myseed)
  SimulationCI21(l,nsim,myseed)
  SimulationCI22(l,nsim,myseed)
}
