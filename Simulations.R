pack <- c("MASS", 
          "nloptr",
          "numDeriv",
          "xtable",   
          "foreach",
          "doParallel",
          "rlist",
          "copula",
          "VineCopula",
          "rvinecopulib",
          "EnvStats")

invisible(lapply(pack, library, character.only=T))

setwd("C:/Users/u0147039/Desktop/PhD/Papers/An instrumental variable approach under dependent censoring/R files/R scripts")

source("Functions.R")

setwd("C:/Users/u0147039/Desktop/PhD/Papers/An instrumental variable approach under dependent censoring/R files/Simulations")

namescoef =  c("beta_{T,0}","beta_{T,1}","alpha_T","lambda_T","beta_{C,0}","beta_{C,1}","alpha_C","lambda_C","sigma_T","sigma_C","rho")
totparl = length(namescoef)-3
parl = totparl/2
parlgamma = (parl-1)

nsim = 1000
myseed = 876661

samsize= c(250,500,1000)

# continuous Z

parN = list(beta=c(1.5,0.6,0.4,0.3),eta=c(1.6,0.4,-0.3,-0.2),sd=c(1.1,1.4,0.75),gamma=c(0.5,-0.4,1))    #55% censoring

for(l in samsize)
{
  message("sample size = ",l)
  SimulationCI11(l,nsim,myseed,parN,parl,totparl,parlgamma,namescoef)
  SimulationCI12(l,nsim,myseed,parN,parl,totparl,parlgamma,namescoef)
}

# binary Z

parN = list(beta=c(1.5,0.6,0.4,0.3),eta=c(1.6,0.4,-0.3,-0.2),sd=c(1.1,1.4,0.75),gamma=c(-1,0.6,2.3))   #55-60% censoring

for(l in samsize)
{
  message("sample size = ",l)
  SimulationCI21(l,nsim,myseed,parN,parl,totparl,parlgamma,namescoef)
  SimulationCI22(l,nsim,myseed,parN,parl,totparl,parlgamma,namescoef)
}

# misspecification ( 1 = Gumbel error, 2 = Frank Copula, 3 = probit control function, 4 = V is identically zero)

for(l in 1:4)
{
  message("misspecification = ",l)
  MisspecificationCI22(500,nsim,myseed,parN,parl,totparl,parlgamma,namescoef,l)
}



