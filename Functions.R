dat.sim.reg = function(n,par,iseed,Zbin,Wbin){
  
  set.seed(iseed)
  beta = par[[1]]
  eta = par[[2]]
  sd = par[[3]]
  gamma = par[[4]]
  
  mu = c(0,0)
  sigma = matrix(c(sd[1]^2,sd[1]*sd[2]*sd[3], sd[1]*sd[2]*sd[3], sd[2]^2),ncol=2)
  err = mvrnorm(n, mu =mu , Sigma=sigma)
  
  err1 = err[,1]
  err2 = err[,2]
  
  x0 = rep(1,n)  # to keep the intercept
  
  x1 = rnorm(n,0,1)
  
  if (Wbin==2)
  {
    W = sample(c(0,1), n, replace = TRUE)}
  else if (Wbin==1)
  {
    W = runif(n,0,2)
  }
  
  XandW=as.matrix(cbind(x0,x1,W))
  
  if (Zbin==2)
  {
    V=rlogis(n)
    Z = as.matrix(as.numeric(XandW%*%gamma-V>0))
    realV=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  }
  else if (Zbin==1)
  {
    V=rnorm(n,0,2)
    Z = XandW%*%gamma+V
    realV= Z-(XandW%*%gamma)
  }
  
  Mgen = matrix(c(x0,x1,Z,realV),ncol=parl,nrow=n) 
  T = Mgen%*%beta+err1
  C = Mgen%*%eta+err2
  
  M = matrix(c(x0,x1,Z,W),ncol=parl,nrow=n)    # data matrix
  
  Y = pmin(T,C)
  d1 = as.numeric(Y==T)
  data = cbind(Y,d1,M,realV)
  
  return(data)
}

LikF = function(par,Y,Delta,M){ # joint model with dependent censoring
  M=as.matrix(M)
  k = ncol(M)
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k])
  eta = as.matrix(par[v:l])
  sigma1 = par[l+1]
  sigma2 = par[l+2]
  rho = par[l+3]
  
  z1 = (Y-(M%*%beta))/sigma1
  z2 = ((1-(rho*sigma2/sigma1))*Y-(M%*%eta-rho*(sigma2/sigma1)*(M%*%beta)))/(sigma2*((1-rho^2)^0.5))
  z3 = (Y-(M%*%eta))/sigma2
  z4 = ((1-(rho*sigma1/sigma2))*Y-(M%*%beta-rho*(sigma1/sigma2)*(M%*%eta)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}

LikFG1 = function(par,Y,Delta,M){ # Needed for Hessian gamma when Z is continuous
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  rho = par[l+5]
  gamma = as.matrix(par[(l+6):(l+5+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=Z-XandW%*%gamma
  
  z1 = (Y-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((1-rho*sigma2/sigma1)*Y-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (Y-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((1-rho*sigma1/sigma2)*Y-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}


LikFG2 = function(par,Y,Delta,M){ # needed for Hessian gamma when Z is binary
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  rho = par[l+5]
  gamma = as.matrix(par[(l+6):(l+5+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  
  z1 = (Y-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = ((1-rho*sigma2/sigma1)*Y-((X%*%eta+Z*alphaC+Vest*lambdaC)-rho*(sigma2/sigma1)*(X%*%beta+Z*alphaT+Vest*lambdaT)))/(sigma2*(1-rho^2)^0.5)
  z3 = (Y-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  z4 = ((1-rho*sigma1/sigma2)*Y-((X%*%beta+Z*alphaT+Vest*lambdaT)-rho*(sigma1/sigma2)*(X%*%eta+Z*alphaC+Vest*lambdaC)))/(sigma1*(1-rho^2)^0.5)
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z3)*(1-pnorm(z4)))^(1-Delta)
  p1 = pmax(tot,1e-100)   
  Logn = sum(log(p1)); 
  return(-Logn)
}

LikI = function(par,Y,Delta,M){ # independence model assumption (rho = 0)
  M=as.matrix(M)
  k = ncol(M)
  l = 2*k
  v = k+1
  beta = as.matrix(par[1:k])
  eta = as.matrix(par[v:l])
  sigma1 = par[l+1]
  sigma2 = par[l+2]
  
  z1 = (Y-(M%*%beta))/sigma1
  z2 = (Y-(M%*%eta))/sigma2
  
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta))
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}

LikGamma1 = function(par,Y,M){ # maximum likelihood for Gamma
  M=as.matrix(M)
  gamma= as.matrix(par)
  
  tot = (Y-M%*%gamma)^2
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(Logn)
}

LikGamma2 = function(par,Y,M){ # maximum likelihood for Gamma
  M=as.matrix(M)
  gamma= as.matrix(par)
  
  tot = (plogis(M%*%gamma)^Y)*((1-plogis(M%*%gamma))^(1-Y))
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}

LikIGamma1 = function(par,Y,Delta,M){ # no dependent censoring
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  gamma = as.matrix(par[(l+5):(l+4+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=Z-(XandW%*%gamma)
  
  z1 = (Y-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = (Y-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*(((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta))
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}

LikIGamma2 = function(par,Y,Delta,M){ # no dependent censoring
  M=as.matrix(M)
  k = ncol(M)-2
  l = 2*(k+1)
  v = k+3
  beta = as.matrix(par[1:k])
  alphaT = par[k+1]
  lambdaT = par[k+2]
  eta = as.matrix(par[v:l])
  alphaC = par[l+1]
  lambdaC = par[l+2]
  sigma1 = par[l+3]
  sigma2 = par[l+4]
  gamma = as.matrix(par[(l+5):(l+4+parlgamma)])
  
  X=as.matrix(M[,1:k])
  Z=as.matrix(M[,k+1])
  W=as.matrix(M[,k+2])
  XandW=as.matrix(cbind(X,W))
  Vest=(1-Z)*((1+exp(XandW%*%gamma))*log(1+exp(XandW%*%gamma))-(XandW%*%gamma)*exp(XandW%*%gamma))-Z*((1+exp(-(XandW%*%gamma)))*log(1+exp(-(XandW%*%gamma)))+(XandW%*%gamma)*exp(-(XandW%*%gamma)))
  
  z1 = (Y-(X%*%beta+Z*alphaT+Vest*lambdaT))/sigma1
  z2 = (Y-(X%*%eta+Z*alphaC+Vest*lambdaC))/sigma2
  
  tot = (((1/sigma1)*dnorm(z1)*(1-pnorm(z2)))^Delta)*((1/sigma2)*dnorm(z2)*(1-pnorm(z1)))^(1-Delta)
  p1 = pmax(tot,1e-100)
  Logn = sum(log(p1)); 
  return(-Logn)
}

SimulationCI11 = function(n,nsim,iseed)
{
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim)
  {
    
    if ( round(i %% (nsim/10)) == 0)
    {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,1,1)
    
    
    Y = data[,1]
    Delta = data[,2]
    X = data[,(4:parl)]
    Z = data[,parl+1]
    W = data[,parl+2]
    XandW = cbind(data[,3],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    M = cbind(data[,3:(1+parl)],V)
    MnoV = data[,3:(2+parl)]
    MrealV = cbind(data[,3:(1+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[1]
    
    init = c(rep(0,totparl),1,1) # Starting values
    
    # Independent model for starting values sigma
    
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5),ub=c(rep(Inf,totparl),Inf,Inf),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V
    
    ME = M[,-ncol(M)]
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    initE = c(initE,0)
    
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-1),ub=c(rep(Inf,(totparl-2)),Inf,Inf,1),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    
    t_s1 = 1/parhatE[totparl-1]*se1[totparl]
    t_s2 = 1/parhatE[totparl-1]*se1[totparl]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
    ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1])))     # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l),ncol=1), matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u),ncol=1)) 
    
    # Model with estimated V
    
    initd = c(parhat1,0)
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    WMI = ginv(WM)
    
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    psii = -WMI%*%mi
    
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l),ncol=1),matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u), ncol=1))
    
    # Model with real V
    
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l),ncol=1),matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma1,parhatGI,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = Hgamma[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      giI = rbind(giI,c(J1I))
    }
    
    giI = t(giI)
    
    partvarI = giI + VargammaI%*%psii
    
    Epartvar2I = (partvarI%*%t(partvarI))
    
    totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
    
    seI = sqrt(abs(diag(totvarexI)))
    
    # Delta method variance
    
    se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
    se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
    st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
    
    # Back transform
    
    s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI) 
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI),ncol=1),matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
    
  }
  
  print(per/(n*nsim))     #percentage of censoring
  
  ## Results of model with estimated V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+3)],2,sd)
  MSD  = apply(results[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with no V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+1),byrow=TRUE)
  
  Bias = apply(results1[,1:(totparl+1)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+1)],2,sd)
  MSD  = apply(results1[,(totparl+2):(2*totparl+2)],2, mean)
  RMSE = sqrt(apply((results1[,1:(totparl+1)]-par0m)^2,2,mean))
  
  
  CP = rep(0,(totparl+1))
  datacp = results1[,(2*totparl+3):(4*totparl+4)]
  for(i in 1:(totparl+1)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with real V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results2[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+3)],2,sd)
  MSD  = apply(results2[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results2[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results2[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V but independence
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2])
  par0m = matrix(par0,nsim,(totparl+2),byrow=TRUE)
  
  Bias = apply(results3[,1:(totparl+2)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+2)],2,sd)
  MSD  = apply(results3[,(totparl+3):(2*totparl+4)],2, mean)
  RMSE = sqrt(apply((results3[,1:(totparl+2)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+2)
  datacp = results3[,(2*totparl+5):(4*totparl+8)]
  for(i in 1:(totparl+2))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  xtab = xtable(sum)
  digits(xtab) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab,file=paste0("estV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("noV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("realV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = namescoef[1:length(namescoef)-1]
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("IndEstV11_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}

SimulationCI12 = function(n,nsim,iseed)
{
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim)
  {
    
    if ( round(i %% (nsim/10)) == 0)
    {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,1,2)
    
    Y = data[,1]
    Delta = data[,2]
    X = data[,(4:parl)]
    Z = data[,parl+1]
    W = data[,parl+2]
    XandW = cbind(data[,3],X,W)
    
    gammaest <- lm(Z~X+W)$coefficients
    V <- Z-(XandW%*%gammaest)
    
    M = cbind(data[,3:(1+parl)],V)
    MnoV = data[,3:(2+parl)]
    MrealV = cbind(data[,3:(1+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[1]
    
    init = c(rep(0,totparl),1,1) # Starting values
    
    # Independent model for starting values sigma
    
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5),ub=c(rep(Inf,totparl),Inf,Inf),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V
    
    ME = M[,-ncol(M)]
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    initE = c(initE,0)
    
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-1),ub=c(rep(Inf,(totparl-2)),Inf,Inf,1),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    
    t_s1 = 1/parhatE[totparl-1]*se1[totparl]
    t_s2 = 1/parhatE[totparl-1]*se1[totparl]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
    ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1])))     # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l),ncol=1), matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u),ncol=1)) 
    
    # Model with estimated V
    
    initd = c(parhat1,0)
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    
    Hgamma = hessian(LikFG1,parhatG,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    sumsecder = c(rep(0,ncol(prodvec)))
    
    for (i in 1:length(sumsecder)) {
      sumsecder[i]= -sum(prodvec[,i])
    }
    
    WM = sumsecder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-sumsecder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    WMI = ginv(WM)
    
    mi = c()
    
    for(i in 1:n){
      newrow<-V[i]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    mi=t(mi)
    
    psii = -WMI%*%mi
    
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l),ncol=1),matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u), ncol=1))
    
    # Model with real V
    
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l),ncol=1),matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma1,parhatGI,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = Hgamma[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      giI = rbind(giI,c(J1I))
    }
    
    giI = t(giI)
    
    partvarI = giI + VargammaI%*%psii
    
    Epartvar2I = (partvarI%*%t(partvarI))
    
    totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
    
    seI = sqrt(abs(diag(totvarexI)))
    
    # Delta method variance
    
    se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
    se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
    st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
    
    # Back transform
    
    s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI) 
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI),ncol=1),matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
    
  }
  
  print(per/(n*nsim))     #percentage of censoring
  
  ## Results of model with estimated V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+3)],2,sd)
  MSD  = apply(results[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with no V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+1),byrow=TRUE)
  
  Bias = apply(results1[,1:(totparl+1)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+1)],2,sd)
  MSD  = apply(results1[,(totparl+2):(2*totparl+2)],2, mean)
  RMSE = sqrt(apply((results1[,1:(totparl+1)]-par0m)^2,2,mean))
  
  
  CP = rep(0,(totparl+1))
  datacp = results1[,(2*totparl+3):(4*totparl+4)]
  for(i in 1:(totparl+1)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with real V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results2[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+3)],2,sd)
  MSD  = apply(results2[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results2[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results2[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V but independence
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2])
  par0m = matrix(par0,nsim,(totparl+2),byrow=TRUE)
  
  Bias = apply(results3[,1:(totparl+2)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+2)],2,sd)
  MSD  = apply(results3[,(totparl+3):(2*totparl+4)],2, mean)
  RMSE = sqrt(apply((results3[,1:(totparl+2)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+2)
  datacp = results3[,(2*totparl+5):(4*totparl+8)]
  for(i in 1:(totparl+2))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  xtab = xtable(sum)
  digits(xtab) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab,file=paste0("estV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("noV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("realV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = namescoef[1:length(namescoef)-1]
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("IndEstV12_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}

SimulationCI21 = function(n,nsim,iseed)
{
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim)
  {
    
    if ( round(i %% (nsim/10)) == 0)
    {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,2,1)
    
    Y = data[,1]
    Delta = data[,2]
    X = data[,(4:parl)]
    Z = data[,parl+1]
    W = data[,parl+2]
    XandW = cbind(data[,3],X,W)
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    M = cbind(data[,3:(1+parl)],V)
    MnoV = data[,3:(2+parl)]
    MrealV = cbind(data[,3:(1+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[1]
    
    init = c(rep(0,totparl),1,1) # Starting values
    
    # Independent model for starting values
    
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5),ub=c(rep(Inf,totparl),Inf,Inf),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V
    
    ME = M[,-ncol(M)]
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    initE = c(initE,0)
    
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-1),ub=c(rep(Inf,(totparl-2)),Inf,Inf,1),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    
    t_s1 = 1/parhatE[totparl-1]*se1[totparl]
    t_s2 = 1/parhatE[totparl-1]*se1[totparl]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
    ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1]))) # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l),ncol=1), matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u),ncol=1)) 
    
    # Model with estimated V
    
    initd = c(parhat1,0)
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    secder=t(-dlogis(XandW%*%gammaest))%*%prodvec
    
    WM = secder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-secder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    WMI = ginv(WM)
    
    diffvec = Z-plogis(XandW%*%gammaest)
    
    mi = c()
    
    for(i in 1:n){
      newrow<-diffvec[i,]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    psii = -WMI%*%t(mi)
    
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l),ncol=1),matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u), ncol=1))
    
    # Model with real V
    
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l),ncol=1),matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma2,parhatGI,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = Hgamma[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      giI = rbind(giI,c(J1I))
    }
    
    giI = t(giI)
    
    partvarI = giI + VargammaI%*%psii
    
    Epartvar2I = (partvarI%*%t(partvarI))
    
    totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
    
    seI = sqrt(abs(diag(totvarexI)))
    
    # Delta method variance
    
    se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
    se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
    st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
    
    # Back transform
    
    s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI) 
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI),ncol=1),matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
    
  }
  
  print(per/(n*nsim))     #percentage of censoring
  
  ## Results of model with estimated V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+3)],2,sd)
  MSD  = apply(results[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with no V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+1),byrow=TRUE)
  
  Bias = apply(results1[,1:(totparl+1)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+1)],2,sd)
  MSD  = apply(results1[,(totparl+2):(2*totparl+2)],2, mean)
  RMSE = sqrt(apply((results1[,1:(totparl+1)]-par0m)^2,2,mean))
  
  
  CP = rep(0,(totparl+1))
  datacp = results1[,(2*totparl+3):(4*totparl+4)]
  for(i in 1:(totparl+1)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with real V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results2[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+3)],2,sd)
  MSD  = apply(results2[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results2[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results2[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V but independence
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2])
  par0m = matrix(par0,nsim,(totparl+2),byrow=TRUE)
  
  Bias = apply(results3[,1:(totparl+2)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+2)],2,sd)
  MSD  = apply(results3[,(totparl+3):(2*totparl+4)],2, mean)
  RMSE = sqrt(apply((results3[,1:(totparl+2)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+2)
  datacp = results3[,(2*totparl+5):(4*totparl+8)]
  for(i in 1:(totparl+2))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  xtab = xtable(sum)
  digits(xtab) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab,file=paste0("estV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("noV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("realV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = namescoef[1:length(namescoef)-1]
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("IndEstV21_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}


SimulationCI22 = function(n,nsim,iseed)
{
  sum = c()
  sum1 = c()
  sum2 = c()
  sum3 = c()
  per=0
  results = c()
  results1 = c()
  results2 = c()
  results3 = c()
  
  for (i in 1:nsim)
  {
    
    if ( round(i %% (nsim/10)) == 0)
    {cat((i/nsim)*100,"%", "\n", sep="")}
    
    data = dat.sim.reg(n,parN,iseed+i,2,2)
    
    Y = data[,1]
    Delta = data[,2]
    X = data[,(4:parl)]
    Z = data[,parl+1]
    W = data[,parl+2]
    XandW = cbind(data[,3],X,W)
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    M = cbind(data[,3:(1+parl)],V)
    MnoV = data[,3:(2+parl)]
    MrealV = cbind(data[,3:(1+parl)],data[,ncol(data)])
    
    per=per+table(Delta)[1]
    
    init = c(rep(0,totparl),1,1) # Starting values
    
    # Independent model for starting values
    
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5),ub=c(rep(Inf,totparl),Inf,Inf),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V
    
    ME = M[,-ncol(M)]
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    initE = c(initE,0)
    
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-1),ub=c(rep(Inf,(totparl-2)),Inf,Inf,1),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    
    t_s1 = 1/parhatE[totparl-1]*se1[totparl]
    t_s2 = 1/parhatE[totparl-1]*se1[totparl]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
    ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1]))) # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l),ncol=1), matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u),ncol=1)) 
    
    # Model with estimated V
    
    initd = c(parhat1,0)
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    secder=t(-dlogis(XandW%*%gammaest))%*%prodvec
    
    WM = secder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-secder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    WMI = ginv(WM)
    
    diffvec = Z-plogis(XandW%*%gammaest)
    
    mi = c()
    
    for(i in 1:n){
      newrow<-diffvec[i,]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    psii = -WMI%*%t(mi)
    
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l),ncol=1),matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u), ncol=1))
    
    # Model with real V
    
    parhatre = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=MrealV,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                      eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    Hre = hessian(LikF,parhatre,Y=Y,Delta=Delta,M=MrealV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    HreI = ginv(Hre)
    
    sere = sqrt(abs(diag(HreI)))
    
    # Delta method variance
    
    sere_s1 = 1/parhatre[totparl+1]*sere[totparl+1]
    sere_s2 = 1/parhatre[totparl+2]*sere[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1re_l = log(parhatre[totparl+1])-1.96*sere_s1 ;  st1re_u = log(parhatre[totparl+1])+1.96*sere_s1 
    st2re_l = log(parhatre[totparl+2])-1.96*sere_s2 ;  st2re_u = log(parhatre[totparl+2])+1.96*sere_s2 
    
    # Back transfrom
    
    s1re_l = exp(st1re_l); s1re_u = exp(st1re_u); s2re_l = exp(st2re_l); s2re_u = exp(st2re_u) 
    
    # Confidence interval for rho
    
    ztre = 0.5*(log((1+parhatre[totparl+3])/(1-parhatre[totparl+3])))     # Fisher's z transform
    sere_z = (1/(1-parhatre[totparl+3]^2))*sere[totparl+3]
    ztre_l = ztre-1.96*(sere_z)
    ztre_u = ztre+1.96*(sere_z)
    
    # Back transform
    
    rre_l = (exp(2*ztre_l)-1)/(exp(2*ztre_l)+1)      
    rre_u = (exp(2*ztre_u)-1)/(exp(2*ztre_u)+1)
    
    EC3 = cbind(matrix(c(parhatre[1:totparl]-1.96*(sere[1:totparl]),s1re_l,s2re_l,rre_l),ncol=1),matrix(c(parhatre[1:totparl]+1.96*(sere[1:totparl]),s1re_u,s2re_u,rre_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma2,parhatGI,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = Hgamma[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      giI = rbind(giI,c(J1I))
    }
    
    giI = t(giI)
    
    partvarI = giI + VargammaI%*%psii
    
    Epartvar2I = (partvarI%*%t(partvarI))
    
    totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
    
    seI = sqrt(abs(diag(totvarexI)))
    
    # Delta method variance
    
    se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
    se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
    st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
    
    # Back transform
    
    s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI) 
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI),ncol=1),matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results2 = rbind(results2,c(parhatre,sere,c(t(EC3))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
    
  }
  
  print(per/(n*nsim))     #percentage of censoring
  
  ## Results of model with estimated V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results[,1:(totparl+3)],2,sd)
  MSD  = apply(results[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with no V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0 = par0[-(parl)]
  par0 = par0[-(2*parl-1)]
  par0m = matrix(par0,nsim,(totparl+1),byrow=TRUE)
  
  Bias = apply(results1[,1:(totparl+1)]-par0m,2,mean)
  ESE = apply(results1[,1:(totparl+1)],2,sd)
  MSD  = apply(results1[,(totparl+2):(2*totparl+2)],2, mean)
  RMSE = sqrt(apply((results1[,1:(totparl+1)]-par0m)^2,2,mean))
  
  
  CP = rep(0,(totparl+1))
  datacp = results1[,(2*totparl+3):(4*totparl+4)]
  for(i in 1:(totparl+1)){
    index = c(2*i-1,2*i)
    CP[i] = sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary1 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Model with real V
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]])
  par0m = matrix(par0,nsim,(totparl+3),byrow=TRUE)
  
  Bias = apply(results2[,1:(totparl+3)]-par0m,2,mean)
  ESE = apply(results2[,1:(totparl+3)],2,sd)
  MSD  = apply(results2[,(totparl+4):(2*totparl+6)],2, mean)
  RMSE = sqrt(apply((results2[,1:(totparl+3)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+3)
  datacp = results2[,(2*totparl+7):(4*totparl+12)]
  for(i in 1:(totparl+3))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary2 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  ## Results of model with estimated V but independence
  
  par0 = c(parN[[1]],parN[[2]],parN[[3]][1],parN[[3]][2])
  par0m = matrix(par0,nsim,(totparl+2),byrow=TRUE)
  
  Bias = apply(results3[,1:(totparl+2)]-par0m,2,mean)
  ESE = apply(results3[,1:(totparl+2)],2,sd)
  MSD  = apply(results3[,(totparl+3):(2*totparl+4)],2, mean)
  RMSE = sqrt(apply((results3[,1:(totparl+2)]-par0m)^2,2,mean))
  
  CP = rep(0,totparl+2)
  datacp = results3[,(2*totparl+5):(4*totparl+8)]
  for(i in 1:(totparl+2))
  {
    index=c(2*i-1,2*i)
    CP[i]=sum(datacp[,index[1]]<=par0[i] & datacp[,index[2]]>=par0[i])/nsim
  } 
  
  summary3 = cbind(Bias,ESE,MSD,RMSE,CP) 
  
  sum = summary
  sum1 = summary1
  sum2 = summary2
  sum3 = summary3
  
  ## Results of model with estimated V
  
  colnames(sum) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum) = namescoef
  xtab = xtable(sum)
  digits(xtab) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab,file=paste0("estV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with no V
  
  colnames(sum1)=c("Bias","ESD","ASE","RMSE","CR")
  namescoefr=namescoef[-(parl)]
  namescoefr=namescoefr[-(2*parl-1)]
  rownames(sum1)=namescoefr
  xtab1 = xtable(sum1)
  digits(xtab1) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab1,file=paste0("noV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with real V
  
  colnames(sum2) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum2) = namescoef
  xtab2 = xtable(sum2)
  digits(xtab2) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab2,file=paste0("realV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab2, add.to.row=addtorow, include.colnames=TRUE)
  
  ## Results of model with estimated V but independence
  
  colnames(sum3) = c("Bias","ESD","ASE","RMSE","CR")
  rownames(sum3) = namescoef[1:length(namescoef)-1]
  xtab3 = xtable(sum3)
  digits(xtab3) = rep(3,6)
  header= c("sample size",n)
  addtorow = list()
  addtorow$pos = list(-1)
  addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
  
  print.xtable(xtab3,file=paste0("IndEstV22_",n,".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
  print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}

Application22 = function(data)
{
  results = c()
  results1 = c()
  results3 = c()
  
    n = nrow(data)
    Y = data[,1]
    Delta = data[,2]
    X = data[,(4:parl)]
    Z = data[,parl+1]
    W = data[,parl+2]
    XandW = as.matrix(cbind(data[,3],X,W))
    
    gammaest <- nloptr(x0=rep(0,parlgamma),eval_f=LikGamma2,Y=Z,M=XandW,lb=c(rep(-Inf,parlgamma)),ub=c(rep(Inf,parlgamma)),
                       eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    V <- (1-Z)*((1+exp(XandW%*%gammaest))*log(1+exp(XandW%*%gammaest))-(XandW%*%gammaest)*exp(XandW%*%gammaest))-Z*((1+exp(-(XandW%*%gammaest)))*log(1+exp(-(XandW%*%gammaest)))+(XandW%*%gammaest)*exp(-(XandW%*%gammaest)))
    
    M = cbind(data[,3:(1+parl)],V)
    MnoV = data[,3:(2+parl)]
    
    init = c(rep(0,totparl),1,1) # Starting values
    
    # Independent model for starting values
    
    parhat1 = nloptr(x0=c(init),eval_f=LikI,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5),ub=c(rep(Inf,totparl),Inf,Inf),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    # Model with no V
    
    ME = M[,-ncol(M)]
    initE = parhat1[-parl]
    initE = initE[-(2*parl-1)]
    initE = c(initE,0)
    
    parhatE = nloptr(x0=initE,eval_f=LikF,Y=Y,Delta=Delta,M=ME,lb=c(rep(-Inf,(totparl-2)),1e-05,1e-5,-1),ub=c(rep(Inf,(totparl-2)),Inf,Inf,1),
                     eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    H1 = hessian(LikF,parhatE,Y=Y,Delta=Delta,M=ME,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    H1I = ginv(H1)
    se1 = sqrt(abs(diag(H1I)));
    
    # Delta method variance (makes sure no negative values in CI for variance)
    
    t_s1 = 1/parhatE[totparl-1]*se1[totparl]
    t_s2 = 1/parhatE[totparl-1]*se1[totparl]
    
    # Conf. interval for transf. sigma's
    
    ms1_l = log(parhatE[totparl-1])-1.96*t_s1 ;  ms1_u = log(parhatE[totparl-1])+1.96*t_s1 
    ms2_l = log(parhatE[totparl])-1.96*t_s2 ;  ms2_u = log(parhatE[totparl])+1.96*t_s2 
    
    # Back transform
    
    S1_l = exp(ms1_l); S1_u = exp(ms1_u); S2_l = exp(ms2_l); S2_u = exp(ms2_u) 
    
    # Confidence interval for rho
    
    z1t = 0.5*(log((1+parhatE[totparl+1])/(1-parhatE[totparl+1]))) # Fisher's z transform
    se1_z = (1/(1-parhatE[totparl+1]^2))*se1[totparl+1]
    z1t_l = z1t-1.96*(se1_z)
    z1t_u = z1t+1.96*(se1_z)
    
    # Back transform
    
    r1_l = (exp(2*z1t_l)-1)/(exp(2*z1t_l)+1)      
    r1_u = (exp(2*z1t_u)-1)/(exp(2*z1t_u)+1)
    
    EC2 = cbind(matrix(c(parhatE[1:(totparl-2)]-1.96*(se1)[1:(totparl-2)],S1_l,S2_l,r1_l),ncol=1), matrix(c(parhatE[1:(totparl-2)]+1.96*(se1)[1:(totparl-2)],S1_u,S2_u,r1_u),ncol=1)) 
    
    # Model with estimated V
    
    initd = c(parhat1,0)
    
    parhat = nloptr(x0=initd,eval_f=LikF,Y=Y,Delta=Delta,M=M,lb=c(rep(-Inf,totparl),1e-05,1e-5,-1),ub=c(rep(Inf,totparl),Inf,Inf,1),
                    eval_g_ineq=NULL,opts = list(algorithm = "NLOPT_LN_BOBYQA","ftol_abs"=1.0e-30,"maxeval"=100000,"xtol_abs"=rep(1.0e-30)))$solution
    
    parhatG = c(parhat,as.vector(gammaest))
    
    Hgamma = hessian(LikFG2,parhatG,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    H = Hgamma[1:length(initd),1:length(initd)]
    HI = ginv(H)
    
    Vargamma = Hgamma[1:length(initd),(length(initd)+1):(length(initd)+parlgamma)]
    
    prodvec = XandW[,1]
    
    for (i in 1:parlgamma) {
      for (j in 2:parlgamma) {
        if (i<=j){
          prodvec<-cbind(prodvec,diag(XandW[,i]%*%t(XandW[,j])))
        }
      }
    }
    
    secder=t(-dlogis(XandW%*%gammaest))%*%prodvec
    
    WM = secder[1:parlgamma]
    for (i in 2:parlgamma) {
      newrow<-secder[c(i,(i+2):(i+parlgamma))]
      WM<-rbind(WM,newrow) 
    }
    
    WMI = ginv(WM)
    
    diffvec = Z-plogis(XandW%*%gammaest)
    
    mi = c()
    
    for(i in 1:n){
      newrow<-diffvec[i,]%*%XandW[i,]
      mi = rbind(mi,newrow)
    }
    
    psii = -WMI%*%t(mi)
    
    gi = c()
    
    for (i in 1:n)
    {
      J1 = jacobian(LikF,parhat,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson")
      gi = rbind(gi,c(J1))
    }
    
    gi = t(gi)
    
    partvar = gi + Vargamma%*%psii
    
    Epartvar2 = (partvar%*%t(partvar))
    
    totvarex = HI%*%Epartvar2%*%t(HI)
    
    se = sqrt(abs(diag(totvarex)))
    
    # Delta method variance
    
    se_s1 = 1/parhat[totparl+1]*se[totparl+1]
    se_s2 = 1/parhat[totparl+2]*se[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_l = log(parhat[totparl+1])-1.96*se_s1 ;  st1_u = log(parhat[totparl+1])+1.96*se_s1  
    st2_l = log(parhat[totparl+2])-1.96*se_s2 ;  st2_u = log(parhat[totparl+2])+1.96*se_s2 
    
    # Back transform
    
    s1_l = exp(st1_l); s1_u = exp(st1_u); s2_l = exp(st2_l); s2_u = exp(st2_u) 
    
    # Confidence interval for rho
    
    zt = 0.5*(log((1+parhat[totparl+3])/(1-parhat[totparl+3])))     # Fisher's z transform
    se_z = (1/(1-parhat[totparl+3]^2))*se[totparl+3]
    zt_l = zt-1.96*(se_z)
    zt_u = zt+1.96*(se_z)
    
    # Back transform
    
    r_l = (exp(2*zt_l)-1)/(exp(2*zt_l)+1)      
    r_u = (exp(2*zt_u)-1)/(exp(2*zt_u)+1)
    
    EC1 = cbind(matrix(c(parhat[1:totparl]-1.96*(se[1:totparl]),s1_l,s2_l,r_l),ncol=1),matrix(c(parhat[1:totparl]+1.96*(se[1:totparl]),s1_u,s2_u,r_u), ncol=1))
    
    # Model with estimated V but assuming independence
    
    parhatGI = c(parhat1,as.vector(gammaest))
    
    HgammaI = hessian(LikIGamma2,parhatGI,Y=Y,Delta=Delta,M=MnoV,method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)) 
    
    HInd = HgammaI[1:(length(initd)-1),1:(length(initd)-1)]
    HIInd = ginv(HInd)
    
    VargammaI = Hgamma[1:(length(initd)-1),(length(initd)):(length(initd)+parlgamma-1)]
    
    giI = c()
    
    for (i in 1:n)
    {
      J1I = jacobian(LikI,parhat1,Y=Y[i],Delta=Delta[i],M=t(M[i,]),method="Richardson",method.args=list(eps=1e-4, d=0.0001, zer.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
      giI = rbind(giI,c(J1I))
    }
    
    giI = t(giI)
    
    partvarI = giI + VargammaI%*%psii
    
    Epartvar2I = (partvarI%*%t(partvarI))
    
    totvarexI = HIInd%*%Epartvar2I%*%t(HIInd)
    
    seI = sqrt(abs(diag(totvarexI)))
    
    # Delta method variance
    
    se_s1I = 1/parhat1[totparl+1]*seI[totparl+1]
    se_s2I = 1/parhat1[totparl+2]*seI[totparl+2]
    
    # Conf. interval for transf. sigma's
    
    st1_lI = log(parhat1[totparl+1])-1.96*se_s1I ;  st1_uI = log(parhat1[totparl+1])+1.96*se_s1I  
    st2_lI = log(parhat1[totparl+2])-1.96*se_s2I ;  st2_uI = log(parhat1[totparl+2])+1.96*se_s2I 
    
    # Back transform
    
    s1_lI = exp(st1_lI); s1_uI = exp(st1_uI); s2_lI = exp(st2_lI); s2_uI = exp(st2_uI) 
    
    EC4 = cbind(matrix(c(parhat1[1:totparl]-1.96*(seI[1:totparl]),s1_lI,s2_lI),ncol=1),matrix(c(parhat1[1:totparl]+1.96*(seI[1:totparl]),s1_uI,s2_uI), ncol=1))
    
    results = rbind(results,c(parhat,se,c(t(EC1))))
    results1 = rbind(results1,c(parhatE,se1,c(t(EC2))))
    results3 = rbind(results3,c(parhat1,seI,c(t(EC4))))
    
    ## Results of model with estimated V
    
    est = results[1:(totparl+3)]
    SE = results[(totparl+4):(2*totparl+6)]
    CIL = EC1[,1]
    CIU = EC1[,2]
    pvalue = 2*pmin((1-pnorm(est/SE)),pnorm(est/SE))
    
    summary = cbind(est,SE,CIL,CIU,pvalue) 
    
    ## Model with no V
    
    est = results1[1:(totparl+1)]
    SE = results1[(totparl+2):(2*totparl+2)]
    CIL = EC2[,1]
    CIU = EC2[,2]
    pvalue = 2*pmin((1-pnorm(est/SE)),pnorm(est/SE))
    
    summary1 = cbind(est,SE,CIL,CIU,pvalue) 
    
    ## Results of model with estimated V but independence
    
    est = results3[1:(totparl+2)]
    SE = results3[(totparl+3):(2*totparl+4)]
    CIL = EC4[,1]
    CIU = EC4[,2]
    pvalue = 2*pmin((1-pnorm(est/SE)),pnorm(est/SE))
    
    summary3 = cbind(est,SE,CIL,CIU,pvalue) 
    
    ## Results of model with estimated V
    
    colnames(summary) = c("Estimate","SD","lower","upper","p-value")
    rownames(summary) = namescoef
    xtab = xtable(summary)
    digits(xtab) = rep(3,6)
    header= c("sample size",n,"Results_2-step_Estimation")
    addtorow = list()
    addtorow$pos = list(-1)
    addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
    print.xtable(xtab,file=paste0("Results_2-step_Estimation",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
    print(xtab, add.to.row=addtorow, include.colnames=TRUE)
    
    ## Results of model with no V
    
    colnames(summary1)=c("Estimate","SD","lower","upper","p-value")
    namescoefr=namescoef[-(parl)]
    namescoefr=namescoefr[-(2*parl-1)]
    rownames(summary1)=namescoefr
    xtab1 = xtable(summary1)
    digits(xtab1) = rep(3,6)
    header= c("sample size",n,"Results_No_Confounding")
    addtorow = list()
    addtorow$pos = list(-1)
    addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
    
    print.xtable(xtab1,file=paste0("Results_No_Confounding",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
    print(xtab1, add.to.row=addtorow, include.colnames=TRUE)
    
    ## Results of model with estimated V but independence
    
    colnames(summary3) = c("Estimate","SD","lower","upper","p-value")
    rownames(summary3) = namescoef[1:length(namescoef)-1]
    xtab3 = xtable(summary3)
    digits(xtab3) = rep(3,6)
    header= c("sample size",n,"Results_No_Dependence")
    addtorow = list()
    addtorow$pos = list(-1)
    addtorow$command = paste0(paste0('& \\multicolumn{1}{c}{', header, '}', collapse=''), '\\\\')
    
    print.xtable(xtab3,file=paste0("Results_No_Dependence",".txt"),add.to.row=addtorow,append=TRUE,table.placement="!")
    print(xtab3, add.to.row=addtorow, include.colnames=TRUE)
}

coloredhist=function(days,delta,nbreaks){
  
  n = nbreaks
  
  histbreaks<-hist(days, breaks = n, plot = F)$breaks
  
  sngmhist = as.data.frame(cbind(delta, (round(days/histbreaks[2])+1)))
  sngmhist <- sngmhist[order(sngmhist[,2]),]
  total<-as.data.frame(table(sngmhist[,2]))
  sngmhist0 = sngmhist[sngmhist[,1]==0,]
  
  censoringmatrix = matrix(0,nrow=nrow(sngmhist0),ncol = length(histbreaks))
  
  for (i in 1:nrow(sngmhist0)) {
    censoringmatrix[i,(sngmhist0[i,2])]= 1
  }
  
  censoringbin = c(rep(0,ncol(censoringmatrix)))
  
  for (i in 1:length(censoringbin)) {
    censoringbin[i]=sum(censoringmatrix[,i])
  }
  
  censoringperc = c(rep(0,length(censoringbin)))
  
  for (i in 1:length(censoringperc)) {
    if (censoringbin[i]==0){
      censoringperc[i] = 0
    }
    else {
      censoringperc[i] = censoringbin[i]/total[total$Var1==i,][,2]
    }
  }
  
  unqcens= as.data.frame(c(0:100)/100)
  colorrange<-colorRampPalette(c("grey100","grey0"))
  unqcens$colour = colorrange(101)
  
  col = c(rep(" ", length(censoringperc)))
  
  for (i in 1:length(censoringperc)) {
        col[i] = unqcens[((round(censoringperc[i],2)*100)+1),2]
  }
  return(col)
}

