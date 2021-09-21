
# Estimate parameters using MCEM algorithm
# -----------------------------------------------------------------------------
MCEM_Spatial = function(y,x,cens,LI,LS,coords,init.phi,init.nugget,type.sc,kappa,
                        lower,upper,MaxIter,nMin,nMax,tol,show.SE){
  if (sum(cens)==0) {
    t1 = Sys.time()
    output = Spatial_model(y,x,coords,init.phi,init.nugget,lower,upper,type.sc,kappa,MaxIter,tol,show.SE)
    t2 = Sys.time()
    ptime = difftime(t2,t1,units="m")

  } else {
    t1 = Sys.time()
    output = MCEMspatial(y,x,cens,LI,LS,coords,init.phi,init.nugget,lower,upper,type.sc,kappa,
                        MaxIter,nMin,nMax,tol,show.SE)
    t2 = Sys.time()
    ptime = difftime(t2,t1,units="m")
  }
  likeli = log_likelihood(y,cens,LI,LS,x,output$beta,output$sigma2,output$phi,output$tau2,
                              coords,type.sc,kappa)
  if (show.SE){
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, tau2=output$tau2, EY=output$EY, EYY=output$EYY, SE=output$SE,
                  InfMat=output$InfMat, loglik=likeli$loglik, AIC=likeli$AIC, BIC=likeli$BIC,
                  Iterations=output$Iterations, ptime=ptime, X=x, coord=coords, type=type.sc,
                  kappa=kappa)
  } else {
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, tau2=output$tau2, EY=output$EY, EYY=output$EYY, loglik=likeli$loglik,
                  AIC=likeli$AIC, BIC=likeli$BIC, Iterations=output$Iterations, ptime=ptime, X=x,
                  coord=coords, type=type.sc, kappa=kappa)
  }
  return (output)
}


# Estimate parameters using SAEM algorithm
# -----------------------------------------------------------------------------
SAEM_Spatial = function(y,x,cens,LI,LS,coords,init.phi,init.nugget,type.sc,kappa,
                        lower,upper,MaxIter,M,pc,tol,show.SE){
  if (sum(cens)==0) {
    t1 = Sys.time()
    output = Spatial_model(y,x,coords,init.phi,init.nugget,lower,upper,type.sc,kappa,MaxIter,tol,show.SE)
    t2 = Sys.time()
    ptime = difftime(t2,t1,units="m")

  } else {
    t1 = Sys.time()
    output = SAEMspatial(y,x,cens,LI,LS,coords,init.phi,init.nugget,lower,upper,type.sc,
                         kappa,MaxIter,pc,M,tol,show.SE)
    t2 = Sys.time()
    ptime = difftime(t2,t1,units="m")
  }
  likeli = log_likelihood(y,cens,LI,LS,x,output$beta,output$sigma2,output$phi,output$tau2,
                          coords,type.sc,kappa)
  if (show.SE){
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, tau2=output$tau2, EY=output$EY, EYY=output$EYY, SE=output$SE,
                  InfMat=output$InfMat, loglik=likeli$loglik, AIC=likeli$AIC, BIC=likeli$BIC,
                  Iterations=output$Iterations, ptime=ptime, X=x, coord=coords, type=type.sc,
                  kappa=kappa)
  } else {
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, tau2=output$tau2, EY=output$EY, EYY=output$EYY, loglik=likeli$loglik,
                  AIC=likeli$AIC, BIC=likeli$BIC, Iterations=output$Iterations, ptime=ptime, X=x,
                  coord=coords, type=type.sc, kappa=kappa)
  }
  return (output)
}


# Estimate parameters using EM algorithm
# -----------------------------------------------------------------------------
EM_Spatial = function(y,x,cens,LI,LS,coords,init.phi,init.nugget,type.sc,kappa,
                      lower,upper,MaxIter,tol,show.SE){
  if (sum(cens)==0) {
    t1 = Sys.time()
    output = Spatial_model(y,x,coords,init.phi,init.nugget,lower,upper,type.sc,kappa,MaxIter,tol,show.SE)
    t2 = Sys.time()
    ptime = difftime(t2,t1,units="m")

  } else {
    t1 = Sys.time()
    output = EMspatial(y,x,cens,LI,LS,coords,init.phi,init.nugget,lower,upper,type.sc,kappa,MaxIter,tol,show.SE)
    t2 = Sys.time()
    ptime = difftime(t2,t1,units="m")
  }
  likeli = log_likelihood(y,cens,LI,LS,x,output$beta,output$sigma2,output$phi,output$tau2,
                          coords,type.sc,kappa)
  if (show.SE){
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, tau2=output$tau2, EY=output$EY, EYY=output$EYY, SE=output$SE,
                  InfMat=output$InfMat, loglik=likeli$loglik, AIC=likeli$AIC, BIC=likeli$BIC,
                  Iterations=output$Iterations, ptime=ptime, X=x, coord=coords, type=type.sc,
                  kappa=kappa)
  } else {
    output = list(Theta=output$Theta, theta=output$theta, beta=output$beta, sigma2=output$sigma2,
                  phi=output$phi, tau2=output$tau2, EY=output$EY, EYY=output$EYY, loglik=likeli$loglik,
                  AIC=likeli$AIC, BIC=likeli$BIC, Iterations=output$Iterations, ptime=ptime, X=x,
                  coord=coords, type=type.sc, kappa=kappa)
  }
  return (output)
}


# Log-likelihood function
# -----------------------------------------------------------------------------
log_likelihood = function(y,cc,lower,upper,x,beta,sigma2,phi,tau2,coords,type,kappa){

  mean = x%*%beta
  distM = crossdist(coords)
  Variance = varianceMat(phi,tau2,sigma2,kappa,distM,type)$Sigma

  if (sum(cc)==0) {
    logver = dmvnorm(c(y),c(mean),Variance,TRUE,TRUE)
  } else {
    invObs = solve(Variance[cc==0,cc==0])
    meanC = mean[cc==1] + Variance[cc==1,cc==0]%*%invObs%*%(y[cc==0]-mean[cc==0])
    meanC = as.vector(meanC)
    varC = Variance[cc==1,cc==1] - Variance[cc==1,cc==0]%*%invObs%*%Variance[cc==0,cc==1]
    if (sum(cc) < 10){
      logden2 = log(pmvnorm(lower[cc==1],upper[cc==1],meanC,NULL,varC)[1])
    } else {
      logden2 = (pmvn(lower[cc==1],upper[cc==1],meanC,varC,uselog2=TRUE)[[1]])/log2(exp(1))
    }
    logver = dmvnorm(c(y[cc==0]),c(mean[cc==0]),as.matrix(Variance[cc==0,cc==0]),TRUE,TRUE) + logden2
  }
  results = list(loglik=logver, AIC=-2*logver+2*(length(beta)+3), BIC=-2*logver+(length(beta)+3)*log(length(y)))
  return (results)

}

# Convergence plot
# -----------------------------------------------------------------------------
plot.convergence = function(model){

  Theta = model$Theta
  X = model$X
  q = length(model$beta)
  Iter = nrow(Theta) - 1

  myplot = vector("list",q+3)
  if (all(X[,1]==1)){ namesE <- c(seq(0,(q-1)),0,0,0) } else { namesE <- c(seq(1,q),0,0,0) }
  Theta1 = rbind(namesE,Theta)
  listabeta = as.list(data.frame(Theta1[,1:q]))
  myplot[1:q] = lapply(listabeta, function(.x) ggplot(data.frame(.x[-1]),aes(x=seq(0,Iter),y=.x[-1])) +
                          geom_line(color="dark blue") + labs(x="Iteration", y=bquote(hat(beta)[.(.x[1])])))
  rm(Theta1)
  myplot[[q+1]] = ggplot(data.frame(Theta),aes(x=seq(0,Iter),y=Theta[,q+1])) + geom_line(color="dark blue") + labs(x="Iteration", y=bquote(hat(sigma)^2))
  myplot[[q+2]] = ggplot(data.frame(Theta),aes(x=seq(0,Iter),y=Theta[,q+2])) + geom_line(color="dark blue") + labs(x="Iteration", y=bquote(hat(phi)))
  myplot[[q+3]] = ggplot(data.frame(Theta),aes(x=seq(0,Iter),y=Theta[,q+3])) + geom_line(color="dark blue") + labs(x="Iteration", y=bquote(hat(tau)^2))
  grid.arrange(grobs=myplot,ncol=3)

}

# Prediction in new locations
# -----------------------------------------------------------------------------
predict.new = function(model,x.new,coord.new){

  mediaO = model$X%*%model$beta
  mediaP = x.new%*%model$beta
  coordAll = rbind(model$coord,coord.new)
  pInd = c(rep(0,length(mediaO)),rep(1,length(mediaP)))
  distM = crossdist(coordAll)
  Variance = varianceMat(model$phi,model$tau2,model$sigma2,model$kappa,distM,model$type)$Sigma
  invObs = solve(Variance[pInd==0,pInd==0])
  yPred = mediaP + Variance[pInd==1,pInd==0]%*%invObs%*%(model$EY - mediaO)
  varPred = Variance[pInd==1,pInd==1] - Variance[pInd==1,pInd==0]%*%invObs%*%Variance[pInd==0,pInd==1]
  return (data.frame(xcoord=coord.new[,1],ycoord=coord.new[,2],predValues=yPred, sdPred=matrix(sqrt(diag(varPred)))))

}

# Generating spatial censored data
# -----------------------------------------------------------------------------
random.cens = function(beta,sigma2,phi,tau,x,coords,type.cens,pcens,npred,type.S,kappa){

  x = as.matrix(x)
  n = nrow(as.matrix(x))
  distM = crossdist(coords)
  varianceC = varianceMat(phi,tau,sigma2,kappa,distM,type.S)$Sigma
  mediaC = x%*%beta
  y = t(rmvnorm(1,mediaC,varianceC))

  # Censored data
  y1 = y[1:(n-npred),]
  if (type.cens=="left"){
    cutoff = quantile(y1,pcens)
    cc = (y1 <= cutoff) + 0
    y1[cc==1] = cutoff
    LI = LS = y1
    LI[cc==1] = -Inf
  } else {
    if (type.cens=="right"){
      cutoff = quantile(y1,1-pcens)
      cc = (y1 >= cutoff) + 0
      y1[cc==1] = cutoff
      LI = LS = y1
      LS[cc==1] = Inf
    }
  } # End if
  TrainingData = data.frame(y1,cc,LI,LS,coords[1:(n-npred),],x[1:(n-npred),])
  names(TrainingData) = c("yobs","cens","LI","LS","xcoord","ycoord",paste("X",1:ncol(x),sep=""))

  if (npred > 0){
    y2 = y[(n-npred+1):n,1]
    TestData = data.frame(y2,coords[(n-npred+1):n,],x[(n-npred+1):n,])
    names(TestData) = c("yobs","xcoord","ycoord",paste("X",1:ncol(x),sep=""))
    return(list(TrainingData=TrainingData,TestData=TestData))
  } else {
    return(TrainingData)
  }

}


# Effective range for spatial correlation
#------------------------------------------------------------------------------
matern.eq <- function(x, cor, phi, kappa){
  mat.eq <- 1/(2^(kappa-1)*gamma(kappa))*(x/phi)^(kappa)*besselK(x/phi,kappa) - cor
  return (mat.eq)
}

Effective.range <- function(cor, phi, kappa, Sp.model) {

  if (Sp.model=="exponential"){
    Eff.range <- -phi*log(cor)
  }

  if (Sp.model=="gaussian"){
    Eff.range <- sqrt(-phi^2*log(cor))
  }

  if (Sp.model=="matern"){
    Eff.range <- uniroot(matern.eq,lower=0.0001,upper=(10000*phi),cor=cor,phi=phi,kappa=kappa)$root
  }

  if (Sp.model=="pow.exp"){
    Eff.range <- (-phi^kappa*log(cor))^(1/kappa)
  }

  Eff.range = as.numeric(Eff.range)
  return (Eff.range)
}
