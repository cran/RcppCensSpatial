#' Distance Matrix Computation
#'
#' This function computes the Euclidean distance matrix for a set of coordinates.
#'
#' @param coords 2D spatial coordinates.
#'
#' @return The function returns the \eqn{n x n} distance matrix.
#'
#' @author Katherine L. Valeriano, Alejandro Ordonez, Christian E. Galarza and Larissa A. Matos.
#'
#' @examples
#' n = 100
#' set.seed(1000)
#' x = round(runif(n,0,10), 5)     # X coordinate
#' y = round(runif(n,0,10), 5)     # Y coordinate
#' Mdist = dist2Dmatrix(cbind(x, y))

dist2Dmatrix = function(coords){

  if (!is.null(coords)){
    coords = as.matrix(coords)
    if (!all(c(is.finite(coords)))) stop ("coords must contain only finite values.")
    if (ncol(coords) != 2) stop("coords must contain 2 columns.")
    if (nrow(coords) <= 1) stop("coords must be a matrix.")
  } else { stop("coords must be specified.") }

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  distancesM = crossdist(coords)
  out.ST = (distancesM + t(distancesM))/2
  return(out.ST)

}



#' Covariance Matrix for Spatial Models
#'
#' This function computes the spatial variance-covariance matrix considering exponential, gaussian, matern,
#' or power exponential correlation functions.
#'
#' @param phi spatial scaling parameter.
#' @param tau2 nugget effect parameter.
#' @param sigma2 partial sill parameter.
#' @param dist \eqn{n x n} distance matrix.
#' @param type type of spatial correlation function: '\code{exponential}', '\code{gaussian}',
#' '\code{matern}', and '\code{pow.exp}' for exponential, gaussian, matern, and power exponential, respectively.
#' @param kappa parameter for all spatial correlation functions. For exponential and gaussian
#' \eqn{\kappa=0}, for power exponential \eqn{0 < \kappa <= 2}, and for matern correlation function
#' \eqn{\kappa > 0}.
#'
#' @details The spatial covariance matrix is given by
#'
#' \eqn{\Sigma = [Cov(s_i, s_j )] = \sigma^2 R(\phi) + \tau^2 I_n},
#'
#' where \eqn{\sigma^2 > 0} is the partial sill, \eqn{\phi > 0} is the spatial scaling parameter,
#' \eqn{\tau^2} is known as the nugget effect in the geostatistical framework, \eqn{R(\phi)} is the
#' \eqn{n x n} correlation matrix computed from the correlation function, and \eqn{I_n} is the
#' \eqn{n x n} identity matrix.
#'
#' The spatial correlation functions available are:
#' \describe{
#' \item{\strong{Exponential}:}{\eqn{Corr(d) = exp(-d/\phi)},}
#'
#' \item{\strong{Gaussian}:}{\eqn{Corr(d) = exp(-(d/\phi)^2)},}
#'
#' \item{\strong{Matern}:}{\eqn{Corr(d) = 1/(2^(\kappa-1)\Gamma(\kappa))(d/\phi)^\kappa K_\kappa(d/\phi)},}
#'
#' \item{\strong{Power exponential}:}{\eqn{Corr(d) = exp(-(d/\phi)^\kappa)},}
#'
#' where \eqn{d >= 0} is the Euclidean distance between two observations, \eqn{\Gamma(.)} is the gamma
#' function, \eqn{\kappa} is the smoothness parameter, and \eqn{K_\kappa(.)} is the modified Bessel
#' function of the second kind of order \eqn{\kappa}.
#' }
#'
#' @return The function returns the \eqn{n x n} spatial covariance matrix.
#'
#' @author Katherine L. Valeriano, Alejandro Ordonez, Christian E. Galarza and Larissa A. Matos.
#'
#' @seealso \code{\link{EM.sclm}}, \code{\link{SAEM.sclm}}, \code{\link{MCEM.sclm}}, \code{\link{dist2Dmatrix}}
#'
#' @examples
#' # Initial parameter values
#' phi = 5;  tau2 = 0.80;  sigma2 = 2
#' n = 20
#' set.seed(1000)
#' x = round(runif(n,0,10), 5)     # X coordinate
#' y = round(runif(n,0,10), 5)     # Y coordinate
#' Ms = dist2Dmatrix(cbind(x, y))
#' Cov = CovMat(phi, tau2, sigma2, Ms, "exponential", 0)

CovMat = function(phi, tau2, sigma2, dist, type="exponential", kappa=0){

  if (length(c(phi))>1 | !is.numeric(phi)) stop("phi must be specified.")
  if (phi <= 0) stop("The spatial parameter (phi) must be non-negative.")

  if (length(c(tau2))>1 | !is.numeric(tau2)) stop("tau2 must be specified.")
  if (tau2 <= 0) stop("The nugget effect (tau2) must be non-negative.")

  if (length(c(sigma2))>1 | !is.numeric(sigma2)) stop("sigma2 must be specified.")
  if (sigma2 <= 0) stop("The partial sill (sigma2) must be non-negative.")

  dist = as.matrix(dist)
  if (!all(c(is.finite(dist)))) stop ("dist must contain only finite values.")
  if (ncol(dist) != nrow(dist)) stop("Distance matrix must be specified.")
  if (!isSymmetric(dist)) stop("Distance matrix must be symmetric.")

  if (is.null(type)) stop("type must be specified.")
  if (type!="matern" & type !="gaussian" & type != "pow.exp" & type != "exponential"){
    stop("type should be one of matern, gaussian, pow.exp, exponential.")}

  if (type!="exponential" & type!="gaussian"){
    if (length(c(kappa))>1 | !is.numeric(kappa)) stop("kappa must be specified.")
    if (type=="pow.exp" & (kappa > 2| kappa<=0)) stop("kappa must be a real in (0,2].")
    if (type=="matern" & kappa <= 0) stop("kappa must be a real number in (0,Inf).")
    if (type=="matern" & is.infinite(kappa)) stop("kappa must be a real number in (0,Inf).")
  }

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  Covariance = varianceMat(phi, tau2, sigma2, kappa, dist, type)$Sigma
  out.ST = (Covariance + t(Covariance))/2
  return(out.ST)

}



#' Prediction in Spatial Model with censored/missing responses
#'
#' This function performs spatial prediction in a set of new \code{S} spatial locations.
#'
#' @param object object of class '\code{sclm}' given as output of \code{\link{EM.sclm}}, \code{\link{SAEM.sclm}}
#' or \code{\link{MCEM.sclm}} function.
#' @param locPre matrix of coordinates for which prediction is performed.
#' @param xPre matrix of covariates for which prediction is performed.
#' @param ... further arguments passed to or from other methods.
#'
#' @details This function predicts using the Mean Squared Error (MSE) criterion, which takes the conditional
#' expectation E(Y|X) as the best linear predictor.
#'
#' @return The function returns a data frame with:
#' \item{xcoord}{x coordinates.}
#' \item{ycoord}{y coordinates.}
#' \item{predValues}{predicted values.}
#' \item{sdPred}{predicted standard deviations.}
#'
#' @author Katherine L. Valeriano, Alejandro Ordonez, Christian E. Galarza and Larissa A. Matos.
#'
#' @seealso \code{\link{EM.sclm}}, \code{\link{SAEM.sclm}}, \code{\link{MCEM.sclm}}
#'
#' @examples
#' \donttest{
#' n = 120
#' set.seed(1000)
#' coords = round(matrix(runif(2*n,0,15),n,2),5)
#' x = cbind(rbinom(n,1,0.50), rnorm(n), rnorm(n))
#' data = rCensSp(c(1,4,-1),2,3,0.50,x,coords,"left",0.10,20,"exponential",0)
#'
#' # Estimation
#' data1 = data$TrainingData
#' # EM algorithm
#' fit1 = EM.sclm(y=data1$yobs, x=data1[,7:9], cens=data1$cens,LI=data1$LI,
#'              LS=data1$LS, coords=data1[,5:6], init.phi=2.50, init.nugget=1,
#'              type="exponential", show.SE=TRUE, error=1e-4)
#' # SAEM algorithm
#' fit2 = SAEM.sclm(y=data1$yobs,x=data1[,7:9],cens=data1$cens,LI=data1$LI,
#'              LS=data1$LS, coords=data1[,5:6], init.phi=2.50, init.nugget=1,
#'              type="exponential", show.SE=TRUE, error=1e-4)
#' # MCEM algorithm
#' fit3 = MCEM.sclm(y=data1$yobs,x=data1[,7:9],cens=data1$cens,LI=data1$LI,
#'              LS=data1$LS, coords=data1[,5:6], init.phi=2.50, init.nugget=1,
#'              type="exponential", MaxIter=300, show.SE=TRUE, error=1e-4)
#' c(fit1$theta)
#' c(fit2$theta)
#' c(fit3$theta)
#'
#' # Prediction
#' data2 = data$TestData
#' pred1 = predict(fit1, data2[,2:3], data2[,4:6])
#' pred2 = predict(fit2, data2[,2:3], data2[,4:6])
#' pred3 = predict(fit3, data2[,2:3], data2[,4:6])
#'
#' # Cross-validation
#' mean((data2$yobs - pred1$predValues)^2)
#' mean((data2$yobs - pred2$predValues)^2)
#' mean((data2$yobs - pred3$predValues)^2)}
#' @export
predict.sclm = function(object, locPre, xPre, ...){

  if (is.null(object)) stop("object must be specified.")

  if (!is.null(locPre) & !is.null(xPre)){
    locPre = as.matrix(locPre)
    xPre = as.matrix(xPre)
    if (!all(c(is.finite(locPre)))) stop("locPre must contain only finite values.")
    if (!all(c(is.finite(xPre)))) stop("xPre must contain only finite values.")
    if(ncol(xPre)!=length(c(object$beta))) stop("Non-conformable dimensions between xPred and beta.")
    if (nrow(locPre)!=nrow(xPre) | ncol(locPre)!=2) stop("Non-conformable dimensions between locPre and xPre.")
  } else { stop("locPre and xPre must be specified.") }

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  ypred = predict.new(object,xPre,locPre)
  out.ST = ypred

  return(out.ST)

}



#' @export
summary.sclm = function(object, ...){
  cat('----------------------------------------------------------------\n')
  cat('           Censored Linear Spatial Regression Model   \n')
  cat('----------------------------------------------------------------\n')
  cat("Call:\n")
  print(object$call)
  cat('\nEstimated parameters:\n')
  # Estimates
  if (all(object$X[,1]==1)) {
    namesx = paste0('\u03b2',0)
    if (ncol(object$X) > 1){ for (i in 2:ncol(object$X)){ namesx = cbind(namesx, paste0('\u03b2',i-1)) } }
  } else {
    namesx = paste0('\u03b2',1)
    if (ncol(object$X) > 1){ for (i in 2:ncol(object$X)){ namesx = cbind(namesx, paste0('\u03b2',i)) } }
  }
  greeks = c(sigma='\u03c3\u00B2', phi='\u03D5', tau='\u03C4\u00B2')
  lab = c(namesx,greeks)

  if (object$show.SE) {
    tab = round(rbind(c(object$theta),c(object$SE)),4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(c(object$theta)),4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  print(tab)
  cat('\n')
  cat(paste('The effective range is',round(object$range,4),'spatial units.\n'))
  cat('\nModel selection criteria:\n')
  critFin <- c(object$loglik, object$AIC, object$BIC)
  critFin <- round(t(as.matrix(critFin)),2)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC"))
  print(critFin)
  cat('\nDetails:\n')
  cat('Number of censored/missing values =',sum(object$cens),'\n')
  cat('Convergence reached? =',(object$Iterations < object$MaxIter),'\n')
  cat('Iterations =',object$Iterations,'/',object$MaxIter,'\n')
  cat('Processing time =',round(object$ptime,4),units(object$ptime),'\n')
}



#' @export
print.sclm = function(x, ...){
  cat('----------------------------------------------------------------\n')
  cat('           Spatial Censored Linear Regression Model   \n')
  cat('----------------------------------------------------------------\n')
  cat("Call:\n")
  print(x$call)
  cat('\nEstimated parameters:\n')
  # Estimates
  if (all(x$X[,1]==1)) {
    namesx = paste0('\u03b2',0)
    if (ncol(x$X) > 1){ for (i in 2:ncol(x$X)){ namesx = cbind(namesx, paste0('\u03b2',i-1)) } }
  } else {
    namesx = paste0('\u03b2',1)
    if (ncol(x$X) > 1){ for (i in 2:ncol(x$X)){ namesx = cbind(namesx, paste0('\u03b2',i)) } }
  }
  greeks = c(sigma='\u03c3\u00B2', phi='\u03D5', tau='\u03C4\u00B2')
  lab = c(namesx,greeks)

  if (x$show.SE) {
    tab = round(rbind(c(x$theta),c(x$SE)),4)
    colnames(tab) = lab
    rownames(tab) = c("","s.e.")
  } else {
    tab = round(rbind(c(x$theta)),4)
    colnames(tab) = lab
    rownames(tab) = c("")
  }
  print(tab)
  cat('\n')
  cat(paste('The effective range is',round(x$range,4),'spatial units.\n'))
  cat('\nModel selection criteria:\n')
  critFin <- c(x$loglik, x$AIC, x$BIC)
  critFin <- round(t(as.matrix(critFin)),2)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC"))
  print(critFin)
  cat('\nDetails:\n')
  cat('Number of censored/missing values =',sum(x$cens),'\n')
  cat('Convergence reached? =',(x$Iterations < x$MaxIter),'\n')
  cat('Iterations =',x$Iterations,'/',x$MaxIter,'\n')
  cat('Processing time =',round(x$ptime,4),units(x$ptime),'\n')
}



#' @export
plot.sclm = function(x, ...){
  plot.convergence(x)
}
