#' Censored Spatial Model Estimation via SAEM Algorithm
#'
#' This function returns the maximum likelihood (ML) estimates of the unknown parameters in Gaussian spatial models
#' with censored/missing responses via the SAEM algorithm. It supports left, right,
#' interval, or missing values in the dependent variable. It also computes the observed information
#' matrix using the method developed by \insertCite{louis1982finding;textual}{RcppCensSpatial}.
#'
#' @param y vector of responses.
#' @param x design matrix.
#' @param cens vector of censoring indicators. For each observation: \code{1} if censored/missing
#' and \code{0} otherwise.
#' @param LI lower limit of detection. For each observation: if non-censored \code{=y},
#' if left-censored/missing \code{=-Inf}, or \code{=LOD} if right/interval-censored.
#' @param LS upper limit of detection. For each observation: if non-censored \code{=y},
#' if right-censored/missing \code{=Inf}, or \code{=LOD} if left/interval-censored.
#' @param coords 2D spatial coordinates.
#' @param init.phi initial value for the spatial scaling parameter.
#' @param init.nugget initial value for the nugget effect parameter.
#' @param type type of spatial correlation function: '\code{exponential}', '\code{gaussian}',
#' '\code{matern}', and '\code{pow.exp}' for exponential, gaussian, matern, and power exponential, respectively.
#' @param kappa parameter for all spatial correlation functions. See \code{\link{CovMat}}.
#' @param lower,upper vectors of lower and upper bounds for the optimization method. If unspecified, the default is
#' \code{c(0.01,0.01)} for lower and \code{c(30,30)} for upper.
#' @param MaxIter maximum number of iterations of the SAEM algorithm. By default \code{=300}.
#' @param M number of Monte Carlo samples for stochastic approximation. By default \code{=20}.
#' @param pc percentage of iterations of the SAEM algorithm with no-memory. By default \code{=0.25}.
#' @param error maximum convergence error. By default \code{=1e-5}.
#' @param show.SE \code{TRUE} or \code{FALSE}. It indicates if the standard errors should be estimated. By default \code{=TRUE}.
#'
#' @details The spatial Gaussian model is given by
#'
#' \eqn{Y = X\beta + \xi},
#'
#' where \eqn{Y} is the \eqn{n x 1} vector of response, \eqn{X} is the \eqn{n x q} design matrix,
#' \eqn{\beta} is the \eqn{q x 1} vector of regression coefficients to be estimated, and \eqn{\xi}
#' is the error term which is normally distributed with zero-mean and covariance matrix
#' \eqn{\Sigma=\sigma^2 R(\phi) + \tau^2 I_n}. We assume that \eqn{\Sigma} is non-singular and
#' \eqn{X} has full rank \insertCite{diggle2007springer}{RcppCensSpatial}.
#'
#' The estimation process was performed via the SAEM \insertCite{delyon1999convergence}{RcppCensSpatial} algorithm.
#' The spatial SAEM algorithm was previously proposed by \insertCite{lachos2017influence;textual}{RcppCensSpatial}
#' and \insertCite{ordonez2018geostatistical;textual}{RcppCensSpatial} and is available in package \code{CensSpatial}.
#' The difference between this package to \code{CensSpatial} is that the random observations are sampled
#' through the slice sampling algorithm available in package \code{relliptical} and the optimization procedure
#' by the \code{roptim} package.
#'
#' This model is also a particular case of the Spatio-temporal model defined by \insertCite{valeriano2021likelihood;textual}{RcppCensSpatial},
#' when the number of temporal observations is equal to one. The computing codes of the Spatio-temporal
#' SAEM algorithm are available in the package \code{StempCens}.
#'
#' @note
#' The SAEM final estimates correspond to the estimates obtained at the last iteration of the
#' algorithm.
#'
#' To fit a regression model for non-censored data, just set \code{cens} as a vector of zeros.
#'
#' Functions \code{print}, \code{summary}, and \code{plot} work for objects of class \code{sclm}.
#'
#' @return The function returns an object of class \code{sclm} which is a list given by:
#'
#' \item{Theta}{estimated parameters in all iterations, \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
#' \item{theta}{final estimation of \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
#' \item{beta}{estimated \eqn{\beta}.}
#' \item{sigma2}{estimated \eqn{\sigma^2}.}
#' \item{phi}{estimated \eqn{\phi}.}
#' \item{tau2}{estimated \eqn{\tau^2}.}
#' \item{EY}{stochastic approximation of the first moment for the truncated normal distribution.}
#' \item{EYY}{stochastic approximation of the second moment for the truncated normal distribution.}
#' \item{SE}{vector of standard errors of \eqn{\theta = (\beta, \sigma^2, \phi, \tau^2)}.}
#' \item{InfMat}{observed information matrix.}
#' \item{loglik}{log-likelihood for the SAEM method.}
#' \item{AIC}{Akaike information criterion.}
#' \item{BIC}{Bayesian information criterion.}
#' \item{Iterations}{number of iterations needed to converge.}
#' \item{ptime}{processing time.}
#' \item{range}{the effective range.}
#'
#' @author Katherine L. Valeriano, Alejandro Ordonez, Christian E. Galarza and Larissa A. Matos.
#'
#' @seealso \code{\link{EM.sclm}}, \code{\link{MCEM.sclm}}, \code{\link{predict.sclm}}
#'
#' @examples
#' # Simulated example: 10% of right-censored observations
#' n = 50   # Test with another values for n
#' set.seed(1000)
#' coords = round(matrix(runif(2*n,0,15),n,2),5)
#' x = cbind(rbinom(n,1,0.50), rnorm(n), rnorm(n))
#' data = rCensSp(c(1,4,-2),2,3,0.50,x,coords,"right",0.10,0,"matern",2)
#'
#' fit = SAEM.sclm(y=data$yobs, x=data[,7:9], cens=data$cens, LI=data$LI,
#'              LS=data$LS, coords=data[,5:6], init.phi=2, init.nugget=1,
#'              type="matern", kappa=2, MaxIter = 20, error=1e-4)
#' summary(fit)
#'
#' \donttest{
#' # Simulated example: censored and missing observations
#' n = 200
#' set.seed(123)
#' coords = round(matrix(runif(2*n,0,20),n,2),5)
#' x = cbind(1, rnorm(n), rexp(n))
#' data = rCensSp(c(1,4,-1),2,4,0.50,x,coords,"left",0.10,0,"exponential",0)
#' data$yobs[c(10,20)] = NA;   data$cens[c(10,20)] = 1
#' data$LI[c(10,20)] = -Inf;   data$LS[c(10,20)] = Inf
#'
#' fit2 = SAEM.sclm(y=data$yobs, x=data[,7:9], cens=data$cens, LI=data$LI,
#'               LS=data$LS, coords=data[,5:6], init.phi=2, init.nugget=1,
#'               type="exponential", MaxIter = 300, error=1e-4)
#' fit2$theta  # Estimates
#' fit2$SE     # Standard error
#' fit2$InfMat # Information matrix
#' plot(fit2)}
#' @references
#' \insertAllCited

SAEM.sclm = function(y,x,cens,LI,LS,coords,init.phi,init.nugget,type="exponential",kappa=0,
                          lower=c(0.01,0.01),upper=c(30,30),MaxIter=300,M=20,pc=0.25,
                          error=1e-5,show.SE=TRUE){

  #---------------------------------------------------------------------#
  #                              Validations                            #
  #---------------------------------------------------------------------#
  y = as.matrix(y);     n = nrow(y)
  x = as.matrix(x);     cens = as.matrix(cens)
  LI = as.matrix(LI);   LS = as.matrix(LS)
  coords = as.matrix(coords)

  if (!is.numeric(y)) stop("y must be a numeric vector.")
  if (ncol(y) > 1) stop("y must have just one column.")

  if (!all(c(is.finite(x)))) stop("x must contain only finite values.")
  if (nrow(x) != n) stop("Non-conformable dimensions between x and y.")

  if (nrow(cens) != n) stop("Non-conformable dimensions between cens and y.")
  if (ncol(cens) > 1) stop("cens must have just one column.")
  if (!all(c(cens)==0 | c(cens)==1)) stop("cens must contain only 0 or 1.")
  if (!all(cens[is.na(y)]==1)) stop("cens must be equal 1 when y is NA.")

  if (nrow(LI)!=n | nrow(LS)!=n) stop("Non-conformable dimensions on LI or LS.")
  if (ncol(LI)>1 | ncol(LS)>1) stop("LI and LS must have just one column.")
  if (any(is.na(LI)) | !is.numeric(LI)) stop("LI is not specified or contains NA.")
  if (any(is.na(LS)) | !is.numeric(LS)) stop("LS is not specified or contains NA.")
  if (any(LI[cens==1]>=LS[cens==1])) stop("LI must be lower or equal than LS.")
  if (sum(is.na(y)) > 0){ LI[is.na(y)] = -Inf;  LS[is.na(y)] = Inf }

  if (!all(c(is.finite(coords)))) stop("coords must contain only finite values.")
  if (nrow(coords)!=n | ncol(coords)!=2) stop("Non-conformable dimensions between coords and y.")

  if (length(c(init.phi))>1 | !is.numeric(init.phi)) stop("Initial value for phi must be provided.")
  if (init.phi <= 0) stop("init.phi must be non-negative.")

  if (length(c(init.nugget))>1 | !is.numeric(init.nugget)) stop("Initial value for nugget effect must be provided.")
  if (init.nugget <= 0) stop("init.nugget must be non-negative.")

  if (is.null(type)) stop("type must be specified.")
  if (type!="matern" & type!="gaussian" & type!="pow.exp" & type!="exponential"){
    stop("type should be one of matern, gaussian, pow.exp or exponential.")}

  if (type!="exponential" & type!="gaussian"){
    if (length(c(kappa))>1 | !is.numeric(kappa)) stop("kappa must be specified.")
    if (type=="pow.exp" & (kappa>2| kappa<=0)) stop("kappa must be a real in (0,2].")
    if (type=="matern" & kappa<=0) stop("kappa must be a real number in (0,Inf).")
    if (type=="matern" & is.infinite(kappa)) stop("kappa must be a real number in (0,Inf).")
  }

  if (length(c(lower))!=2 | length(c(upper))!=2) stop("lower and upper must be vectors of length 2.")
  if (any(is.na(lower)) | !is.numeric(lower)) stop("lower is not specified or contains NA.")
  if (any(is.na(upper)) | !is.numeric(upper)) stop("upper is not specified or contains NA.")
  if (any(lower>=upper)) stop("lower must be smaller than upper.")
  if (any(lower<=0)) stop("Values in lower must be non-negative.")
  if (!all(is.finite(upper))) stop("upper must contain only finite values.")

  if (length(c(MaxIter))>1 | !is.numeric(MaxIter)) stop("MaxIter must be a positive integer.")
  if (MaxIter<=0 | MaxIter%%1!=0) stop("MaxIter must be a positive integer.")

  if (length(c(M))>1 | !is.numeric(M)) stop("M must be a positive integer.")
  if (M<=0 | M%%1!=0) stop("M must be a positive integer.")

  if (length(c(pc))>1 | !is.numeric(pc)) stop("pc must belong to the interval [0,1].")
  if (pc<0 | pc>1) stop("pc must belong to the interval [0,1].")

  if (length(c(error))>1 | !is.numeric(error)) stop("error must be specified.")
  if (error<=0 | error>1) stop("error must belong to the interval (0,1].")

  if (!is.logical(show.SE)) stop("show.SE must be logical (TRUE/FALSE).")

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  out.Sp = SAEM_Spatial(y,x,cens,LI,LS,coords,init.phi,init.nugget,type,kappa,
                                         lower,upper,MaxIter,M,pc,error,show.SE)

  out.Sp$call = match.call()
  out.Sp$range  = Effective.range(0.05,out.Sp$phi,kappa,type)
  out.Sp$show.SE = show.SE
  out.Sp$cens = cens
  out.Sp$MaxIter = MaxIter

  class(out.Sp) <- "sclm"
  return(out.Sp)

}
