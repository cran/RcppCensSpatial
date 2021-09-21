#' Censored Spatial Data Simulation
#'
#' This function simulates censored spatial data with a linear structure for an established censoring rate.
#'
#' @param beta linear regression parameters.
#' @param sigma2 partial sill parameter.
#' @param phi spatial scaling parameter.
#' @param nugget nugget effect parameter.
#' @param x design matrix.
#' @param coords 2D spatial coordinates.
#' @param cens '\code{left}' or '\code{right}' censoring. By default \code{='left'}.
#' @param pcens desired censoring rate. By default \code{=0.10}.
#' @param npred number of simulated data used for cross-validation (Prediction). By default \code{=0}.
#' @param cov.model type of spatial correlation function: '\code{exponential}', '\code{gaussian}',
#' '\code{matern}', and '\code{pow.exp}' for exponential, gaussian, matern, and power exponential, respectively.
#' @param kappa parameter for all spatial correlation functions. For exponential and gaussian
#' \eqn{\kappa=0}, for power exponential \eqn{0 < \kappa <= 2}, and for matern correlation function
#' \eqn{\kappa > 0}.
#'
#' @return If \code{npred > 0}, it returns a list with two datasets: \code{TrainingData} and \code{TestData};
#' otherwise, it returns a data frame with the simulated data.
#'
#' \code{TrainingData}
#' \item{yobs}{generated response vector.}
#' \item{cens}{censoring indicator.}
#' \item{LI}{lower censoring bound.}
#' \item{LS}{upper censoring bound.}
#' \item{xcoord}{x coordinates.}
#' \item{ycoord}{y coordinates.}
#' \item{X}{design matrix.}
#'
#' \code{TestData}
#' \item{yobs}{generated response vector.}
#' \item{xcoord}{x coordinates.}
#' \item{ycoord}{y coordinates.}
#' \item{X}{design matrix.}
#'
#' @author Katherine L. Valeriano, Alejandro Ordonez, Christian E. Galarza and Larissa A. Matos.
#'
#' @examples
#' n = 100
#' set.seed(1000)
#' coords = round(matrix(runif(2*n,0,15),n,2),5)
#' x = cbind(1, rnorm(n))
#' data = rCensSp(c(5,2),2,4,0.70,x,coords,"left",0.10,10,"gaussian",0)
#' data$TrainingData
#' data$TestData

rCensSp = function(beta,sigma2,phi,nugget,x,coords,cens="left",pcens=0.10,npred=0,
                   cov.model="exponential",kappa=0){

  beta = c(beta)
  if (!is.numeric(beta)) stop("beta must be a numeric vector.")
  if (!all(is.finite(beta))) stop("beta must contain only finite values.")

  if (length(c(sigma2))>1 | !is.numeric(sigma2)) stop("sigma2 must be a non-negative number.")
  if (sigma2 <= 0) stop("sigma2 must be non-negative.")

  if (length(c(phi))>1 | !is.numeric(phi)) stop("phi must be a non-negative number.")
  if (phi <= 0) stop("phi must be non-negtaive.")

  if (length(c(nugget))>1 | !is.numeric(nugget)) stop("nugget must be a non-negative number.")
  if (nugget <= 0) stop("nugget must be non-negative.")

  x = as.matrix(x)
  if (!all(c(is.finite(x)))) stop("x must contain only finite values.")
  if (ncol(x) != length(c(beta))) stop("Non-conformable dimensions between x and beta.")

  coords = as.matrix(coords)
  if (!all(c(is.finite(coords)))) stop ("coords must contain only finite values.")
  if (ncol(coords) != 2) stop("coords must contain 2 columns.")
  if (nrow(coords) != nrow(x)) stop("Non-conformable dimensions between coords and x.")

  if (is.null(cens)) stop("cens must be specified.")
  if (cens!="left" & cens!="right") stop("cens should be one of left or right.")

  if (!is.numeric(pcens) | length(c(pcens))>1) stop ("pcens must be a real number in (0,1).")
  if (pcens<=0 | pcens>=1) stop("pcens must be a real number in (0,1).")

  if (!is.numeric(npred) | length(c(npred))>1) stop ("npred must be a positive integer.")
  if (npred<0 | npred>=nrow(x) | npred%%1!=0) stop("npred must be a positive integer.")

  if (is.null(cov.model)) stop("cov.model must be specified.")
  if (cov.model!="matern" & cov.model!="gaussian" & cov.model!="pow.exp" & cov.model!="exponential"){
    stop("cov.model should be one of matern, gaussian, pow.exp or exponential.")}

  if (cov.model!="exponential" & cov.model!="gaussian"){
    if (length(c(kappa))>1 | !is.numeric(kappa)) stop("kappa must be specified.")
    if (cov.model=="pow.exp" & (kappa>2| kappa<=0)) stop("kappa must be a real in (0,2].")
    if (cov.model=="matern" & kappa<=0) stop("kappa must be a real number in (0,Inf).")
    if (cov.model=="matern" & is.infinite(kappa)) stop("kappa must be a real number in (0,Inf).")
  }

  #---------------------------------------------------------------------#
  #                                Outputs                              #
  #---------------------------------------------------------------------#

  out.gen = random.cens(beta,sigma2,phi,nugget,x,coords,cens,pcens,npred,cov.model,kappa)
  return(out.gen)

}
