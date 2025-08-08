#' Simulate dose-response data
#' 
#' Functions to generate a dataset with simulated normally distributed, binary or count responses, according to a pre-specified
#' dose-response model (or mean vector)
#'
#' @name genDFdata
#' @param model Character string giving the name of a model function. 
#' @param argsMod vector with arguments for the model function
#' @param doses vector of dose levels to be used
#' @param n vector of group sample sizes
#' @param mu mean vector to use if no model is specified
#'
#' @return A data frame with two columns corresponding to dose and simulated response
NULL

#' @rdname genDFdata
#' @param sigma standard deviation. Can be a single value or a vector to simulate heteroscedasticity across doses
#'
#' @export
#'
#' @examples
#' # Simulate normally distributed data under Emax model
#' genDFdataNorm("emax", c(e0 = 0.2, eMax = 1, ed50 = 0.05), 
#'               doses = c(0,0.05,0.2,0.6,1), n = 20, sigma = 1)
#' # use fixed mean vector
#' genDFdataNorm(mu = 1:5, doses = 0:4, n = c(20, 20, 10, 5, 1), sigma = 0.2)
#' # simulate heteroscedastic data
#' genDFdataNorm("emax", c(e0 = 0.2, eMax = 1, ed50 = 0.05), 
#'               doses = c(0,0.05,0.2,0.6,1), n = 20, sigma = seq(0.6, 1, 0.1))
#' 
genDFdataNorm <- function(model, argsMod, doses, n, sigma, mu = NULL){
  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) 
    n <- rep(n, nD)
  dose <- rep(dose, n)
  if(!missing(model)){
    args <- c(list(dose), argsMod)
    mu <- do.call(model, args)
  } else if(!is.null(mu)){
    if(length(doses) != length(mu)){
      stop("'mu' needs to be of the same length as doses")
    }
    mu <- rep(mu,  n)
  } else {
    stop("either 'model' or 'mu' needs to be specified")
  }
  if(length(sigma) > 1){
    if((length(doses) != length(sigma)))
      stop("'sigma' needs to be of length 1 or the same length as doses")
    sigma <- rep(sigma, n)
  }
  data.frame(dose = dose, 
             resp = mu + rnorm(sum(n), sd = sigma))
}
 
#' @rdname genDFdata
#' @param scale specifies scale on which the dose-response model or mean vector is specified. The defaults are logit scale for 
#'        binary data and log scale for count data.
#'
#' @export
#'
#' @examples
#' # Simulate binary data under Emax model
#' genDFdataBin("emax", c(e0 = log(0.4/0.6), eMax = log(0.6/0.4), ed50 = 0.05), 
#'              doses = c(0,0.05,0.2,0.6,1), n = 20)
#' # use fixed mean vector on identity scale
#' genDFdataBin(mu = c(0.2, 0.25, 0.4, 0.5, 0.6), 
#'              doses = 0:4, n = c(20, 20, 10, 5, 1), scale = "identity")
#' 
genDFdataBin <- function(model, argsMod, doses, n, mu = NULL, scale = c("logit", "identity", "probit", "cloglog")){
  scale <- match.arg(scale)
  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) 
    n <- rep(n, nD)
  dose <- rep(dose,  n)
  if(!missing(model)){
    args <- c(list(dose), argsMod)
    mu <- do.call(model, args)
  } else if(!is.null(mu)){
    if(length(doses) != length(mu)){
      stop("'mu' needs to be of the same length as doses")
    }
    mu <- rep(mu,  n)
  } else {
    stop("either 'model' or 'mu' needs to be specified")
  }
  if(scale == "logit")
    mu <- exp(mu)/(1 + exp(mu))
  if(scale == "probit")
    mu <- pnorm(mu)
  if(scale == "cloglog")
    mu <- 1-exp(-exp(mu))
  
  data.frame(dose = dose, 
             resp = rbinom(sum(n), 1, mu))
}

#' @rdname genDFdata
#' @param theta parameter controlling the overdispersion of the negative binomial (see [MASS::rnegbin()]). If not 
#' specified poisson distribution will be used.
#'
#' @export
#'
#' @examples
#' # Simulate count data under Emax model
#' genDFdataCount("emax", c(e0 = log(0.4), eMax = log(0.6), ed50 = 0.05), 
#'                doses = c(0,0.05,0.2,0.6,1), n = 20)
#' # negative binomially distributed responses
#' genDFdataCount("emax", c(e0 = log(0.4), eMax = log(0.6), ed50 = 0.05), 
#'                doses = c(0,0.05,0.2,0.6,1), n = 20, theta = 0.2)
#' # use fixed mean vector on identify scale
#' genDFdataCount(mu = c(0.2, 0.25, 0.4, 0.5, 0.6), 
#'                doses = 0:4, n = c(20, 20, 10, 5, 1), scale = "identity")
genDFdataCount <- function(model, argsMod, doses, n, mu = NULL, theta = NULL, scale = c("log", "identity")){
  scale <- match.arg(scale)
  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) 
    n <- rep(n, nD)
  dose <- rep(dose,  n)
  if(!missing(model)){
    args <- c(list(dose), argsMod)
    mu <- do.call(model, args)
  } else if(!is.null(mu)){
    if(length(doses) != length(mu)){
      stop("'mu' needs to be of the same length as doses")
    }
    mu <- rep(mu,  n)
  } else {
    stop("either 'model' or 'mu' needs to be specified")
  }
  if(scale == "log")
    mu <- exp(mu)
  if(!missing(theta)){
    out <- data.frame(dose = dose, resp = MASS::rnegbin(sum(n), mu, theta))
  }
  else{
    out <- data.frame(dose = dose, resp = rpois(sum(n), mu))
  }
  out
}