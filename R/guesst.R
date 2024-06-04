
#' Calculate guesstimates based on prior knowledge
#'
#' Calculates guesstimates for standardized model parameter(s) using the general approach described in Pinheiro et al.
#' (2006).
#'
#' Calculates guesstimates for the parameters \eqn{\theta_2}{theta2} of the standardized model function based on the
#' prior expected percentage of the maximum effect at certain dose levels. Note that this function should be used
#' together with the \code{\link{plot.Mods}} function to ensure that the guesstimates are reflecting the prior beliefs.
#'
#' For the logistic and sigmoid emax models at least two pairs (d,p) need to be specified.
#'
#' For the beta model the dose at which the maximum effect occurs (dMax) has to be specified in addition to the (d,p)
#' pair.
#'
#' For the exponential model the maximum dose administered (Maxd) needs to be specified in addition to the (d,p) pair.
#'
#' For the quadratic model one (d,p) pair is needed. It is advisable to specify the location of the maximum within the
#' dose range with this pair.
#'
#' For the emax, sigmoid Emax and logistic model one can choose between a local and an asymptotic version. In the local
#' version one explicitly forces the standardized model function to pass through the specified points (d,p). For the
#' asymptotic version it assumed that the standardized model function is equal to 1 at the largest dose (this is the
#' approach described in Pinheiro et al. (2006)). If the local version is used, convergence problems with the underlying
#' nonlinear optimization can occur.
#'
#' @param d Vector containing dose value(s).
#' @param p Vector of expected percentages of the maximum effect achieved at d.
#' @param model Character string. Should be one of "emax", "exponential", "quadratic", "betaMod", "sigEmax", "logistic".
#' @param less Logical, only needed in case of quadratic model.  Determines if d is smaller (\samp{less=TRUE}) or larger
#'   (\samp{less=FALSE}) than dopt (see Pinheiro et al. (2006) for details).
#' @param local Logical indicating whether local or asymptotic version of guesstimate should be derived (defaults to
#'   \samp{FALSE}).  Only needed for emax, logistic and sigEmax model.  When \samp{local=TRUE} the maximum dose must be
#'   provided via \samp{Maxd}.
#' @param dMax Dose at which maximum effect occurs, only needed for the beta model
#' @param Maxd Maximum dose to be administered in the trial
#' @param scal Scale parameter, only needed for the beta model
#' @return Returns a numeric vector containing the guesstimates.
#' @seealso \code{\link{emax}}, \code{\link{logistic}}, \code{\link{betaMod}}, \code{\link{sigEmax}},
#'   \code{\link{quadratic}}, \code{\link{exponential}}, \code{\link{plot.Mods}}
#' @references Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An R
#' Package for the Design and Analysis of Dose-Finding Studies, \emph{Journal
#' of Statistical Software}, \bold{29}(7), 1--23
#'
#'   Pinheiro, J. C., Bretz, F., and Branson, M. (2006). Analysis of dose-response studies - modeling approaches,
#'   \emph{in} N. Ting (ed.), \emph{Dose Finding in Drug Development}, Springer, New York, pp. 146--171
#' @examples
#'
#' ## Emax model
#' ## Expected percentage of maximum effect: 0.8 is associated with
#' ## dose 0.3 (d,p)=(0.3, 0.8), dose range [0,1]
#' emx1 <- guesst(d=0.3, p=0.8, model="emax")
#' emax(0.3,0,1,emx1)
#' ## local approach
#' emx2 <- guesst(d=0.3, p=0.8, model="emax", local = TRUE, Maxd = 1)
#' emax(0.3,0,1,emx2)/emax(1,0,1,emx2)
#' ## plot models
#' models <- Mods(emax=c(emx1, emx2), doses=c(0,1))
#' plot(models)
#'
#' ## Logistic model
#' ## Select two (d,p) pairs (0.2, 0.6) and (0.2, 0.95)
#' lgc1 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "logistic")
#' logistic(c(0.2,0.6), 0, 1, lgc1[1], lgc1[2])
#' ## local approach
#' lgc2 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "logistic",
#'                local = TRUE, Maxd = 1)
#' r0 <- logistic(0, 0, 1, lgc2[1], lgc2[2])
#' r1 <- logistic(1, 0, 1, lgc2[1], lgc2[2])
#' (logistic(c(0.2,0.6), 0, 1, lgc2[1], lgc2[2])-r0)/(r1-r0)
#' ## plot models
#' models <- Mods(logistic = rbind(lgc1, lgc2), doses=c(0,1))
#' plot(models)
#'
#' ## Beta Model
#' ## Select one pair (d,p): (0.4,0.8)
#' ## dose, where maximum occurs: 0.8
#' bta <- guesst(d=0.4, p=0.8, model="betaMod", dMax=0.8, scal=1.2, Maxd=1)
#' ## plot
#' models <- Mods(betaMod = bta, doses=c(0,1), addArgs = list(scal = 1.2))
#' plot(models)
#'
#' ## Sigmoid Emax model
#' ## Select two (d,p) pairs (0.2, 0.6) and (0.2, 0.95)
#' sgE1 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "sigEmax")
#' sigEmax(c(0.2,0.6), 0, 1, sgE1[1], sgE1[2])
#' ## local approach
#' sgE2 <- guesst(d = c(0.2, 0.6), p = c(0.2, 0.95), "sigEmax",
#'                local = TRUE, Maxd = 1)
#' sigEmax(c(0.2,0.6), 0, 1, sgE2[1], sgE2[2])/sigEmax(1, 0, 1, sgE2[1], sgE2[2])
#' models <- Mods(sigEmax = rbind(sgE1, sgE2), doses=c(0,1))
#' plot(models)
#'
#' ## Quadratic model
#' ## For the quadratic model it is assumed that the maximum effect occurs at
#' ## dose 0.7
#' quad <- guesst(d = 0.7, p = 1, "quadratic")
#' models <- Mods(quadratic = quad, doses=c(0,1))
#' plot(models)
#'
#' ## exponential model
#' ## (d,p) = (0.8,0.5)
#' expo <- guesst(d = 0.8, p = 0.5, "exponential", Maxd=1)
#' models <- Mods(exponential = expo, doses=c(0,1))
#' plot(models)
#' @export
guesst <- function(d, p, model = c("emax", "exponential", "logistic", "quadratic",
                   "betaMod", "sigEmax"), less = TRUE,  local = FALSE, dMax, Maxd, scal){
  model <- match.arg(model)
  if(any(p <= 0) | any(p > 1)) stop("must have 0 < p <= 1")
  if(model == "emax"){
    if(!local){
      return(c(ed50 = mean(d * (1 - p)/p)))
    } else {
      if (any(p <= d/Maxd))
        stop("must have p > d/Maxd, for local version")
      val <- (d/p-d)/(1-d/(Maxd*p))
      return(c(ed50=mean(val)))
    }
  }
  if(model == "exponential"){
    if(any(p >= d/Maxd)) stop("must have p < d/Maxd")
    init <- d/log(1 + p)
    fooexp <- function(delta, d, p, Maxd){
      sum((exponential(d, 0, 1, delta)/
           exponential(Maxd, 0, 1, delta) - p)^2)           
    }
    val <- optimize(fooexp, c(0, 2*Maxd), d=d, p=p, Maxd=Maxd)$minimum
    return(c(delta = mean(val)))
  }
  if(model == "logistic"){
    if(length(d) == 1) {
      stop("logistic model needs at least two pairs (d,p)")
    }
    logit <- function(p) log(p/(1-p)) 
    if(length(d) == 2) {
      ed50  <- diff(rev(d)*logit(p))/diff(logit(p))
      delta <- diff(d)/diff(logit(p))
      res <- c(ed50 = ed50, delta = delta)
    } else {
      m <- lm(logit(p)~d)
      par <- coef(m)
      names(par) <- NULL
        res <- c(ed50 = -par[1]/par[2], delta = 1/par[2])
    }
    if(local){
      foolog <- function(par, d, p, Maxd) {
        e0 <- logistic(0,0,1,par[1],par[2])
        sum(((logistic(d,0,1,par[1],par[2]) - e0)/
             (logistic(Maxd,0,1,par[1],par[2])-e0)-p)^2)
      }
      res <- try(optim(par=res, fn=foolog, d=d, p=p, Maxd=Maxd))
      if(res$convergence > 0)
        stop("cannot find guesstimates for specified values")
      else res <- res$par
    }
    if(res[1] < 0)
      message("Message: specified values lead to negative ed50, which should be positive")
    return(res)
  }
  if(model == "quadratic"){
    aux <- sqrt(1 - p)
    if (less){
      return(c(delta = mean(-(1 - aux)/(2 * d))))
    } else {
      return(c(delta = mean(-(1 + aux)/(2 * d))))
    }
  }
  if(model == "betaMod"){
    if(scal <= dMax)
      stop("scal needs to be larger than dMax to calculate guesstimate")
    if(dMax > Maxd)
      stop("dose with maximum effect (dMax) needs to be smaller than maximum dose (Maxd)")
    k <- dMax/(scal-dMax)
    val <- d^k*(scal-d)/(dMax^k*(scal-dMax))
    beta <- log(p)/(log(val))
    return(c(delta1 = mean(k*beta), delta2 = mean(beta)))
  }
  if(model == "sigEmax"){
    if(length(d) == 1) {
      stop("sigmoid Emax model needs at least two pairs (d,p)")
    }
    if(length(d) == 2){
      num <- log((p[1]*(1-p[2]))/(p[2]*(1-p[1])))
      h <- num/log(d[1]/d[2])
      ed50 <- ((1-p[1])/p[1])^(1/h)*d[1]
      res <- c(ed50=ed50, h=h)
    } else {
      y <- log((1-p)/p)
      x <- log(d)
      par <- coef(lm(y~x))
      names(par) <- NULL
      res <- c(ed50 = exp(par[1]/-par[2]), delta = -par[2])
    }
    if(local) {
      fooSE <- function(par, d, p, Maxd) {
        sum((sigEmax(d,0,1,par[1],par[2])/
             sigEmax(Maxd,0,1,par[1],par[2])-p)^2)
      }
      res <- try(optim(par=res, fn=fooSE, d=d, p=p, Maxd=Maxd))
      if(res$convergence > 0)
        stop("cannot find guesstimates for specified values")
      else res <- res$par
    }
    if(res[1] < 0)
      message("Message: specified values lead to negative ed50, which should be positive")
    return(res)
  }
}
