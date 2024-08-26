
## model functions
#' @rdname drmodels
#' @param dose Dose variable
#' @param e0 For most models placebo effect. For logistic model left-asymptote
#' parameter, corresponding to a basal effect level (not the placebo effect)
#' @param eMax Beta Model: Maximum effect within dose-range\cr Emax, sigmoid
#' Emax, logistic Model: Asymptotic maximum effect
#' @param ed50 Dose giving half of the asymptotic maximum effect
#' @usage NULL
#' @export 
emax <-  function(dose, e0, eMax, ed50){
  e0 + eMax*dose/(ed50 + dose)
}

#' @rdname drmodels
#' @inheritParams emax
#' @param ...  Just included for convenience in the gradient functions, so that
#' for example \code{quadratic(dose, e0=0, b1=1, b2=3)} will not throw an error
#' (although the gradient of the quadratic model is independent of e0, b1 and
#' b2).
#' @usage NULL
#' @export 
emaxGrad <- function(dose, eMax, ed50, ...){
  cbind(e0=1, eMax=dose/(ed50 + dose), ed50=-eMax * dose/(dose + ed50)^2)
}

#' @rdname drmodels
#' @inheritParams emax
#' @param h Hill parameter, determining the steepness of the model at the ED50
#' @usage NULL
#' @export 
sigEmax <- function(dose, e0, eMax, ed50, h){
  e0 + eMax * 1 /(1 + (ed50/dose)^h)
}

#' @rdname drmodels
#' @inheritParams sigEmax
#' @usage NULL
#' @export 
sigEmaxGrad <- function(dose, eMax, ed50, h, ...){
  lg2 <- function(x) {l<-x; l[x==0] <- 0; l[x!=0] <- log(x[x!=0]); l}
  a <-  1 / (1 + (dose/ed50)^h)
  g1 <- 1 / (1 + (ed50/dose)^h)
  g2 <- -(h * eMax / ed50) * g1 * a
  g3 <- eMax * lg2(dose / ed50) * g1 * a
  cbind(e0=1, eMax=g1, ed50=g2, h=g3)
}

#' @rdname drmodels
#' @inheritParams emax
#' @param e1 Slope parameter for exponential model
#' @param delta Exponential model: Parameter, controlling the convexity of the
#' model.\cr Linear and linlog model: Slope parameter\cr Logistic model:
#' Parameter controlling determining the steepness of the curve
#' @usage NULL
#' @export 
exponential <- function(dose, e0, e1, delta){
  e0 + e1*(exp(dose/delta) - 1)
}

#' @rdname drmodels
#' @inheritParams exponential
#' @usage NULL
#' @export 
exponentialGrad <- function(dose, e1, delta, ...){
  cbind(e0=1, e1=exp(dose/delta)-1, delta=-exp(dose/delta) * dose * e1/delta^2)
}

#' @rdname drmodels
#' @inheritParams emax
#' @param b1 first parameter of quadratic model
#' @param b2 second parameter of quadratic model (controls, whether model is
#' convex or concave)
#' @usage NULL
#' @export 
quadratic <- function(dose, e0, b1, b2){
  e0 + b1 * dose + b2 * dose^2
}

#' @rdname drmodels
#' @inheritParams quadratic
#' @usage NULL
#' @export 
quadraticGrad <- function(dose, ...){
  cbind(e0=1, b1 = dose, b2 = dose^2)
}

#' @rdname drmodels
#' @inheritParams emax
#' @param delta1 delta1 parameter for beta model
#' @param delta2 delta2 parameter for beta model
#' @param scal Scale parameter (treated as a fixed value, not estimated)
#' @usage NULL
#' @export 
betaMod <- function(dose, e0, eMax, delta1, delta2, scal){
  xlogx <- function(x) if(x == 0) 0 else x * log(x) # will not be called with vector x
  logMaxDens <- xlogx(delta1) + xlogx(delta2) - xlogx(delta1 + delta2)
  dose <- dose/scal
  e0 + eMax/exp(logMaxDens) * (dose^delta1) * (1 - dose)^delta2
}

#' @rdname drmodels
#' @inheritParams betaMod
#' @usage NULL
#' @export 
betaModGrad <- function(dose, eMax, delta1, delta2, scal, ...){
  lg2 <- function(x) {l<-x; l[x==0] <- 0; l[x!=0] <- log(x[x!=0]); l}
  xlogx <- function(x) if(x == 0) 0 else x * log(x) # will not be called with vector x
  dose <- dose/scal
  if(any(dose > 1)) {
    stop("doses cannot be larger than scal in betaModel")
  }
  logMaxDens <- xlogx(delta1) + xlogx(delta2) - xlogx(delta1 + delta2)
  g1 <- ((dose^delta1) * (1 - dose)^delta2)/exp(logMaxDens)
  g2 <- g1 * eMax * (lg2(dose) + lg2(delta1 + delta2) - lg2(delta1))
  g3 <- g1 * eMax * (lg2(1 - dose) + lg2(delta1 + delta2) - lg2(delta2))
  cbind(e0=1, eMax=g1, delta1=g2, delta2=g3)
}

#' @rdname drmodels
#' @inheritParams exponential
#' @usage NULL
#' @export 
linear <- function(dose, e0, delta){
  e0 + delta * dose
}

#' @rdname drmodels
#' @inheritParams linear
#' @usage NULL
#' @export 
linearGrad <- function(dose, ...){
  cbind(e0=1, delta=dose)
}

#' @rdname drmodels
#' @inheritParams exponential
#' @param off Offset value to avoid problems with dose=0 (treated as a fixed
#' value, not estimated)
#' @usage NULL
#' @export 
linlog <- function(dose, e0, delta, off = 1){
  linear(log(dose + off), e0, delta)
}

#' @rdname drmodels
#' @inheritParams linlog
#' @usage NULL
#' @export 
linlogGrad <- function(dose, off, ...){
  cbind(e0=1, delta=log(dose+off))
}

#' @rdname drmodels
#' @inheritParams emax
#' @param delta Exponential model: Parameter, controlling the convexity of the
#' model.\cr Linear and linlog model: Slope parameter\cr Logistic model:
#' Parameter controlling determining the steepness of the curve
#' @usage NULL
#' @export 
logistic <- function(dose, e0, eMax, ed50, delta){
  e0 + eMax/(1 + exp((ed50 - dose)/delta))
}

#' @rdname drmodels
#' @inheritParams logistic
#' @usage NULL
#' @export 
logisticGrad <- function(dose, eMax, ed50, delta, ...){
  den <- 1 + exp((ed50 - dose)/delta)
  g1 <- -eMax * (den - 1)/(delta * den^2)
  g2 <- eMax * (den - 1) * (ed50 - dose)/(delta^2 * den^2)
  cbind(e0=1, eMax=1/den, ed50=g1, delta=g2)
}

#' @rdname drmodels
#' @inheritParams emax
#' @param resp Response values at the nodes for the linInt model
#' @param nodes Interpolation nodes for the linear interpolation for the linInt
#' model (treated as a fixed value, not estimated)
#' @usage NULL
#' @export 
linInt <- function(dose, resp, nodes){
  if(length(nodes) != length(resp))
    stop("\"nodes\" and \"resp\" need to be of same length in \"linInt\"")
  approx(x=nodes, y=resp, xout = dose)$y
}

#' @rdname drmodels
#' @inheritParams linInt
#' @usage NULL
#' @export 
linIntGrad <- function(dose, resp, nodes, ...){
  knts <- c(nodes[1], nodes, nodes[length(nodes)])
  splines::splineDesign(knots=knts, ord=2, x=dose)
}
