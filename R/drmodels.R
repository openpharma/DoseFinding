
## model functions
#' @rdname drmodels
#' @export 
linear <- function(dose, e0, delta){
  e0 + delta * dose
}

#' @rdname drmodels
#' @export 
linlog <- function(dose, e0, delta, off = 1){
  linear(log(dose + off), e0, delta)
}

#' @rdname drmodels
#' @export 
emax <-  function(dose, e0, eMax, ed50){
  e0 + eMax*dose/(ed50 + dose)
}

#' @rdname drmodels
#' @export 
quadratic <- function(dose, e0, b1, b2){
  e0 + b1 * dose + b2 * dose^2
}

#' @rdname drmodels
#' @export 
exponential <- function(dose, e0, e1, delta){
  e0 + e1*(exp(dose/delta) - 1)
}

#' @rdname drmodels
#' @export 
logistic <- function(dose, e0, eMax, ed50, delta){
  e0 + eMax/(1 + exp((ed50 - dose)/delta))
}

#' @rdname drmodels
#' @export 
betaMod <- function(dose, e0, eMax, delta1, delta2, scal){
  xlogx <- function(x) if(x == 0) 0 else x * log(x) # will not be called with vector x
  logMaxDens <- xlogx(delta1) + xlogx(delta2) - xlogx(delta1 + delta2)
  dose <- dose/scal
  e0 + eMax/exp(logMaxDens) * (dose^delta1) * (1 - dose)^delta2
}

#' @rdname drmodels
#' @export 
sigEmax <- function(dose, e0, eMax, ed50, h){
  e0 + eMax * 1 /(1 + (ed50/dose)^h)
}

#' @rdname drmodels
#' @export 
linInt <- function(dose, resp, nodes){
  if(length(nodes) != length(resp))
    stop("\"nodes\" and \"resp\" need to be of same length in \"linInt\"")
  approx(x=nodes, y=resp, xout = dose)$y
}

## gradients of built-in model functions
#' @rdname drmodels
#' @export 
linearGrad <- function(dose, ...){
  cbind(e0=1, delta=dose)
}

#' @rdname drmodels
#' @export 
linlogGrad <- function(dose, off, ...){
  cbind(e0=1, delta=log(dose+off))
}

#' @rdname drmodels
#' @export 
quadraticGrad <- function(dose, ...){
  cbind(e0=1, b1 = dose, b2 = dose^2)
}

#' @rdname drmodels
#' @export 
emaxGrad <- function(dose, eMax, ed50, ...){
  cbind(e0=1, eMax=dose/(ed50 + dose), ed50=-eMax * dose/(dose + ed50)^2)
}

#' @rdname drmodels
#' @export 
exponentialGrad <- function(dose, e1, delta, ...){
  cbind(e0=1, e1=exp(dose/delta)-1, delta=-exp(dose/delta) * dose * e1/delta^2)
}

#' @rdname drmodels
#' @export 
logisticGrad <- function(dose, eMax, ed50, delta, ...){
  den <- 1 + exp((ed50 - dose)/delta)
  g1 <- -eMax * (den - 1)/(delta * den^2)
  g2 <- eMax * (den - 1) * (ed50 - dose)/(delta^2 * den^2)
  cbind(e0=1, eMax=1/den, ed50=g1, delta=g2)
}

#' @rdname drmodels
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
#' @export 
linIntGrad <- function(dose, resp, nodes, ...){
  knts <- c(nodes[1], nodes, nodes[length(nodes)])
  splines::splineDesign(knots=knts, ord=2, x=dose)
}
