\name{DR-Models}
\alias{drmodels}
\alias{betaMod}
\alias{emax}
\alias{sigEmax}
\alias{exponential}
\alias{logistic}
\alias{linear}
\alias{linlog}
\alias{quadratic}
\alias{linInt}
\alias{betaModGrad}
\alias{emaxGrad}
\alias{sigEmaxGrad}
\alias{exponentialGrad}
\alias{logisticGrad}
\alias{linearGrad}
\alias{linlogGrad}
\alias{quadraticGrad}
\alias{linIntGrad}
\title{ Built-in dose-response models in DoseFinding }
\description{
  Dose-response model functions and gradients.
  
  Below are the definitions of the model functions:
  
  \bold{Emax model}
  \deqn{
    f(d,\theta)=E_0+E_{max}\frac{d}{ED_{50}+d}}{f(d,theta)=E0+Emax d/(ED50 + d).}

  \bold{Sigmoid Emax Model}
  \deqn{
    f(d,\theta)=E_0+E_{max}\frac{d^h}{ED^h_{50}+d^h}}{f(d,theta)=E0+Emax d^h/(ED50^h + d^h).}
  
  \bold{Exponential Model}
  \deqn{
    f(d,\theta)=E_0+E_1(\exp(d/\delta)-1)}{f(d,theta)=E0+E1 (exp(d/delta)-1).}

  \bold{Beta model}
  \deqn{
    f(d,\theta)=E_0+E_{max}B(\delta_1,\delta_2)(d/scal)^{\delta_1}(1-d/scal)^{\delta_2}
  }{f(d,theta)=E0+Emax B(delta1,delta2)(d/scal)^delta1(1-d/scal)^delta2}
  here
  \deqn{B(\delta_1,\delta_2)=(\delta_1+\delta_2)^{\delta_1+\delta_2}/(\delta_1^{\delta_1}
    \delta_2^{\delta_2})}{B(delta1,delta2)=(delta1+delta2)^(delta1+delta2)/(delta1^delta1
    delta2^delta2).}
  and \eqn{scal}{scal} is a fixed dose scaling parameter.
  
  \bold{Linear Model}
 \deqn{
  f(d,\theta)=E_0+\delta d}{f(d,theta)=E0+delta d.}

  \bold{Linear in log Model}
  \deqn{
    f(d,\theta)=E_0+\delta \log(d + off)}{f(d,theta)=E0+delta log(d +
  off),} here \eqn{off}{off} is a fixed offset parameter.

  \bold{Logistic Model} 
  \deqn{
    f(d, \theta) = E_0 + E_{\max}/\left\{1 + \exp\left[ \left(ED_{50} - d
    \right)/\delta \right] \right\}}{f(d,theta)=E0+Emax/(1 + exp((ED50-d)/delta)).}

  \bold{Quadratic Model}
 \deqn{
  f(d,\theta)=E_0+\beta_1d+\beta_2d^2}{f(d,theta)=E0+beta1 d+beta2 d^2.}

\bold{Linear Interpolation model}\cr
  The linInt model provides linear interpolation at the values
  defined by the nodes vector. In virtually all situations the nodes
  vector is equal to the doses used in the analysis. For example the
  \code{\link{Mods}} and the \code{\link{fitMod}} function automatically
  use the doses that are used in the context of the function call as
  nodes. The guesstimates specified in the \code{\link{Mods}} function
  need to be the treatment effects at the active doses standardized to
  the interval [0,1] (see the examples in the \code{\link{Mods}} function).
}
\usage{
emax(dose, e0, eMax, ed50)
emaxGrad(dose, eMax, ed50, ...)

sigEmax(dose, e0, eMax, ed50, h)
sigEmaxGrad(dose, eMax, ed50, h, ...)

exponential(dose, e0, e1, delta)
exponentialGrad(dose, e1, delta, ...)

quadratic(dose, e0, b1, b2)
quadraticGrad(dose, ...)

betaMod(dose, e0, eMax, delta1, delta2, scal)
betaModGrad(dose, eMax, delta1, delta2, scal, ...)

linear(dose, e0, delta)
linearGrad(dose, ...)

linlog(dose, e0, delta, off = 1)
linlogGrad(dose, off, ...)

logistic(dose, e0, eMax, ed50, delta)
logisticGrad(dose, eMax, ed50, delta, ...)

linInt(dose, resp, nodes)
linIntGrad(dose, resp, nodes, ...)
}
\arguments{
  \item{dose}{ Dose variable }
  \item{e0}{ For most models placebo effect. For logistic model 
    left-asymptote parameter, corresponding to a basal effect level 
     (not the placebo effect) }
  \item{eMax}{ Beta Model: Maximum effect within dose-range\cr
               Emax, sigmoid Emax, logistic Model: Asymptotic maximum effect }
  \item{ed50}{ Dose giving half of the asymptotic maximum effect }
  \item{h}{ Hill parameter, determining the steepness of the model at
    the ED50 }
  \item{e1}{ Slope parameter for exponential model}
  \item{delta}{ Exponential model: Parameter, controlling the convexity
    of the model.\cr
    Linear and linlog model: Slope parameter\cr
    Logistic model: Parameter controlling determining the steepness of the curve}
  \item{delta1}{ delta1 parameter for beta model }
  \item{delta2}{ delta2 parameter for beta model}
  \item{b1}{ first parameter of quadratic model }
  \item{b2}{ second parameter of quadratic model (controls, whether
    model is convex or concave) }
  \item{resp}{ Response values at the nodes for the linInt model}
  \item{off}{ Offset value to avoid problems with dose=0 (treated as a
    fixed value, not estimated) }
  \item{scal}{ Scale parameter (treated as a
    fixed value, not estimated)}
  \item{nodes}{ Interpolation nodes for the linear interpolation for the
    linInt model (treated as a fixed value, not estimated)}
  \item{...}{
    Just included for convenience in the gradient functions, so that for example
    \code{quadratic(dose, e0=0, b1=1, b2=3)} will not throw an error
    (although the gradient of the quadratic model is independent of e0,
    b1 and b2).
  }

}
\details{
  
The \bold{Emax model} is used to represent monotone, concave
dose-response shapes.  To distinguish it from the more general sigmoid
emax model it is sometimes also called hyperbolic emax model.

The \bold{sigmoid Emax} model is an extension of the (hyperbolic) Emax
model by introducing an additional parameter h, that determines the
steepness of the curve at the ed50 value. The sigmoid Emax model
describes monotonic, sigmoid dose-response relationships. In the
toxicology literature this model is also called four-parameter
log-logistic (4pLL) model.

The \bold{quadratic} model is intended to capture a possible
non-monotonic dose-response relationship.

The \bold{exponential model} is intended to capture a possible
sub-linear or a convex dose-response relationship.

The \bold{beta model} is intended to capture non-monotone dose-response
relationships and is more flexible than the quadratic model.  The kernel
of the beta model function consists of the kernel of the density
function of a beta distribution on the interval [0,scal]. The parameter
scal is not estimated but needs to be set to a value larger than the
maximum dose. It can be set in most functions (\samp{fitMod},
\samp{Mods}) via the \samp{addArgs} argument, when omitted a value of
\samp{1.2*(maximum dose)} is used as default, where the maximum dose is
inferred from other input to the respective function.

The \bold{linear in log-dose} model is intended to capture concave
shapes. The parameter \code{off} is not estimated in the code but set to
a pre-specified value. It can be set in most functions (\samp{fitMod},
\samp{Mods}) via the \samp{addArgs} argument, when omitted a value of
\samp{0.01*(maximum dose)} is used as default, where the maximum dose is
inferred from other input to the respective function.

The \bold{logistic model} is intended to capture general monotone, 
sigmoid dose-response relationships. The logistic model and the sigmoid
Emax model are closely related: The sigmoid Emax model is a logistic
model in log(dose).

The \bold{linInt model} provids linear interpolation of the means at the
doses. This can be used as a "nonparametric" estimate of the
dose-response curve, but is probably most interesting for specifying a
"nonparametric" truth during planning and assess how well parametric
models work under a nonparametric truth. For the function \samp{Mods}
and \samp{fitMod} the interpolation \samp{nodes} are selected equal to
the dose-levels specified.

}
\value{
 Response value for model functions or matrix containing the gradient evaluations.
}
\references{ 
MacDougall, J. (2006). Analysis of dose-response studies - Emax model,\emph{in} N. Ting (ed.),
\emph{Dose Finding in Drug Development}, Springer, New York, pp. 127--145

Pinheiro, J. C., Bretz, F. and Branson, M. (2006). Analysis of dose-response studies - modeling 
approaches, \emph{in} N. Ting (ed.). \emph{Dose Finding in Drug Development}, Springer, New York,
pp. 146--171
}
\examples{
## some quadratic example shapes
quadModList <- Mods(quadratic = c(-0.5, -0.75, -0.85, -1), doses = c(0,1))
plot(quadModList)

## some emax example shapes
emaxModList <- Mods(emax = c(0.02,0.1,0.5,1), doses = c(0,1))
plot(emaxModList)
## example for gradient
emaxGrad(dose = (0:4)/4, eMax = 1, ed50 = 0.5)

## some sigmoid emax example shapes
sigEmaxModList <- Mods(sigEmax = rbind(c(0.05,1), c(0.15,3), c(0.4,8),
                       c(0.7,8)), doses = c(0,1))
plot(sigEmaxModList)
sigEmaxGrad(dose = (0:4)/4, eMax = 1, ed50 = 0.5, h = 8)

## some exponential example shapes
expoModList <- Mods(exponential = c(0.1,0.25,0.5,2), doses=c(0,1))
plot(expoModList)
exponentialGrad(dose = (0:4)/4, e1 = 1, delta = 2)

## some beta model example shapes
betaModList <- Mods(betaMod = rbind(c(1,1), c(1.5,0.75), c(0.8,2.5),
                    c(0.4,0.9)), doses=c(0,1), addArgs=list(scal = 1.2))
plot(betaModList)
betaModGrad(dose = (0:4)/4, eMax = 1, delta1 = 1, delta2 = 1, scal = 5)

## some logistic model example shapes
logistModList <- Mods(logistic = rbind(c(0.5,0.05), c(0.5,0.15),
                      c(0.2,0.05), c(0.2,0.15)), doses=c(0,1))
plot(logistModList)
logisticGrad(dose = (0:4)/4, eMax = 1, ed50 = 0.5, delta = 0.05)

## some linInt shapes
genModList <- Mods(linInt = rbind(c(0.5,1,1),
                      c(0,1,1), c(0,0,1)), doses=c(0,0.5,1,1.5))
plot(genModList)
linIntGrad(dose = (0:4)/4, resp=c(0,0.5,1,1,1), nodes=(0:4)/4)

}

\seealso{ \code{\link{fitMod}} }


