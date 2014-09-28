\name{guesst}
\alias{guesst}
\title{ Calculate guesstimates based on prior knowledge }
\description{
Calculates guesstimates for standardized model parameter(s) 
using the general approach described in Pinheiro et al. (2006). 
}
\usage{
guesst(d, p, model = c("emax", "exponential", "logistic", "quadratic",
       "betaMod", "sigEmax"), less = TRUE,  local = FALSE,
       dMax, Maxd, scal)
}
\arguments{
  \item{d}{ Vector containing dose value(s). }
  \item{p}{ Vector of expected percentages of the maximum effect achieved at d. }
  \item{model}{ Character string. Should be one of "emax", "exponential",
               "quadratic", "betaMod", "sigEmax", "logistic". }
  \item{less}{ Logical, only needed in case of quadratic model.
            Determines if d is smaller (\samp{less=TRUE}) or larger (\samp{less=FALSE})
            than dopt (see Pinheiro et al. (2006) for details).}
  \item{local}{ Logical indicating whether local or asymptotic version
           of guesstimate should be derived (defaults to \samp{FALSE}).
           Only needed for emax, logistic and sigEmax model. 
           When \samp{local=TRUE} the maximum dose must be provided via 
           \samp{Maxd}.}
  \item{dMax}{ Dose at which maximum effect occurs, only needed for the beta model }
  \item{Maxd}{ Maximum dose to be administered in the trial }
  \item{scal}{ Scale parameter, only needed for the beta model }
}
\details{
Calculates guesstimates for the parameters \eqn{\theta_2}{theta2} of the
standardized model function based on the prior expected percentage of
the maximum effect at certain dose levels. Note that this function
should be used together with the \code{\link{plot.Mods}} function to
ensure that the guesstimates are reflecting the prior beliefs.

For the logistic and sigmoid emax models at least two pairs (d,p) need to
be specified. 

For the beta model the dose at which the maximum effect occurs (dMax)
has to be specified in addition to the (d,p) pair.

For the exponential model the maximum dose administered (Maxd) needs to 
be specified in addition to the (d,p) pair.

For the quadratic model one (d,p) pair is needed. It is advisable to
specify the location of the maximum within the dose range with this
pair.

For the emax, sigmoid Emax and logistic model one can choose between a local
and an asymptotic version. In the local version one explicitly forces the
standardized model function to pass through the specified points (d,p). For the
asymptotic version it assumed that the standardized model function is equal to 1
at the largest dose (this is the approach described in Pinheiro et al. (2006)).
If the local version is used, convergence problems
with the underlying nonlinear optimization can occur.
}
\value{
Returns a numeric vector containing the guesstimates.
}

\references{ 
Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An            
R Package for the Design and Analysis of Dose-Finding        
Studies, \emph{Journal of Statistical Software}, \bold{29}(7), 1--23

Pinheiro, J. C., Bretz, F., and Branson, M. (2006). Analysis of dose-response studies - modeling
approaches, \emph{in} N. Ting (ed.), \emph{Dose Finding in Drug Development}, Springer, New York,
pp. 146--171
}

\examples{
## Emax model
## Expected percentage of maximum effect: 0.8 is associated with
## dose 0.3 (d,p)=(0.3, 0.8), dose range [0,1]
emx1 <- guesst(d=0.3, p=0.8, model="emax")
## local approach
emx2 <- guesst(d=0.3, p=0.8, model="emax", local = TRUE, Maxd = 1)
## plot models
models <- Mods(emax=c(emx1, emx2), doses=c(0,1))
plot(models)

## Logistic model
## Select two (d,p) pairs (0.2, 0.5) and (0.6, 0.95)
lgc1 <- guesst(d = c(0.2, 0.5), p = c(0.6, 0.95), "logistic")
## local approach
lgc2 <- guesst(d = c(0.2, 0.5), p = c(0.6, 0.95), "logistic", 
               local = TRUE, Maxd = 1)
## plot models
models <- Mods(logistic = rbind(lgc1, lgc2), doses=c(0,1))
plot(models)

## Beta Model
## Select one pair (d,p): (0.5,0.5)
## dose, where maximum occurs: 0.8
bta <- guesst(d=0.5, p=0.5, model="betaMod", dMax=0.8, scal=1.2, Maxd=1)
## plot
models <- Mods(betaMod = bta, doses=c(0,1), addArgs = list(scal = 1.2))
plot(models)

## Sigmoid Emax model
## Select two (d,p) pairs (0.2, 0.5) and (0.6, 0.95)
sgE1 <- guesst(d = c(0.2, 0.5), p = c(0.6, 0.95), "sigEmax")
## local approach
sgE2 <- guesst(d = c(0.2, 0.5), p = c(0.6, 0.95), "sigEmax",
               local = TRUE, Maxd = 1)
models <- Mods(sigEmax = rbind(sgE1, sgE2), doses=c(0,1))
plot(models)

## Quadratic model
## For the quadratic model it is assumed that the maximum effect occurs at 
## dose 0.7
quad <- guesst(d = 0.7, p = 1, "quadratic")
models <- Mods(quadratic = quad, doses=c(0,1))
plot(models)

## exponential model
## (d,p) = (0.8,0.5)
expo <- guesst(d = 0.8, p = 0.5, "exponential", Maxd=1)
models <- Mods(exponential = expo, doses=c(0,1))
plot(models)
}
\seealso{ \code{\link{emax}}, \code{\link{logistic}}, \code{\link{betaMod}},
          \code{\link{sigEmax}}, 
          \code{\link{quadratic}}, \code{\link{exponential}},
          \code{\link{plot.Mods}}
          }

