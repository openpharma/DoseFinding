\name{Target Doses}
\alias{TD}
\alias{ED}
\title{
  Calculate dose estimates for a fitted dose-response model (via
  \code{\link{fitMod}} or \code{\link{bFitMod}}) or a \code{\link{Mods}}
  object.
}
\usage{
TD(object, Delta, TDtype = c("continuous", "discrete"),
   direction = c("increasing", "decreasing"), doses)

ED(object, p, EDtype = c("continuous", "discrete"), doses)
}
\description{
  The TD (target dose) is defined as the dose that achieves a target
  effect of Delta over placebo (if there are multiple such doses, the
  smallest is chosen):

  \deqn{TD_\Delta = \min \{x|f(x) > f(0)+\Delta\}}{TD = min {x|f(x) > f(0)+Delta}}

  If a decreasing trend is beneficial the definition of the TD is

  \deqn{TD_\Delta = \min \{x|f(x) < f(0)-\Delta\}}{TD = min {x|f(x) < f(0)-Delta}}

  When \eqn{\Delta}{Delta} is the clinical relevance threshold, then the
  TD is similar to the usual definition of the minimum effective dose (MED).
  
  The ED (effective dose) is defined as the dose that achieves a certain
  percentage p of the full effect size (within the observed dose-range!)
  over placebo (if there are multiple such doses, the smallest is
  chosen).  \deqn{ED_p=\min\{x|f(x) > f(0) + p(f(dmax)-f(0))}{ EDp=min
  {x|f(x) > f(0) + p(f(dmax)-f(0))}}
  
  Note that this definition of the EDp is different from traditional
  definition based on the Emax model, where the EDp is defined relative
  to the \emph{asymptotic} maximum effect (rather than the maximum
  effect in the observed dose-range).
  
}
\arguments{
  \item{object}{
    An object of class c(Mods, fullMod), DRMod or bFitMod
  }
  \item{Delta, p}{
    Delta: The target effect size use for the target dose (TD) (Delta should be
    > 0).\cr
    p: The percentage of the dose to use for the effective dose.
  }
  \item{TDtype, EDtype}{
    character that determines, whether the dose should be treated as a
    continuous variable when calculating the TD/ED or whether the TD/ED
    should be calculated based on a grid of doses specified in
    \samp{doses}
  }
  \item{direction}{
    Direction to be used in defining the TD. This depends on whether an
    increasing or decreasing of the response variable is beneficial.
  }
  \item{doses}{
    Dose levels to be used, this needs to include placebo, \samp{TDtype}
    or \samp{EDtype} are equal to \samp{"discrete"}.
  }
}
\value{
  Returns the dose estimate
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{Mods}}, \code{\link{fitMod}}, \code{\link{bFitMod}}, \code{\link{drmodels}}
}
\examples{
## example for creating a "full-model" candidate set placebo response
## and maxEff already fixed in Mods call
doses <- c(0, 10, 25, 50, 100, 150)
fmodels <- Mods(linear = NULL, emax = 25,
                   logistic = c(50, 10.88111), exponential = 85,
                   betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
                   linInt = rbind(c(0, 1, 1, 1, 1),
                                  c(0, 0, 1, 1, 0.8)), 
                   doses=doses, placEff = 0, maxEff = 0.4,
                   addArgs=list(scal=200))

## calculate doses giving an improvement of 0.3 over placebo
TD(fmodels, Delta=0.3)
## discrete version
TD(fmodels, Delta=0.3, TDtype = "discrete", doses=doses)
## doses giving 50\% of the maximum effect
ED(fmodels, p=0.5)
ED(fmodels, p=0.5, EDtype = "discrete", doses=doses)

plot(fmodels, plotTD = TRUE, Delta = 0.3)
}
