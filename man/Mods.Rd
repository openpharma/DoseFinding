\name{Mods}
\alias{Mods}
\alias{getResp}
\alias{plot.Mods}
\title{
  Define dose-response models
}
\description{
  The Mods functions allows to define a set of dose-response models.
  The function is used as input object for a number of other different
  functions.

  The dose-response models used in this package (see
  \code{\link{drmodels}} for details) are of form
  
  \deqn{f(d) = \theta_0+\theta_1 f^0(d,\theta_2)}{f(d) = theta0+theta1 f0(d,theta2)}

  where the parameter \eqn{\theta_2}{theta2} is the only non-linear
  parameter and can be one- or two-dimensional, depending on the used
  model.
  
  One needs to hand over the effect at placebo and the maximum effect in
  the dose range, from which \eqn{\theta_0,\theta_1}{theta0,theta1} are
  then back-calculated, the output object is of class
  \samp{"Mods"}. This object can form the input for other functions
  to extract the mean response (\samp{getResp}) or target doses
  (\code{\link{TD}} and \code{\link{ED}}) corresponding to the models. It
  is also needed as input to the functions \code{\link{powMCT}},
  \code{\link{optDesign}}
  
  Some models, for example the beta model (\samp{scal}) and the linlog
  model (\samp{off}) have parameters that are not estimated by the code,
  they need to be specified via the \samp{addArgs} argument.

  NOTE: If a decreasing effect is beneficial for the considered response
  variable it needs to specified here, either by using \samp{direction =
  "decreasing"} or by specifying a negative "maxEff" argument.
}
\usage{
Mods(..., doses, placEff = 0, maxEff,
     direction = c("increasing", "decreasing"), addArgs=NULL,
     fullMod = FALSE)

getResp(fmodels, doses)

\method{plot}{Mods}(x, nPoints = 200, superpose = FALSE, 
     xlab = "Dose", ylab = "Model means", modNams = NULL, 
     plotTD = FALSE, Delta, ...)

}
\arguments{
  \item{...}{
    In function Mods:\cr
    Dose-response model names with parameter values specifying the guesstimates for
    the \eqn{\theta_2}{theta2} parameters. See \code{\link{drmodels}} for a complete
    list of dose-response models implemented. See below for an example
    specification.\cr
    \cr
    In function plot.Mods:\cr
    Additional arguments to the \samp{xyplot} call.
  }
  \item{doses}{
    Dose levels to be used, this needs to include placebo.
  }
  \item{addArgs}{
    List containing two entries named "scal" and "off" for the
    "betaMod" and "linlog". When addArgs is NULL the following
    defaults are used \samp{list(scal = 1.2*max(doses), off =
    0.01*max(doses), nodes = doses)}.
  }
  \item{fullMod}{
    Logical determining, whether the model parameters specified in the Mods
    function (via the ... argument) should be interpreted as
    standardized or the full model parameters.
  }
  \item{placEff, maxEff}{
    Specify used placebo effect and the maximum effect over placebo.
    Either a numeric vector of the same size as the number of candidate
    models or of length one.\cr When these parameters are not specified
    \samp{placEff = 0} is assumed, for \samp{maxEff = 1} is assumed, if
    \samp{direction = "increasing"} and \samp{maxEff = -1} is assumed, for
    \samp{direction = "decreasing"}.
  }
  \item{direction}{
    Character determining whether the beneficial direction is \samp{increasing}
    or \samp{decreasing} with increasing dose levels. This argument is
    ignored if \samp{maxEff} is specified.
  }
  \item{fmodels}{
    An object of class Mods
  }
  \item{Delta}{
    Delta: The target effect size use for the target dose (TD) (Delta should be
    > 0).
  }
  \item{x}{
    Object of class Mods with type Mods
  }
  \item{nPoints}{ Number of points for plotting }
  \item{superpose}{ Logical determining, whether model plots should be superposed }
  \item{xlab, ylab}{ Label for y-axis and x-axis. }
  \item{modNams}{
    When \samp{modNams == NULL}, the names for the panels are determined
    by the underlying model functions, otherwise the contents of
    \samp{modNams} are used.
  }
  \item{plotTD}{
    \samp{plotTD} is a logical determining, whether the TD should be
    plotted. \samp{Delta} is the target effect to estimate for the TD.
  }
}
\value{
  Returns an object of class \samp{"Mods"}. The object contains the
  specified model parameter values and the derived linear parameters (based on
  \samp{"placEff"} and \samp{"maxEff"}) in a list.
}
\references{
  Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies
  combining multiple comparisons and modeling procedures, \emph{Journal of Biopharmaceutical
    Statistics}, \bold{16}, 639--656
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{Mods}}, \code{\link{drmodels}}, \code{\link{optDesign}}, \code{\link{powMCT}}
}
\examples{
## Example on how to specify candidate models

## Suppose one would like to use the following models with the specified
## guesstimates for theta2, in a situation where the doses to be used are
## 0, 0.05, 0.2, 0.6, 1

## Model            guesstimate(s) for theta2 parameter(s) (name)
## linear           -
## linear in log    -
## Emax             0.05 (ED50)
## Emax             0.3 (ED50)
## exponential      0.7 (delta)
## quadratic       -0.85 (delta)
## logistic         0.4  0.09 (ED50, delta)
## logistic         0.3  0.1 (ED50, delta)
## betaMod          0.3  1.3 (delta1, delta2)
## sigmoid Emax     0.5  2 (ED50, h)
## linInt           0.5 0.75 1 1 (perc of max-effect at doses)
## linInt           0.5 1 0.7 0.5 (perc of max-effect at doses)

## for the linInt model one specifies the effect over placebo for
## each active dose.
## The fixed "scal" parameter of the betaMod is set to 1.2
## The fixed "off"  parameter of the linlog is set to 0.1
## These (standardized) candidate models can be specified as follows

models <- Mods(linear = NULL, linlog = NULL, emax = c(0.05, 0.3),
               exponential = 0.7, quadratic = -0.85,
               logistic = rbind(c(0.4, 0.09), c(0.3, 0.1)),
               betaMod = c(0.3, 1.3), sigEmax = c(0.5, 2),
               linInt = rbind(c(0.5, 0.75, 1, 1), c(0.5, 1, 0.7, 0.5)),
               doses = c(0, 0.05, 0.2, 0.6, 1),
               addArgs = list(scal=1.2, off=0.1))
## "models" now contains the candidate model set, as placEff, maxEff and
## direction were not specified a placebo effect of 0 and an effect of 1
## is assumed

## display of specified candidate set
plot(models)

## example for creating a candidate set with decreasing response 
doses <- c(0, 10, 25, 50, 100, 150)
fmodels <- Mods(linear = NULL, emax = 25,
                   logistic = c(50, 10.88111), exponential = 85,
                   betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
                   linInt = rbind(c(0, 1, 1, 1, 1),
                                  c(0, 0, 1, 1, 0.8)), 
                   doses=doses, placEff = 0.5, maxEff = -0.4,
                   addArgs=list(scal=200))
plot(fmodels)
## some customizations (different model names, symbols, line-width)
plot(fmodels, lwd = 3, pch = 3, cex=1.2, col="red",
     modNams = paste("mod", 1:8, sep="-"))

## for a full-model object one can calculate the responses
## in a matrix
getResp(fmodels, doses=c(0, 20, 100, 150))

## calculate doses giving an improvement of 0.3 over placebo
TD(fmodels, Delta=0.3, direction = "decreasing")
## discrete version
TD(fmodels, Delta=0.3, TDtype = "discrete", doses=doses, direction = "decreasing")
## doses giving 50\% of the maximum effect
ED(fmodels, p=0.5)
ED(fmodels, p=0.5, EDtype = "discrete", doses=doses)

plot(fmodels, plotTD = TRUE, Delta = 0.3)

## example for specifying all model parameters (fullMod=TRUE)
fmods <- Mods(emax = c(0, 1, 0.1), linear = cbind(c(-0.4,0), c(0.2,0.1)),
              sigEmax = c(0, 1.1, 0.5, 3),
              doses = 0:4, fullMod = TRUE)
getResp(fmods, doses=seq(0,4,length=11))
## calculate doses giving an improvement of 0.3 over placebo
TD(fmods, Delta=0.3)
## discrete version
TD(fmods, Delta=0.3, TDtype = "discrete", doses=0:4)
## doses giving 50\% of the maximum effect
ED(fmods, p=0.5)
ED(fmods, p=0.5, EDtype = "discrete", doses=0:4)
plot(fmods)
}


