\name{MCPMod}
\alias{MCPMod}
\alias{predict.MCPMod}
\alias{plot.MCPMod}
\title{
  MCPMod - Multiple Comparisons and Modeling
}
\description{
  Tests for a dose-response effect using a model-based multiple contrast
  test (see \code{\link{MCTtest}}), selects one (or several) model(s)
  from the significant shapes, fits them using \code{\link{fitMod}}.
  For details on the method see Bretz et al. (2005).
}
\usage{
MCPMod(dose, resp, data, models, S = NULL, type = c("normal", "general"),
       addCovars = ~1, placAdj = FALSE, selModel = c("AIC", "maxT", "aveAIC"),
       alpha = 0.025, df = NULL, critV = NULL, doseType = c("TD", "ED"),
       Delta, p, pVal = TRUE, alternative = c("one.sided", "two.sided"),
       na.action = na.fail, mvtcontrol = mvtnorm.control(),
       bnds, control = NULL)

\method{predict}{MCPMod}(object,
        predType = c("full-model", "ls-means", "effect-curve"),
        newdata = NULL, doseSeq = NULL, se.fit = FALSE, ...)



\method{plot}{MCPMod}(x, CI = FALSE, level = 0.95,
       plotData = c("means", "meansCI", "raw", "none"),
       plotGrid = TRUE, colMn = 1, colFit = 1, ...)

}
\arguments{
  \item{dose, resp}{
    Either vectors of equal length specifying dose and response values,
    or names of variables in the data frame specified in \samp{data}.
  }
  \item{data}{
    Data frame containing the variables referenced in dose and resp if
    \samp{data} is not specified it is assumed that \samp{dose} and
    \samp{resp} are variables referenced from data (and no vectors)
  }
  \item{models}{
    An object of class \samp{"Mods"}, see \code{\link{Mods}} for
    details
  }
  \item{S}{
    The covariance matrix of \samp{resp} when \samp{type = "general"},
    see Description.
  }
  \item{type}{
    Determines whether inference is based on an ANCOVA model under a
    homoscedastic normality assumption (when \samp{type = "normal"}), or
    estimates at the doses and their covariance matrix and degrees of
    freedom are specified directly in \samp{resp}, \samp{S} and
    \samp{df}. See also \code{\link{fitMod}} and Pinheiro et al. (2013).
  }
  \item{addCovars}{
    Formula specifying additive linear covariates (for \samp{type = "normal"})
}
\item{placAdj}{
  Logical, if true, it is assumed that placebo-adjusted estimates are
  specified in \samp{resp} (only possible for \samp{type = "general"}).
}
\item{selModel}{
  Optional character vector specifying the model selection criterion for
  dose estimation.  Possible values are 
  \itemize{
    \item \code{AIC}: Selects model with smallest AIC (this is the default)
    \item \code{maxT}: Selects the model corresponding to the largest t-statistic.
    \item \code{aveAIC}: Uses a weighted average of the models corresponding to the
    significant contrasts.  The model weights are chosen by the formula:
    \eqn{w_i = \exp(-0.5AIC_i)/\sum_i(\exp(-0.5AIC_i))}{w_i = exp(-0.5AIC_i)/sum(exp(-0.5AIC_i))}
    See Buckland et al. (1997) for details.
  }
  For \samp{type = "general"} the "gAIC" is used.

}
\item{alpha}{
    Significance level for the multiple contrast test
}
\item{df}{
  Specify the degrees of freedom to use in case \samp{type = "general"},
  for the call to \code{\link{MCTtest}} and \code{\link{fitMod}}.
  Infinite degrees of (\samp{df=Inf}) correspond to the multivariate
  normal distribution.
  For type = "normal" the degrees of freedom deduced from the AN(C)OVA fit are
  used and this argument is ignored.
}
\item{critV}{
  Supply a pre-calculated critical value. If this argument is NULL, no
  critical value will be calculated and the test decision is based on
  the p-values. If \samp{critV = TRUE} the critical value will be calculated. 
}
  \item{doseType, Delta, p}{
    \samp{doseType} determines the dose to estimate, ED or TD (see also
    \code{\link{Mods}}), and \samp{Delta} and \samp{p} need to be
    specified depending on whether TD or ED is to be estimated.
    See \code{\link{TD}} and \code{\link{ED}} for details.
  }
\item{pVal}{
    Logical determining, whether p-values should be calculated.
}
\item{alternative}{
  Character determining the alternative for the multiple contrast
  trend test. 
}
\item{na.action}{
  A function which indicates what should happen when
  the data contain NAs.
}
\item{mvtcontrol}{
  A list specifying additional control parameters for the \samp{qmvt}
  and \samp{pmvt} calls in the code, see also \code{\link{mvtnorm.control}}
  for details.
}
\item{bnds}{
  Bounds for non-linear parameters. This needs to be a list with list
  entries corresponding to the selected bounds. The names of the list
  entries need to correspond to the model names. The
  \code{\link{defBnds}} function provides the default selection.
}
\item{control}{
  Control list for the optimization.\cr
  A list with entries: "nlminbcontrol", "optimizetol" and "gridSize".
  
  The entry nlminbcontrol needs to be a list and is passed directly to
  control argument in the nlminb function, that is used internally for
  models with 2 nonlinear parameters (e.g. sigmoid Emax or beta model).
    
  The entry optimizetol is passed directly to the tol argument of the
  optimize function, which is used for models with 1 nonlinear
  parameters (e.g. Emax or exponential model).

  The entry gridSize needs to be a list with entries dim1 and dim2
  giving the size of the grid for the gridsearch in 1d or 2d models.
}
\item{object, x}{
  MCPMod object
}
\item{predType, newdata, doseSeq, se.fit, ...}{
  predType determines whether predictions are returned for the full
  model (including potential covariates), the ls-means (SAS type) or the
  effect curve (difference to placebo).

  newdata gives the covariates to use in producing the predictions (for
  \samp{predType = "full-model"}), if missing the covariates used for fitting
  are used.

  doseSeq dose-sequence on where to produce predictions (for
  \samp{predType = "effect-curve"} and \samp{predType = "ls-means"}). If
  missing the doses used for fitting are used.
  
  se.fit: logical determining, whether the standard error should be
  calculated.

  \ldots: Additional arguments, for plot.MCPMod these are passed to
  plot.DRMod.
}
\item{CI, level, plotData, plotGrid, colMn, colFit}{
  Arguments for plot method: \samp{CI} determines whether confidence
  intervals should be plotted. \samp{level} determines the level of the
  confidence intervals. \samp{plotData} determines how the data are plotted:
  Either as means or as means with CI, raw data or none. In case of
  \samp{type = "normal"} and covariates the ls-means are displayed, when
  \samp{type = "general"} the option "raw" is not available.
  \samp{colMn} and \samp{colFit} determine the colors of fitted model
  and the raw means.
}

}
\value{
  An object of class \samp{MCPMod}, which contains the fitted
  \samp{MCTtest} object as well as the \samp{DRMod} objects and
  additional information (model selection criteria, dose estimates,
  selected models).
}
\references{
Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining
multiple comparisons and modeling techniques in dose-response
studies, \emph{Biometrics}, \bold{61}, 738--748

Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies
combining multiple comparisons and modeling procedures, \emph{Journal of Biopharmaceutical
Statistics}, \bold{16}, 639--656

Pinheiro, J. C., Bretz, F., and Branson, M. (2006). Analysis of dose-response studies - modeling 
approaches, \emph{in} N. Ting (ed.). \emph{Dose Finding in Drug Development}, Springer, New York,
pp. 146--171

Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014)
Model-based dose finding under model uncertainty using general
parametric models, \emph{Statistics in Medicine}, \bold{33},
1646--1661

Schorning, K., Bornkamp, B., Bretz, F., & Dette, H. (2016). Model
selection versus model averaging in dose finding
studies. \emph{Statistics in Medicine}, \bold{35}, 4021--4040

Bretz, F., Pinheiro, J. C., and Branson, M. (2004), On a hybrid method
in dose-finding studies, \emph{Methods of Information in Medicine},
\bold{43}, 457--460

Buckland, S. T., Burnham, K. P. and Augustin, N. H. (1997). Model selection an integral part
of inference, \emph{Biometrics}, \bold{53}, 603--618

Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley.
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{MCTtest}}, \code{\link{fitMod}}, \code{\link{drmodels}}
}
\examples{
data(biom)
## first define candidate model set (only need "standardized" models)
models <- Mods(linear = NULL, emax=c(0.05,0.2), linInt=c(1, 1, 1, 1),
               doses=c(0,0.05,0.2,0.6,1))
## perform MCPMod procedure
MM <- MCPMod(dose, resp, biom, models, Delta=0.5)
## a number of things can be done with an MCPMod object
MM # print method provides basic information
summary(MM) # more information
## predict all significant dose-response models
predict(MM, se.fit=TRUE, doseSeq=c(0,0.2,0.4, 0.9, 1),
        predType="ls-means")
## display all model functions 
plot(MM, plotData="meansCI", CI=TRUE)

## now perform model-averaging
MM2 <- MCPMod(dose, resp, biom, models, Delta=0.5, selModel = "aveAIC")
sq <- seq(0,1,length=11)
pred <- predict(MM, doseSeq=sq, predType="ls-means")
modWeights <- MM2$selMod
## model averaged predictions
pred <- do.call("cbind", pred)\%*\%modWeights
## model averaged dose-estimate
TDEst <- MM2$doseEst\%*\%modWeights

## now an example using a general fit and fitting based on placebo
## adjusted first-stage estimates
data(IBScovars)
## ANCOVA fit model including covariates
anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
vCov <- vcov(anovaMod)[2:5,2:5]
dose <- sort(unique(IBScovars$dose))[-1] # no estimate for placebo
## candidate models
models <- Mods(emax = c(0.5, 1), betaMod=c(1,1), doses=c(0,4))
## hand over placebo-adjusted estimates drFit to MCPMod
MM3 <- MCPMod(dose, drFit, S=vCov, models = models, type = "general",
              placAdj = TRUE, Delta=0.2)
plot(MM3, plotData="meansCI")

## The first example, but with critical value handed over
## this is useful, e.g. in simulation studies
MM4 <- MCPMod(dose, resp, biom, models, Delta=0.5, critV = 2.31)
}


