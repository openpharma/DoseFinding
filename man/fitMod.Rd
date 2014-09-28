\name{fitMod}
\alias{fitMod}
\alias{coef.DRMod}
\alias{vcov.DRMod}
\alias{predict.DRMod}
\alias{plot.DRMod}
\alias{logLik.DRMod}
\alias{AIC.DRMod}
\alias{gAIC}
\alias{gAIC.DRMod}
\title{
  Fit non-linear dose-response model
}
\description{
  Fits a dose-response model. Built-in dose-response models are
  "linlog", "linear", "quadratic", "emax", "exponential", "sigEmax",
  "betaMod" and "logistic" (see \code{\link{drmodels}}).

  When \samp{type = "normal"} ordinary least squares is used and
  additional additive covariates can be specified in
  \samp{addCovars}. The underlying assumption is hence normally
  distributed data and homoscedastic variance.

  For \samp{type = "general"} a generalized least squares criterion is
  used \deqn{
  (f(dose,\theta)-resp)'S^{-1}(f(dose,\theta)-resp)}{(f(dose,theta)-resp)'S^{-1}(f(dose,theta)-resp)}
  and an inverse weighting matrix is specified in \samp{S}, \samp{type = "general"} 
  is primarily of interest, when fitting a model to
  AN(C)OVA type estimates obtained in a first stage fit, then
  \samp{resp} contains the estimates and \samp{S} is the estimated
  covariance matrix for the estimates in \samp{resp}. Statistical
  inference (e.g. confidence intervals) rely on asymptotic normality
  of the first stage estimates, which makes this method of interest only
  for sufficiently large sample size for the first stage fit. A modified
  model-selection criterion can be applied to these model fits (see also
  Pinheiro et al. 2013 for details).

  For details on the implemented numerical optimizer see the Details
  section below.
}
\usage{
fitMod(dose, resp, data = NULL, model, S = NULL, type = c("normal", "general"),
       addCovars = ~1, placAdj = FALSE, bnds, df = NULL,
       start = NULL, na.action = na.fail, control = NULL,
       addArgs = NULL)

\method{coef}{DRMod}(object, sep = FALSE, ...)

\method{predict}{DRMod}(object, predType = c("full-model", "ls-means", "effect-curve"),
        newdata = NULL, doseSeq = NULL, se.fit = FALSE, ...)

\method{vcov}{DRMod}(object, ...)

\method{plot}{DRMod}(x, CI = FALSE, level = 0.95,
     plotData = c("means", "meansCI", "raw", "none"),
     plotGrid = TRUE, colMn = 1, colFit = 1, ...)

\method{logLik}{DRMod}(object, ...)

\method{AIC}{DRMod}(object, ..., k = 2)

\method{gAIC}{DRMod}(object, ..., k = 2)

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
\item{model}{
  The dose-response model to be used for fitting the data. Built-in
  models are "linlog", "linear", "quadratic", "emax", "exponential",
  "sigEmax", "betaMod" and "logistic" (see \link{drmodels}).
}
  \item{S}{
    The inverse weighting matrix used in case, when \samp{type =
      "general"}, see Description. For later inference statements (vcov
      or predict methods) it is assumed this is the estimated covariance
      of the estimates in the first stage fit.
  }
\item{type}{
    Determines whether inference is based on an ANCOVA model under a
    homoscedastic normality assumption (when \samp{type = "normal"}), or
    estimates at the doses and their covariance matrix and degrees of
    freedom are specified directly in \samp{resp}, \samp{S} and
    \samp{df}. See also the Description above and Pinheiro et al. (2013).
  }
  \item{addCovars}{
    Formula specifying additional additive linear covariates (only for \samp{type = "normal"})
  }
\item{placAdj}{
  Logical, if true, it is assumed that placebo-adjusted estimates are
  specified in \samp{resp} (only possible for \samp{type = "general"}).
}
  \item{bnds}{
    Bounds for non-linear parameters. If missing the the default bounds
    from \code{\link{defBnds}} is used.
    
    When the dose-response model has only one non-linear parameter (for
    example Emax or exponential model), \samp{bnds} needs to be a vector
    containing upper and lower bound. For models with two non-linear
    parameters \samp{bnds} needs to be a matrix containing the bounds in
    the rows, see the Description section of \code{\link{defBnds}} for
    details on the formatting of the bounds for the individual models.
  }
  \item{df}{
    Degrees of freedom to use in case of \samp{type = "general"}. If
    this argument is missing \samp{df = Inf} is used. For \samp{type = "normal"} 
    this argument is ignored as the exact degrees of freedom
    can be deduced from the model.
  }
\item{start}{
  Vector of starting values for the nonlinear parameters (ignored
  for linear models). When equal to NULL, a grid optimization is performed and
  the best value is used as starting value for the local optimizer.
}
  \item{na.action}{
    A function which indicates what should happen when the data contain NAs.
}
\item{control}{
  A list with entries: "nlminbcontrol", "optimizetol" and "gridSize".
  
  The entry nlminbcontrol needs to be a list and it is passed directly
  to control argument in the nlminb function, that is used internally
  for models with 2 nonlinear parameters.
    
  The entry optimizetol is passed directly to the tol argument of the
  optimize function, which is used for models with 1 nonlinear
  parameters. 

  The entry gridSize needs to be a list with entries dim1 and dim2
  giving the size of the grid for the gridsearch in 1d or 2d models.
}
\item{addArgs}{
  List containing two entries named "scal" and "off" for the "betaMod"
  and "linlog" model. When addArgs is NULL the following defaults is
  used \samp{list(scal = 1.2*max(doses), off = 0.01*max(doses))}.
}
\item{object, x}{
  DRMod object
}
\item{sep}{
  Logical determining whether all coefficients should be returned in
  one numeric or separated in a list.
}
\item{predType, newdata, doseSeq, se.fit}{
  predType determines whether predictions are returned for the full
  model (including potential covariates), the ls-means (SAS type) or the
  effect curve (difference to placebo).

  newdata gives the covariates to use in producing the predictions (for
  predType = "full-model"), if missing the covariates used for fitting
  are used.

  doseSeq dose-sequence on where to produce predictions (for
  predType = "effect-curve" and predType = "ls-means"). If missing the
  doses used for fitting are used.
  
  se.fit: logical determining, whether the standard error should be
  calculated.
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
\item{k}{
  Penalty to use for model-selection criterion (AIC uses 2, BIC uses
  log(n)).
}
\item{...}{
  Additional arguments for plotting for the \samp{plot} method. For all other cases
  additional arguments are ignored.
}
}
\details{
  Details on numerical optimizer for model-fitting:\cr
  For linear models fitting is done using numerical linear algebra based
  on the QR decomposition.
  For nonlinear models numerical optimization is performed only in the
  nonlinear parameters in the model and optimizing over the linear
  parameters in each iteration (similar as the Golub-Pereyra
  implemented in \code{\link{nls}}). For models with 1 nonlinear
  parameter the \code{\link{optimize}} function is used for 2 nonlinear
  parameters the \code{\link{nlminb}} function is used. The starting
  value is generated using a grid-search (with the grid size specified
  via \samp{control$gridSize}), or can directly be handed over
  via \samp{start}.

  For details on the asymptotic approximation used for \samp{type = "normal"},
  see Seber and Wild (2003, chapter 5). For details on the asymptotic
  approximation used for \samp{type = "general"}, and the gAIC, see
  Pinheiro et al. (2013).
}
\value{
  An object of class DRMod. Essentially a list containing information
  about the fitted model coefficients, the residual sum of squares (or
  generalized residual sum of squares),
}
\references{
  Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014)
  Model-based dose finding under model uncertainty using general
  parametric models, \emph{Statistics in Medicine}, \bold{33},
  1646--1661

  Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley.
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{defBnds}}, \code{\link{drmodels}}
}
\examples{
## Fit the emax model to the IBScovars data set
data(IBScovars)
fitemax <- fitMod(dose, resp, data=IBScovars, model="emax",
                  bnds = c(0.01, 4))

## methods for DRMod objects
summary(fitemax)
## extracting coefficients
coef(fitemax)
## (asymptotic) covariance matrix of estimates
vcov(fitemax)
## predicting
newdat <- data.frame(dose = c(0,0.5,1), gender=factor(1))
predict(fitemax, newdata=newdat, predType = "full-model", se.fit = TRUE)
## plotting 
plot(fitemax, plotData = "meansCI", CI=TRUE)

## now include (additive) covariate gender
fitemax2 <- fitMod(dose, resp, data=IBScovars, model="emax",
                   addCovars = ~gender, bnds = c(0.01, 4))
vcov(fitemax2)
plot(fitemax2)
## fitted log-likelihood
logLik(fitemax2)
## extracting AIC (or BIC)
AIC(fitemax2)

## Illustrating the "general" approach for a binary regression
## produce first stage fit (using dose as factor)
data(migraine)
PFrate <- migraine$painfree/migraine$ntrt
doseVec <- migraine$dose
doseVecFac <- as.factor(migraine$dose)
## fit logistic regression with dose as factor
fitBin <- glm(PFrate~doseVecFac-1, family = binomial,
              weights = migraine$ntrt)
drEst <- coef(fitBin)
vCov <- vcov(fitBin)
## now fit an Emax model (on logit scale)
gfit <- fitMod(doseVec, drEst, S=vCov, model = "emax", bnds = c(0,100),
                type = "general")
## model fit on logit scale
plot(gfit, plotData = "meansCI", CI = TRUE)
## model on probability scale
logitPred <- predict(gfit, predType ="ls-means", doseSeq = 0:200,
                     se.fit=TRUE)
plot(0:200, 1/(1+exp(-logitPred$fit)), type = "l", ylim = c(0, 0.5),
     ylab = "Probability of being painfree", xlab = "Dose")
LB <- logitPred$fit-qnorm(0.975)*logitPred$se.fit
UB <- logitPred$fit+qnorm(0.975)*logitPred$se.fit
lines(0:200, 1/(1+exp(-LB)))
lines(0:200, 1/(1+exp(-UB)))


## now illustrate "general" approach for placebo-adjusted data (on
## IBScovars) note that the estimates are identical to fitemax2 above)
anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
vCov <- vcov(anovaMod)[2:5,2:5]
dose <- sort(unique(IBScovars$dose))[-1]
## now fit an emax model to these estimates
gfit2 <- fitMod(dose, drFit, S=vCov, model = "emax", type = "general",
               placAdj = TRUE, bnds = c(0.01, 2))
## some outputs
summary(gfit2)
coef(gfit2)
vcov(gfit2)
predict(gfit2, se.fit = TRUE, doseSeq = c(1,2,3,4), predType = "effect-curve")
plot(gfit2, CI=TRUE, plotData = "meansCI")
gAIC(gfit2)
}


