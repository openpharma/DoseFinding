\name{MCTtest}
\alias{MCTtest}
\title{
  Performs multiple contrast test
}
\description{
  This function performs a multiple contrast test. The contrasts are
  either directly specified in \samp{contMat} or optimal contrasts
  derived from the \samp{models} argument. The directionality of the
  data (i.e. whether an increase or decrease in the response variable is
  beneficial is inferred from the \samp{models} object, see
  \code{\link{Mods}}).

  For \samp{type = "normal"} an ANCOVA model based on a homoscedastic
  normality assumption (with additive covariates specified in
  \samp{addCovars}) is fitted.

  For \samp{type = "general"} it is assumed multivariate normally
  distributed estimates are specified in \samp{resp} with covariance
  given by \samp{S}, and the contrast test statistic is calculated based
  on this assumption. Degrees of freedom specified in \samp{df}.
}
\usage{
MCTtest(dose, resp, data = NULL, models, S = NULL, type = c("normal", "general"),
        addCovars = ~1, placAdj = FALSE, alpha = 0.025, df = NULL,
        critV = NULL, pVal = TRUE,
        alternative = c("one.sided", "two.sided"), na.action = na.fail,
        mvtcontrol = mvtnorm.control(), contMat = NULL)

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
    An object of class \samp{Mods}, see \code{\link{Mods}} for
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
  \item{alpha}{
    Significance level for the multiple contrast test
}
\item{df}{
  Specify the degrees of freedom to use in case \samp{type = "general"}.
  If this argument is missing \samp{df = Inf} is used (which corresponds
  to the multivariate normal distribution).
  For type = "normal" the degrees of freedom deduced from the AN(C)OVA fit are
  used and this argument is ignored.
}
\item{critV}{
  Supply a pre-calculated critical value. If this argument is NULL, no
  critical value will be calculated and the test decision is based on
  the p-values. If \samp{critV = TRUE} the critical value will be calculated. 
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
\item{contMat}{
  Contrast matrix to apply to the ANCOVA dose-response estimates. The
  contrasts need to be in the columns of the matrix (i.e. the column
  sums need to be 0).
}
}
\details{
  Integrals over the multivariate t and multivariate normal distribution
  are calculated using the \samp{mvtnorm} package.
}
\value{
  An object of class MCTtest, a list containing the output.
}
\references{
  Hothorn, T., Bretz, F., and Westfall, P. (2008). Simultaneous
  Inference in General Parametric Models, \emph{Biometrical Journal},
  \bold{50}, 346--363

  Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014)
  Model-based dose finding under model uncertainty using general
  parametric models, \emph{Statistics in Medicine}, \bold{33},
  1646--1661
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{powMCT}}, \code{\link{optContr}}
}
\examples{
## example without covariates
data(biom)
## define shapes for which to calculate optimal contrasts
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 1, 1, 1), doses = c(0, 0.05, 0.2, 0.6, 1))
m1 <- MCTtest(dose, resp, biom, models=modlist)
## now calculate critical value (but not p-values)
m2 <- MCTtest(dose, resp, biom, models=modlist, critV = TRUE, pVal = FALSE)
## now hand over critical value
m3 <- MCTtest(dose, resp, biom, models=modlist, critV = 2.24)

## example with covariates
data(IBScovars)
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 1, 1, 1), doses = c(0, 1, 2, 3, 4))
MCTtest(dose, resp, IBScovars, models = modlist, addCovars = ~gender)

## example using general approach (fitted on placebo-adjusted scale)
ancMod <- lm(resp~factor(dose)+gender, data=IBScovars)
## extract estimates and information to feed into MCTtest
drEst <- coef(ancMod)[2:5]
vc <- vcov(ancMod)[2:5, 2:5]
doses <- 1:4
MCTtest(doses, drEst, S = vc, models = modlist, placAdj = TRUE,
        type = "general", df = Inf)

## example with general alternatives handed over
data(biom)
## calculate contrast matrix for the step-contrasts
## represent them as linInt models
models <- Mods(linInt=rbind(c(1,1,1,1),
                            c(0,1,1,1),
                            c(0,0,1,1),
                            c(0,0,0,1)),
                doses=c(0,0.05,0.2,0.6,1))
plot(models)
## now calculate optimal contrasts for these means
## use weights from actual sample sizes
weights <- as.numeric(table(biom$dose))
contMat <- optContr(models, w = weights)
## plot contrasts
plot(contMat)
## perform multiple contrast test
MCTtest(dose, resp, data=biom, contMat = contMat)

## example for using the Dunnett contrasts
## Dunnett contrasts
doses <- sort(unique(biom$dose))
contMat <- rbind(-1, diag(4))
rownames(contMat) <- doses
colnames(contMat) <- paste("D", doses[-1], sep="")
MCTtest(dose, resp, data=biom, contMat = contMat)
}


