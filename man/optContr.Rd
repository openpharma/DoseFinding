\name{optContr}
\alias{optContr}
\alias{plot.optContr}
\title{
  Calculate optimal contrasts
}
\description{
  This function calculates a contrast vectors that are optimal for
  detecting certain alternatives. The contrast is optimal in the sense
  of maximizing the non-centrality parameter of the underlying contrast
  test statistic:
  \deqn{\frac{c'\mu}{\sqrt{c'Sc}}}{c'mu/sqrt(c'Sc).}
  Here \eqn{\mu}{mu} is the mean vector under the alternative and \eqn{S}{S} the
  covariance matrix associated with the estimate of \eqn{\mu}{mu}.
  The optimal contrast is given by 
  \deqn{
    c^{opt} \propto S^{-1}\left(\mu - \frac{\mu^{\prime}S^{-1}
      1}{1^\prime S^{-1} 1}\right),}{c propto S^(-1) (mu - mu'S^(-1)1)/(1'S^(-1)1),}
  see Pinheiro et al. (2014).

  Note that the directionality (i.e. whether in "increase" in the
  response variable is beneficial or a "decrease", is inferred from the
  specified \samp{models} object, see \code{\link{Mods}} for details).

  Constrained contrasts (type = "constrained") add the additional
  constraint in the optimization that the sign of the contrast
  coefficient for control and active treatments need to be
  different. The quadratic programming algorithm from the quadprog
  package is used to calculate the contrasts.
}
\usage{
optContr(models, doses, w, S, placAdj = FALSE,
         type = c("unconstrained", "constrained"))

\method{plot}{optContr}(x, superpose = TRUE, xlab = "Dose",
     ylab = NULL, plotType = c("contrasts", "means"), ...)

}
\arguments{
  \item{models}{
    An object of class \samp{Mods} defining the dose-response shapes for
    which to calculate optimal contrasts. 
  }
  \item{doses}{
    Optional argument. If this argument is missing the doses attribute
    in the \samp{Mods} object specified in \samp{models} is used.
  }
\item{w, S}{
  Arguments determining the matrix S used in the formula for the optimal
  contrasts. Exactly one of \samp{w} and \samp{S} has
  to be specified. Note that \samp{w} and \samp{S} only have to be
  specified up to proportionality \cr
  \itemize{
    \item{w}{
      Vector specifying weights for the different doses,
      in the formula for calculation of the optimal contrasts. Specifying a weights
      vector is equivalent to specifying S=diag(1/w)
      (e.g. in a homoscedastic case with unequal sample sizes,
      \samp{w} should be proportional to the group sample sizes).}
    \item{S}
    {Directly specify a matrix proportional to the covariance
    matrix to use.}
  }
}
\item{placAdj}{
  Logical determining, whether the contrasts should be
  applied to placebo-adjusted estimates. If yes the returned
  coefficients are no longer contrasts (i.e. do not sum to 0). However,
  the result of multiplying of this "contrast" matrix with the placebo
  adjusted estimates, will give the same results as multiplying the
  original contrast matrix to the unadjusted estimates.
}
\item{type}{
  For \samp{type = "constrained"} the contrast coefficients of the zero
  dose group are constrained to be different from the coefficients of
  the active treatment groups. So that a weighted sum of the active
  treatments is compared against the zero dose group. For an increasing
  trend the coefficient of the zero dose group is negative and all other
  coefficients have to be positive (for a decreasing trend the other way
  round).
}
\item{x, superpose, xlab, ylab, plotType}{
  Arguments for the plot method for optContr objects. plotType
  determines, whether the contrasts or the underlying (standardized)
  mean matrix should be plotted.
}
\item{...}{
  Additional arguments for plot method
}
}
\value{
  Object of class \samp{optContr}. A list containing entries contMat and
  muMat (i.e. contrast, mean and correlation matrix).
}
\references{
  Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining multiple
  comparisons and modeling techniques in dose-response studies,
  \emph{Biometrics}, \bold{61}, 738--748

  Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014)
  Model-based dose finding under model uncertainty using general
  parametric models, \emph{Statistics in Medicine}, \bold{33},
  1646--1661
}
\author{
Bjoern Bornkamp
}

\seealso{
  \code{\link{MCTtest}}
}
\examples{
doses <- c(0,10,25,50,100,150)
models <- Mods(linear = NULL, emax = 25,
               logistic = c(50, 10.88111), exponential= 85,
               betaMod=rbind(c(0.33,2.31), c(1.39,1.39)),
               doses = doses, addArgs = list(scal = 200))
contMat <- optContr(models, w = rep(50,6))
plot(contMat)

## now we would like the "contrasts" for placebo adjusted estimates
dosPlac <- doses[-1]
## matrix proportional to cov-matrix of plac. adj. estimates for balanced data
S <- diag(5)+matrix(1, 5,5)
## note that we explicitly hand over the doses here
contMat0 <- optContr(models, doses=dosPlac, S = S, placAdj = TRUE)
## -> contMat0 is no longer a contrast matrix (columns do not sum to 0)
colSums(contMat0$contMat)
## calculate contrast matrix for unadjusted estimates from this matrix
## (should be same as above)
aux <- rbind(-colSums(contMat0$contMat), contMat0$contMat)
t(t(aux)/sqrt(colSums(aux^2))) ## compare to contMat$contMat

## now calculate constrained contrasts 
optContr(models, w = rep(50,6), type = "constrained")
optContr(models, doses=dosPlac, S = S, placAdj = TRUE,
         type = "constrained")
}


