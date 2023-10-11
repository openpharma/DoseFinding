\name{critVal}
\alias{critVal}
\title{Calculate critical value for multiple contrast test}
\description{
  Calculation of the critical value for a maximum contrast test. This is
  based on the equicoordinate quantile function of the multivariate
  normal or t distribution as implemented in the \code{qmvt} function
  from the mvtnorm package.
}
\usage{
critVal(corMat, alpha = 0.025, df = NULL,
        alternative = c("one.sided", "two.sided"),
        control = mvtnorm.control())
}
\arguments{
\item{corMat}{
  Correlation matrix of contrasts
}
\item{alpha}{
  Significance level for the multiple contrast test
}
\item{df}{
  Specify the degrees of freedom to use, if this argument is missing
  \samp{df = Inf} is used (which corresponds to the multivariate
  normal distribution).
}
\item{alternative}{
  Character determining the alternative for the multiple contrast
  trend test. 
}
\item{control}{
  A list specifying additional control parameters for the \samp{qmvt}
  and \samp{pmvt} calls in the code, see also \code{\link{mvtnorm.control}}
  for details.
}
}
\seealso{
  \code{\link{powMCT}}, \code{\link{optContr}}, \code{\link{MCTtest}}
}
\examples{
R <- matrix(c(1,0.5,0.5,1), nrow=2)
critVal(R, alpha = 0.05, df = 1)
critVal(R, alpha = 0.05, df = 20)
critVal(R, alpha = 0.05, df = Inf)
}
\author{
  Bjoern Bornkamp
}
