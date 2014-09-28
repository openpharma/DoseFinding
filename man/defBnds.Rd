\name{defBnds}
\alias{defBnds}
\title{
  Calculates default bounds for non-linear parameters in dose-response
  models
}
\description{
  Calculates reasonable bounds for non-linear parameters for the
  built-in non-linear regression model based on the dose range under
  investigation.

  For the logistic model the first row corresponds to the ED50 parameter and the
  second row to the delta parameter. For the sigmoid Emax model the
  first row corresponds to the ED50 parameter and the second row to the
  h parameter, while for the beta model first and second row correspond
  to the delta1 and delta2 parameters. See \code{\link{logistic}},
  \code{\link{sigEmax}} and \code{\link{betaMod}} for details.
}
\usage{
   defBnds(mD, emax = c(0.001, 1.5)*mD,
           exponential = c(0.1, 2)*mD, 
           logistic = matrix(c(0.001, 0.01, 1.5, 1/2)*mD, 2),
           sigEmax = matrix(c(0.001*mD, 0.5, 1.5*mD, 10), 2),
           betaMod = matrix(c(0.05,0.05,4,4), 2))
}
\arguments{
  \item{mD}{
    Maximum dose in the study.
}
  \item{emax, exponential, logistic, sigEmax, betaMod}{
    values for the nonlinear parameters for these model-functions
}
}
\value{
  List containing bounds for the model parameters.
}
\author{
  Bjoern Bornkamp
}
\seealso{
\code{\link{fitMod}}
}
\examples{
  defBnds(mD = 1)
  defBnds(mD = 200)
}


