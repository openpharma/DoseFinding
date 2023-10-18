\name{mvpostmix}
\alias{mvpostmix}
\title{Prior to posterior updating for a multivariate normal mixture}
\description{
  Calculate conjugate posterior mixture of multivariate normals with known covariance matrix 
}
\usage{
mvpostmix(priormix, mu_hat, S_hat)
}
\arguments{
  \item{priormix}{
    Prior multivariate normal mixture given as a list of length 3. The
    first list entry contains the mixture weights, the second component
    the mean vectors and the third component of the list the covariance matrices.
  }
  \item{mu_hat}{
    estimated mean response for each dose
  }
  \item{S_hat}{
    estimated covariance matrix
  }
}
\value{
  Returns a posterior multivariate normal mixture as a list of length 3,
  containing mixture weights, mean vectors and covariance matrices.
}
\references{
  Bernardo, J. M., and Smith, A. F. (1994). Bayesian theory. John Wiley & Sons.
}
\author{
  Marius Thomas
}
