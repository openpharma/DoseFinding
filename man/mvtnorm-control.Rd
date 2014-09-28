\name{mvtnorm.control}
\alias{mvtnorm.control}
\title{ Control options for pmvt and qmvt functions }
\description{
  Returns a list (an object of class "GenzBretz") with control
  parameters for the \samp{pmvt} and \samp{qmvt} functions from the
  \samp{mvtnorm} package. Note that the DoseFinding package always uses
  "GenzBretz" algorithm. See the mvtnorm documentation for more
  information.
}
\usage{
mvtnorm.control(maxpts = 30000, abseps = 0.001, releps = 0,
           interval = NULL)
}
\arguments{
  \item{maxpts}{ Maximum number of function values as integer. }
  \item{abseps}{ Absolute error tolerance as double.  }
  \item{releps}{ Relative error tolerance as double.  }
  \item{interval}{
    Interval to be searched, when the quantile is calculated.
  }

}




