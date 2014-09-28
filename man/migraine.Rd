\name{migraine}
\alias{migraine}
\docType{data}
\title{ Migraine Dose Response data }
\description{
  Data set obtained from clinicaltrials.gov (NCT00712725).
  This was randomized placebo controlled dose-response trial for
  treatment of acute migraine. The primary endpoint was "pain freedom at
  2 hours postdose" (a binary measurement).
}
\usage{data(migraine)}
\format{
  A data frame with 517 columns corresponding to the patients that
  completed the trial
  \describe{
    \item{\code{dose}}{a numeric vector containing the dose values}
    \item{\code{painfree}}{number of treatment responders}
    \item{\code{ntrt}}{number of subject per treatment group}
  }
}
\source{
  http://clinicaltrials.gov/ct2/show/results/NCT00712725
}
\keyword{datasets}
