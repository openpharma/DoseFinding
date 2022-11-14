\name{DoseFinding-package}
\alias{DoseFinding-package}
\alias{DoseFinding}
\docType{package}
\title{
  Design and Analysis of dose-finding studies
}
\description{
  The DoseFinding package provides functions for the design and analysis
  of dose-finding experiments (for example pharmaceutical Phase II
  clinical trials). It provides functions for: multiple contrast tests
  (\code{MCTtest}), fitting non-linear dose-response models
  (\code{fitMod}), a combination of testing and dose-response modelling
  (\code{MCPMod}), and calculating optimal designs (\code{optDesign}),
  both for normal and general response variable.
}
\details{
\tabular{ll}{
Package: \tab DoseFinding\cr
Type: \tab Package\cr
Version: \tab 1.0-3\cr
Date: \tab 2022-11-14\cr
License: \tab GPL-3
}
The main functions are:\cr
\bold{MCTtest}: Implements a multiple contrast tests\cr
\bold{powMCT}: Power calculations for multiple contrast tests\cr
\bold{fitMod}: Fits non-linear dose-response models\cr
\bold{optDesign}: Calculates optimal designs for dose-response models\cr
\bold{MCPMod}: Performs MCPMod methodology\cr
\bold{sampSize}: General function for sample size calculation\cr

}
\author{
Bjoern Bornkamp, Jose Pinheiro, Frank Bretz

Maintainer: Bjoern Bornkamp <bbnkmp@gmail.com>
}
\references{
  Bornkamp, B., Bretz, F., Dette, H. and Pinheiro,
  J. C. (2011). Response-Adaptive Dose-Finding under model uncertainty,
  \emph{Annals of Applied Statistics}, \bold{5}, 1611--1631
  
  Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An
  R Package for the Design and Analysis of Dose-Finding Studies,
  \emph{Journal of Statistical Software}, \bold{29}(7), 1--23
 
  Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining
  multiple comparisons and modeling techniques in dose-response
  studies, \emph{Biometrics}, \bold{61}, 738--748

  Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro,
  J. C. (2008). Optimal Designs for Dose Finding Studies, \emph{Journal
  of the American Statisical Association}, \bold{103}, 1225--1237
  
  O'Quigley, J., Iasonos, A. and Bornkamp, B. (2017) Handbook of
  methods for designing, monitoring, and analyzing dose-finding trials,
  CRC press, Part 3: Dose-Finding Studies in Phase II
  
  Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and
  analysis of dose finding studies combining multiple comparisons
  and modeling procedures, \emph{Journal of Biopharmaceutical
    Statistics}, \bold{16}, 639--656

  Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014)
  Model-based dose finding under model uncertainty using general
  parametric models, \emph{Statistics in Medicine}, \bold{33},
  1646--1661
  
  Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley
}
\keyword{ package }
\examples{
data(IBScovars)

## perform (model based) multiple contrast test
## define candidate dose-response shapes
models <- Mods(linear = NULL, emax = 0.2, quadratic = -0.17,
               doses = c(0, 1, 2, 3, 4))
## plot models
plot(models)
## perform multiple contrast test
test <- MCTtest(dose, resp, IBScovars, models=models,
                addCovars = ~ gender)

## fit non-linear emax dose-response model
fitemax <- fitMod(dose, resp, data=IBScovars, model="emax",
                  bnds = c(0.01,5))
## display fitted dose-effect curve
plot(fitemax, CI=TRUE, plotData="meansCI")

## Calculate optimal designs for target dose (TD) estimation
doses <- c(0, 10, 25, 50, 100, 150)
fmodels <- Mods(linear = NULL, emax = 25, exponential = 85,
                logistic = c(50, 10.8811),
                doses = doses, placEff=0, maxEff=0.4)
plot(fmodels, plotTD = TRUE, Delta = 0.2)
weights <- rep(1/4, 4)
desTD <- optDesign(fmodels, weights, Delta=0.2, designCrit="TD")
}
