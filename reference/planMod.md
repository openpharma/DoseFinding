# Evaluate performance metrics for fitting dose-response models

This function evaluates, the performance metrics for fitting
dose-response models (using asymptotic approximations or simulations).
Note that some metrics are available via the print method and others
only via the summary method applied to planMod objects. The implemented
metrics are

- Root of the mean-squared error to estimate the placebo-adjusted
  dose-response averaged over the used dose-levels, i.e. a rather
  discrete set (`dRMSE`). Available via the print method of planMod
  objects.

- Root of the mean-squared error to estimate the placebo-adjusted
  dose-response (`cRMSE`) averaged over fine (almost continuous) grid at
  101 equally spaced values between placebo and the maximum dose. NOTE:
  Available via the summary method applied to planMod objects.

- Ratio of the placebo-adjusted mean-squared error (at the observed
  doses) of model-based vs ANOVA approach (`Eff-vs-ANOVA`). This can be
  interpreted on the sample size scale. NOTE: Available via the summary
  method applied to planMod objects.

- Power that the (unadjusted) one-sided `1-alpha` confidence interval
  comparing the dose with maximum effect vs placebo is larger than
  `tau`. By default `alpha = 0.025` and `tau = 0` (`Pow(maxDose)`).
  Available via the print method of planMod objects.

- Probability that the EDp estimate is within the true \[EDpLB, EDpUB\]
  (by default `p=0.5`, `pLB=0.25` and `pUB=0.75`). This metric gives an
  idea on the ability to characterize the increasing part of the
  dose-response curve (`P(EDp)`). Available via the print method of
  planMod objects.

- Length of the quantile range for a target dose (TD or EDp). This is
  calculated by taking the difference of the dUB and dLB quantile of the
  empirical distribution of the dose estimates. (`lengthTDCI` and
  `lengthEDpCI`). It is NOT calculated by calculating confidence
  interval lengths in each simulated data-set and taking the mean. NOTE:
  Available via the summary method of planMod objects.

## Usage

``` r
planMod(
  model,
  altModels,
  n,
  sigma,
  S,
  doses,
  asyApprox = TRUE,
  simulation = FALSE,
  alpha = 0.025,
  tau = 0,
  p = 0.5,
  pLB = 0.25,
  pUB = 0.75,
  nSim = 100,
  cores = 1,
  showSimProgress = TRUE,
  bnds,
  addArgs = NULL
)

# S3 method for class 'planMod'
summary(
  object,
  digits = 3,
  len = 101,
  Delta = NULL,
  p = NULL,
  dLB = 0.05,
  dUB = 0.95,
  ...
)

# S3 method for class 'planMod'
plot(
  x,
  type = c("dose-response", "ED", "TD"),
  p,
  Delta,
  placAdj = FALSE,
  xlab = "Dose",
  ylab = "",
  ...
)
```

## Arguments

- model:

  Character vector determining the dose-response model(s) to be used for
  fitting the data. When more than one dose-response model is provided
  the best fitting model is chosen using the AIC. Built-in models are
  "linlog", "linear", "quadratic", "emax", "exponential", "sigEmax",
  "betaMod" and "logistic" (see
  [drmodels](https://openpharma.github.io/DoseFinding/reference/drmodels.md)).

- altModels:

  An object of class `Mods`, defining the true mean vectors under which
  operating characteristics should be calculated.

- n, sigma, S:

  Either a vector `n` and `sigma` or `S` need to be specified. When `n`
  and `sigma` are specified it is assumed computations are made for a
  normal homoscedastic ANOVA model with group sample sizes given by `n`
  and residual standard deviation `sigma`, i.e. the covariance matrix
  used for the estimates is thus `sigma^2*diag(1/n)` and the degrees of
  freedom are calculated as `sum(n)-nrow(contMat)`. When a single number
  is specified for `n` it is assumed this is the sample size per group
  and balanced allocations are used.  

  When `S` is specified this will be used as covariance matrix for the
  estimates.

- doses:

  Doses to use

- asyApprox, simulation:

  Logicals determining, whether asymptotic approximations or simulations
  should be calculated. If multiple models are specified in `model`
  asymptotic approximations are not available.

- alpha, tau:

  Significance level for the one-sided confidence interval for
  model-based contrast of best dose vs placebo. Tau is the threshold to
  compare the confidence interval limit to. CI(MaxDCont) gives the
  percentage that the bound of the confidence interval was larger than
  tau.

- p, pLB, pUB:

  p determines the type of EDp to estimate. pLB and pUB define the
  bounds for the EDp estimate. The performance metric Pr(Id-ED) gives
  the percentage that the estimated EDp was within the true EDpLB and
  EDpUB.

- nSim:

  Number of simulations

- cores:

  Number of cores to use for simulations. By default 1 cores is used,
  note that cores \> 1 will have no effect Windows, as the mclapply
  function is used internally.

- showSimProgress:

  In case of simulations show the progress using a progress-bar.

- bnds:

  Bounds for non-linear parameters. This needs to be a list with list
  entries corresponding to the selected bounds. The names of the list
  entries need to correspond to the model names. The
  [`defBnds()`](https://openpharma.github.io/DoseFinding/reference/defBnds.md)
  function provides the default selection.

- addArgs:

  See the corresponding argument in function
  [`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md).
  This argument is directly passed to fitMod.

- object, digits:

  object: A planMod object. digits: Digits in summary output

- len:

  Number of equally spaced points to determine the mean-squared error on
  a grid (cRMSE).

- Delta:

  Additional arguments determining what dose estimate to plot, when
  `type = "ED"` or `type = "TD"`

- dLB, dUB:

  Which quantiles to use for calculation of `lengthTDCI` and
  `lengthEDpCI`. By default dLB = 0.05 and dUB = 0.95, so that this
  corresponds to a 90% interval.

- ...:

  Additional arguments (currently ignored)

- x:

  An object of class planMod

- type:

  Type of plot to produce

- placAdj:

  When `type = "dose-response"`, this determines whether dose-response
  estimates are shown on placebo-adjusted or original scale

- xlab, ylab:

  Labels for the plot (ylab only applies for `type = "dose-response"`)

## Details

A plot method exists to summarize dose-response and dose estimations
graphically.

## References

Pinheiro J, Bornkamp B (2017). “Designing phase II dose-finding studies:
sample size, doses, and dose allocation weights.” In *Handbook of
Methods for Designing, Monitoring, and Analyzing Dose-Finding Trials*,
229–246. Chapman and Hall/CRC.

## See also

[`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)

## Author

Bjoern Bornkamp

## Examples

``` r
if (FALSE) { # \dontrun{
doses <- c(0,10,25,50,100,150)
fmodels <- Mods(linear = NULL, emax = 25,
                logistic = c(50, 10.88111), exponential= 85,
                betaMod=rbind(c(0.33,2.31),c(1.39,1.39)),
                doses = doses, addArgs=list(scal = 200),
                placEff = 0, maxEff = 0.4)
sigma <- 1
n <- rep(62, 6)*2

model <- "quadratic"
pObj <- planMod(model, fmodels, n, sigma, doses=doses,
               simulation = TRUE,
               alpha = 0.025, nSim = 200,
               p = 0.5, pLB = 0.25, pUB = 0.75)
print(pObj)
## to get additional metrics (e.g. Eff-vs-ANOVA, cRMSE, lengthTDCI, ...)
summary(pObj, p = 0.5, Delta = 0.3)
plot(pObj)
plot(pObj, type = "TD", Delta=0.3)
plot(pObj, type = "ED", p = 0.5)
} # }
```
