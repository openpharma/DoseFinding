# Fit a dose-response model using Bayesian or bootstrap methods.

For `type = "Bayes"`, MCMC sampling from the posterior distribution of
the dose-response model is done. The function assumes a multivariate
normal distribution for `resp` with covariance matrix `S`, and this is
taken as likelihood function and combined with the prior distributions
specified in prior to form the posterior distribution.

## Usage

``` r
bFitMod(
  dose,
  resp,
  model,
  S,
  placAdj = FALSE,
  type = c("Bayes", "bootstrap"),
  start = NULL,
  prior = NULL,
  nSim = 1000,
  MCMCcontrol = list(),
  control = NULL,
  bnds,
  addArgs = NULL
)

# S3 method for class 'bFitMod'
predict(
  object,
  predType = c("full-model", "effect-curve"),
  summaryFct = function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)),
  doseSeq = NULL,
  lenSeq = 101,
  ...
)

# S3 method for class 'bFitMod'
plot(
  x,
  plotType = c("dr-curve", "effect-curve"),
  quant = c(0.025, 0.5, 0.975),
  plotData = c("means", "meansCI", "none"),
  level = 0.95,
  lenDose = 201,
  ...
)

# S3 method for class 'bFitMod'
coef(object, ...)
```

## Arguments

- dose:

  Numeric specifying the dose variable.

- resp:

  Numeric specifying the response estimate corresponding to the doses in
  `dose`

- model:

  Dose-response model to fit, possible models are "linlog", "linear",
  "quadratic", "emax", "exponential", "sigEmax", "betaMod" and
  "logistic", see
  [`drmodels()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md).

- S:

  Covariance matrix associated with the dose-response estimate specified
  via `resp`

- placAdj:

  Whether or not estimates in "placAdj" are placebo-adjusted (note that
  the linear in log and the logistic model cannot be fitted for
  placebo-adjusted data)

- type:

  Character with allowed values "Bayes" and "bootstrap", Determining
  whether samples are drawn from the posterior, or the bootstrap
  distribution.

- start:

  Optional starting values for the dose-response parameters in the MCMC
  algorithm.

- prior:

  List containing the information regarding the prior distributions for
  `type = "Bayes"`. The list needs to have as many entries as there are
  model parameters. The ordering of the list entries should be the same
  as in the arguments list of the model see (see
  [`drmodels()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)).
  For example for the Emax model the first entry determines the prior
  for e0, the second to eMax and the third to ed50.

  For each list entry the user has the choice to choose from 4 possible
  distributions:

  - `norm`: Vector of length 2 giving mean and standard deviation.

  - `t`: Vector of length 3 giving median, scale and degrees of freedom
    of the t-distribution.

  - `lnorm`: Vector of length 2 giving mean and standard deviation on
    log scale.

  - `beta`: Vector of length 4 giving lower and upper bound of the beta
    prior as well as the alpha and beta parameters of the beta
    distribution

- nSim:

  Desired number of samples to produce with the algorithm

- MCMCcontrol:

  List of control parameters for the MCMC algorithm

  - `thin` Thinning rate. Must be a positive integer.

  - `w` Numeric of same length as number of parameters in the model,
    specifies the width parameters of the slice sampler.

  - `adapt` Logical whether to adapt the `w` (width) parameter of the
    slice sampler in a short trial run. The widths are chosen as IQR/1.3
    of the trial run.

- control:

  Same as the control argument in
  [`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md).

- bnds:

  Bounds for non-linear parameters, in case `type = "bootstrap"`. If
  missing the the default bounds from
  [`defBnds()`](https://openpharma.github.io/DoseFinding/reference/defBnds.md)
  is used.

- addArgs:

  List containing two entries named "scal" and "off" for the "betaMod"
  and "linlog" model. When addArgs is NULL the following defaults are
  used `list(scal = 1.2*max(doses), off = 0.01*max(doses))`

- predType, summaryFct, doseSeq, lenSeq:

  Arguments for the predict method.

  `predType`: predType determines whether predictions are returned for
  the dose-response curve or the effect curve (difference to placebo).

  `summaryFct`: If equal to NULL predictions are calculated for each
  sampled parameter value. Otherwise a summary function is applied to
  the dose-response predictions for each parameter value. The default is
  to calculate 0.025, 0.25, 0.5, 0.75, 0.975 quantiles of the
  predictions for each dose.

  `doseSeq`: Where to calculate predictions. If not specified
  predictions are calculated on a grid of length `lenSeq` between
  minimum and maximum dose.

  `lenSeq`: Length of the default grid where to calculate predictions.

- ...:

  Additional arguments are ignored.

- x, object:

  A bFitMod object

- plotType, quant, plotData, level, lenDose:

  Arguments for plot method.

  `plotType`: Determining whether the dose-response curve or the effect
  curve should be plotted.

  `quant`: Vector of quantiles to display in plot

  `plotData`: Determines how the original data are plotted: Either as
  means or as means with CI or not. The level of the CI is determined by
  the argument `level`.

  `level`: Level for CI, when plotData is equal to `meansCI`.

  `lenDose`: Number of grid values to use for display.

## Value

An object of class bFitMod, which is a list containing the matrix of
posterior simulations plus some additional information on the fitted
model.

## Details

For `type = "bootstrap"`, a multivariate normal distribution for `resp`
with covariance matrix `S` is assumed, and a large number of samples is
drawn from this distribution. For each draw the fitMod function with
`type = "general"` is used to fit the draws from the multivariate normal
distribution.

Componentwise univariate slice samplers are implemented (see Neal 2003)
to sample from the posterior distribution.

## References

Neal RM (2003). “Slice sampling.” *The Annals of Statistics*, **31**(3),
705–767.
[doi:10.1214/aos/1056562461](https://doi.org/10.1214/aos/1056562461) .

## See also

[`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)

## Author

Bjoern Bornkamp

## Examples

``` r
data(biom)
anMod <- lm(resp~factor(dose)-1, data=biom)
drFit <- coef(anMod)
S <- vcov(anMod)
dose <- sort(unique(biom$dose))
## define prior list
## normal prior for E0 (mean=0 and sdev=10)
## normal prior for Emax (mean=0 and sdev=100)
## beta prior for ED50: bounds: [0,1.5] parameters shape1=0.45, shape2=1.7
prior <- list(norm = c(0, 10), norm = c(0,100), beta=c(0,1.5,0.45,1.7))
## now fit an emax model
gsample <- bFitMod(dose, drFit, S, model = "emax",
                   start = c(0, 1, 0.1), nSim = 1000, prior = prior)
## summary information
gsample
#> Dose Response Model
#> 
#> Model: emax 
#> 
#> Summary of posterior draws
#>       mean  sdev   2.5%   25%   50%   75% 97.5% n.eff
#> e0   0.371 0.138 0.0932 0.275 0.373 0.470 0.626   269
#> eMax 0.864 0.333 0.2923 0.647 0.815 1.057 1.588   157
#> ed50 0.407 0.343 0.0241 0.123 0.308 0.611 1.262   172
#> 
#> Fitted to:
#>       0    0.05     0.2     0.6       1 
#> 0.34491 0.45675 0.81032 0.93444 0.94871 
## samples are stored in
head(gsample$samples)
#>             e0      eMax       ed50
#> [1,] 0.2285613 0.6680033 0.03873810
#> [2,] 0.3556436 0.5905146 0.21497721
#> [3,] 0.3281728 0.4185104 0.09686954
#> [4,] 0.5147517 0.5065391 0.26245196
#> [5,] 0.5357584 0.3991353 0.34698651
#> [6,] 0.5983099 0.4332852 0.22866849
## predict 0.025, 0.25, 0.5, 0.75, 0.975 Quantile at 0, 0.5 and 1
predict(gsample, doseSeq = c(0, 0.5, 1))
#>               0       0.5         1
#> [1,] 0.09315894 0.6989780 0.7596054
#> [2,] 0.27503601 0.8082407 0.9003632
#> [3,] 0.37335209 0.8722832 0.9825789
#> [4,] 0.47044910 0.9334431 1.0687566
#> [5,] 0.62629836 1.0591235 1.2376306
## simple plot function
plot(gsample)


## now look at bootstrap distribution
gsample <- bFitMod(dose, drFit, S, model = "emax", type = "bootstrap",
                   nSim = 100, bnds = defBnds(1)$emax)
plot(gsample)


## now fit linear interpolation
prior <- list(norm = c(0,1000), norm = c(0,1000),
              norm = c(0,1000), norm = c(0,1000), norm = c(0,100))
gsample <- bFitMod(dose, drFit, S, model = "linInt",
                   start = rep(1,5), nSim = 1000, prior = prior)
gsample <- bFitMod(dose, drFit, S, model = "linInt", type = "bootstrap",
                   nSim = 100)

## data fitted on placebo adjusted scale
data(IBScovars)
anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
vCov <- vcov(anovaMod)[2:5,2:5]
dose <- sort(unique(IBScovars$dose))[-1]
prior <- list(norm = c(0,100), beta=c(0,6,0.45,1.7))
## Bayes fit
gsample <- bFitMod(dose, drFit, vCov, model = "emax", placAdj=TRUE,
                   start = c(1, 0.1), nSim = 1000, prior = prior)
## bootstrap fit
gsample <- bFitMod(dose, drFit, vCov, model = "emax", placAdj=TRUE,
                   type = "bootstrap", start = c(1, 0.1),
                   nSim = 100, prior = prior, bnds = c(0.01,6))
## calculate target dose estimate
TD(gsample, Delta = 0.2)
#>   [1] 0.328484679 0.021359123 1.100133333          NA 0.488505236 0.307754046
#>   [7] 0.633737653 0.232463613 0.007004876 0.326581628 0.010923631 3.777847585
#>  [13] 0.008681314 0.923299138 0.496244794 0.020377425 0.981195000 0.906321593
#>  [19] 0.897700910 1.134684881 0.040656828 1.796020526 0.385860070 0.007959519
#>  [25] 2.982013818 0.670871832 0.415647030 0.166747195 0.010979531          NA
#>  [31] 0.846949944 0.010298088 0.015255131 1.039726850 0.770168362 0.039217410
#>  [37] 0.043668547 2.186565977          NA 0.353336531 0.223994146 0.020984293
#>  [43] 0.350208918 0.012829304 2.724866852 0.006441993 0.260906993 0.025312635
#>  [49] 0.213119246 1.022828579 0.313769267 0.082124574 0.532987332 0.581898481
#>  [55] 2.131150244 1.255991534 1.284198185 0.598099581 0.564644464 0.014088899
#>  [61] 1.008487835          NA 0.677538538 0.023163561 0.946908952 1.754442525
#>  [67] 0.723084674 0.933466938 0.012367718 0.226667254 1.780935413 0.689908416
#>  [73] 1.760939465 0.012383204          NA 0.437408095 0.055769611          NA
#>  [79]          NA 0.024661051 0.388383728          NA 0.732919397 0.404642512
#>  [85] 0.359144120 0.348400198 0.421696183 8.634140945 0.402564631 0.967146474
#>  [91] 0.924591264 0.478732434 0.014010219 0.208402911 0.006788798 0.395789511
#>  [97] 0.009400354 0.015219517 0.015145893 1.078904099
## now fit linear interpolation
prior <- list(norm = c(0,1000), norm = c(0,1000), norm = c(0,1000), norm = c(0,100))
gsample <- bFitMod(dose, drFit, vCov, model = "linInt", placAdj=TRUE,
                   start = rep(1,4), nSim = 1000, prior = prior)
gsample <- bFitMod(dose, drFit, vCov, model = "linInt", type = "bootstrap",
                   placAdj = TRUE, nSim = 100)
```
