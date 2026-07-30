# Performs Bayesian multiple contrast test

This function performs a Bayesian multiple contrast test using normal
mixture priors for the response on each dose, as proposed in Fleischer
et al. (2022) . For a general description of the multiple contrast test
see
[`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md).

## Usage

``` r
bMCTtest(
  dose,
  resp,
  data = NULL,
  models,
  S = NULL,
  type = c("normal", "general"),
  prior,
  alpha = 0.025,
  na.action = na.fail,
  mvtcontrol = mvtnorm.control(),
  contMat = NULL,
  critV = NULL
)
```

## Arguments

- dose, resp:

  Either vectors of equal length specifying dose and response values, or
  names of variables in the data frame specified in `data`.

- data:

  Data frame containing the variables referenced in dose and resp if
  `data` is not specified it is assumed that `dose` and `resp` are
  variables referenced from data (and no vectors)

- models:

  An object of class `Mods`, see
  [`Mods()`](https://openpharma.github.io/DoseFinding/reference/Mods.md)
  for details

- S:

  The covariance matrix of `resp` when `type = "general"`, see
  Description.

- type:

  Determines whether inference is based on an ANCOVA model under a
  homoscedastic normality assumption (when `type = "normal"`), or
  estimates at the doses and their covariance matrix and degrees of
  freedom are specified directly in `resp`, `S` and `df`. See also
  [`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  and (Pinheiro et al. 2014) .

- prior:

  List of length equal to the number of doses with the prior for each
  arm. Each element needs to be of class `normMix` (See `RBesT` package
  documentation). It is assumed that the i-th component of the prior
  list corresponds to the i-th largest dose. For example the first entry
  in the list is the prior for the placebo group, the second entry the
  prior for the second lowest dose and so on. Internally the priors
  across the different arms are combined (densities multiplied) assuming
  independence. The resulting multivariate normal mixture prior will
  have as many components as the product of the number of components of
  the individual mixture priors. The posterior mixture is part of the
  result object under "posterior".

- alpha:

  Significance level for the frequentist multiple contrast test. If no
  critical values are supplied via `critV` this is used to derive
  critical values for Bayesian decision rule.

- na.action:

  A function which indicates what should happen when the data contain
  NAs.

- mvtcontrol:

  A list specifying additional control parameters for the `qmvt` and
  `pmvt` calls in the code, see also
  [`mvtnorm.control()`](https://openpharma.github.io/DoseFinding/reference/mvtnorm-control.md)
  for details.

- contMat:

  Contrast matrix to apply to the posterior dose-response estimates. The
  contrasts need to be in the columns of the matrix (i.e. the column
  sums need to be 0). If not specified optimal contrasts are calculated
  using
  [`optContr()`](https://openpharma.github.io/DoseFinding/reference/optContr.md).

- critV:

  Supply a critical value for the maximum posterior probability of the
  contrasts being greater than zero that needs to be surpassed to
  establish a non-flat dose-response. If this argument is NULL, this
  will be derived from critical values for frequentist MCP-Mod using the
  provided `alpha`.

## Value

An object of class bMCTtest, a list containing the output.

## Details

If `type = "normal"`, an ANCOVA model based on a homoscedastic normality
assumption is fitted and posteriors for dose-response and contrast
vectors are obtained assuming a known variance.

For `type = "general"` it is assumed multivariate normally distributed
estimates are specified in `resp` with covariance given by `S`, which
define the likelihood. Posteriors for dose-response and contrast vectors
are then obtained assuming a known covariance matrix S

The multiple contrast test decision is based on the maximum posterior
probability of a contrast being greater than zero. Thresholds for the
posterior probability can either be supplied or will be derived from
frequentist critical values. In the latter case the Bayesian test will
give approximately the same results as the frequentist multiple contrast
test if uninformative priors are used.

For the default calculation of optimal contrasts the prior information
is ignored (i.e. contrasts are calculated in the same way as in
[`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md)).
Fleischer et al. (2022) discuss using contrasts that take the prior
effective sample sizes into account, which can be slightly more
favourable for the Bayesian MCT test. Such alternative contrasts can be
directly handed over via the `contMat` argument.

For analysis with covariate adjustment, covariate-adjusted `resp` and
`S` can be supplied together with using `type = "general"`. See
[`vignette("binary_data")`](https://openpharma.github.io/DoseFinding/articles/binary_data.md)
vignette "Design and analysis template MCP-Mod for binary data" for an
example on how to obtain covariate adjusted estimates.

## References

Fleischer F, Bossert S, Deng Q, Loley C, Gierse J (2022). “Bayesian
MCPMod.” *Pharmaceutical Statistics*, **21**(3), 654–670.
[doi:10.1002/pst.2193](https://doi.org/10.1002/pst.2193) .  
  
Pinheiro J, Bornkamp B, Glimm E, Bretz F (2014). “Model-based dose
finding under model uncertainty using general parametric models.”
*Statistics in Medicine*, **33**, 1646-1661.
[doi:10.1002/sim.6052](https://doi.org/10.1002/sim.6052) .

## See also

[`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md),
[`optContr()`](https://openpharma.github.io/DoseFinding/reference/optContr.md)

## Author

Marius Thomas

## Examples

``` r


if (require("RBesT")) {

###############################
## Normal outcome
###############################

data(biom)
## define shapes for which to calculate optimal contrasts
doses <- c(0, 0.05, 0.2, 0.6, 1)
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 1, 1, 1), doses = doses)
## specify an informative prior for placebo, weakly informative for other arms
plc_prior <- mixnorm(inf = c(0.8, 0.4, 0.1), rob = c(0.2, 0.4, 10))
vague_prior <- mixnorm(c(1, 0, 10))
## i-th component of the prior list corresponds to the i-th largest dose
## (e.g. 1st component -> placebo prior; last component prior for top dose)
prior <- list(plc_prior, vague_prior, vague_prior, vague_prior, vague_prior)

m1 <- bMCTtest(dose, resp, biom, models=modlist, prior = prior)
## now supply a critical value (= threshold for maxmimum posterior probability)
m2 <- bMCTtest(dose, resp, biom, models=modlist, prior = prior, critV = 0.99)

####################################
## Binary outcome with covariates
####################################
if (FALSE) { # \dontrun{
## generate data
logit <- function(p) log(p / (1 - p))
inv_logit <- function(y) 1 / (1 + exp(-y))
doses <- c(0, 0.5, 1.5, 2.5, 4)

## set seed and ensure reproducibility across R versions
set.seed(1, kind = "Mersenne-Twister", sample.kind = "Rejection", normal.kind = "Inversion")
group_size <- 100
dose_vector <- rep(doses, each = group_size)
N <- length(dose_vector)
## generate covariates
x1 <- rnorm(N, 0, 1)
x2 <- factor(sample(c("A", "B"), N, replace = TRUE, prob = c(0.6, 0.4)))
## assume approximately logit(10%) placebo and logit(35%) asymptotic response with ED50=0.5
prob <- inv_logit(emax(dose_vector, -2.2, 1.6, 0.5) + 0.3 * x1 + 0.3 * (x2 == "B"))
dat <- data.frame(y = rbinom(N, 1, prob),
                  dose = dose_vector, x1 = x1, x2 = x2)

## specify an informative prior for placebo (on logit scale), weakly informative for other arms
plc_prior <- mixnorm(inf = c(0.8, -2, 0.5), rob = c(0.2, -2, 10))
vague_prior <- mixnorm(c(1, 0, 10))
prior <- list(plc_prior, vague_prior, vague_prior, vague_prior, vague_prior)

## candidate models
mods <- Mods(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1),
             placEff = logit(0.1), maxEff = logit(0.35)-logit(0.1),
             doses = doses)

fit_cov <- glm(y~factor(dose) + 0 + x1 + x2, data = dat, family = binomial)

covariate_adjusted_estimates <- function(mu_hat, S_hat, formula_rhs,
                                         doses, other_covariates, n_sim) {
  ## predict every patient under *every* dose
  oc_rep <- as.data.frame(lapply(other_covariates, function(col) rep(col, times = length(doses))))
  d_rep <- rep(doses, each = nrow(other_covariates))
  pdat <- cbind(oc_rep, dose = d_rep)
  X <- model.matrix(formula_rhs, pdat)
  ## average on probability scale then backtransform to logit scale
  mu_star <- logit(tapply(inv_logit(X %*% mu_hat), pdat$dose, mean))
  ## estimate covariance matrix of mu_star
  pred <- replicate(n_sim, logit(tapply(inv_logit(X %*% drop(mvtnorm::rmvnorm(1, mu_hat, S_hat))),
                                        pdat$dose, mean)))
  return(list(mu_star = as.numeric(mu_star), S_star = cov(t(pred))))
}

ca <- covariate_adjusted_estimates(coef(fit_cov), vcov(fit_cov), ~factor(dose)+0+x1+x2,
                                   doses, dat[, c("x1", "x2")], 1000)
bMCTtest(doses, ca$mu_star, S = ca$S_star, type = "general", models = mods, prior = prior)
} # }
################################################
## example with contrasts handed over
################################################

data(biom)
## define shapes for which to calculate optimal contrasts
doses <- c(0, 0.05, 0.2, 0.6, 1)
modlist <- Mods(emax = 0.05, linear = NULL, sigEmax = c(0.5, 5),
                linInt = c(0, 1, 1, 1), doses = doses)

## specify an informative prior for placebo, weakly informative for other arms
plc_prior <- mixnorm(inf = c(0.8, 0.4, 0.1), rob = c(0.2, 0.4, 10), sigma = 0.7)
vague_prior <- mixnorm(c(1, 0, 10), sigma = 0.7)
prior <- list(plc_prior, vague_prior, vague_prior, vague_prior, vague_prior)

## use prior effective sample sizes to calculate optimal contrasts
prior_ess <- unlist(lapply(prior, ess))
n_grp <- as.numeric(table(biom$dose))
weights <- n_grp + prior_ess
cmat <- optContr(modlist, w = weights)

bMCTtest(dose, resp, biom, models=modlist, prior = prior, contMat = cmat)
}
#> Loading required package: RBesT
#> This is RBesT version 1.10.0 (released 2026-07-01, git-sha 6d43b8e)
#> Using default prior reference scale 0.7
#> Using default prior reference scale 0.7
#> Using default prior reference scale 0.7
#> Using default prior reference scale 0.7
#> Using default prior reference scale 0.7
#> Bayesian MCP-Mod
#> 
#> Contrasts:
#>        emax linear sigEmax linInt
#> 0    -0.870 -0.681  -0.609 -0.762
#> 0.05  0.026 -0.191  -0.210 -0.263
#> 0.2   0.221 -0.060  -0.201  0.342
#> 0.6   0.302  0.291   0.400  0.342
#> 1     0.321  0.641   0.620  0.342
#> 
#> Posterior Mixture Weights:
#> Comp. 1 Comp. 2 
#>   0.995   0.005 
#> 
#> Bayesian t-statistics:
#>         Comp. 1 Comp. 2 posterior probability
#> linInt    4.167   3.372                     1
#> emax      4.111   2.996                     1
#> linear    3.709   3.198                     1
#> sigEmax   3.437   3.095                     1
#> 
#> Critical value (for maximum posterior probability): 0.987 (alpha = 0.025, one-sided) 
```
