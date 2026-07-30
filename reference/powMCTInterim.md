# Calculate Conditional or Predictive Power for Multiple Contrast Test

Calculates the predictive or conditional power for a multiple contrast
test based on interim data, e.g. for a futility interim analysis. This
function can also be applied to longitudinal endpoints, where at the
time of interim analysis incomplete data is available. For more details
see the vignette on longitudinal data analysis with MCP-Mod:
`vignette("Longitudinal Data MCP-Mod", package = "DoseFinding")`.

## Usage

``` r
powMCTInterim(
  contMat,
  mu_0t,
  S_0t,
  S_01,
  alpha = 0.025,
  type = c("predictive", "conditional"),
  mu_assumed = NULL,
  alternative = c("one.sided", "two.sided"),
  control = mvtnorm.control()
)
```

## Arguments

- contMat:

  Contrast matrix to use. The individual contrasts should be saved in
  the columns of the matrix

- mu_0t:

  The first stage estimates

- S_0t:

  The covariance matrix for the first stage estimates

- S_01:

  The covariance matrix anticipated for the estimates at study end

- alpha:

  Significance level to use

- type:

  Whether predictive power (for a flat prior) or conditional power
  should be calculated. For conditional power mu_assumed needs to be
  specified.

- mu_assumed:

  Mean vector to assume for the second stage (only used when type is
  `conditional`). If `NULL` (default), the first stage estimates `mu_0t`
  are used.

- alternative:

  Character determining the alternative for the multiple contrast trend
  test.

- control:

  A list specifying additional control parameters for the `pmvnorm`
  calls in the code, see also `mvtnorm.control` for details.

## Value

Numeric containing the calculated power values

## References

Bornkamp B, Zhou J, Xi D, Cao W (2024). “Futility analyses for the
MCP-Mod methodology based on longitudinal models.” 2406.19965,
<https://arxiv.org/abs/2406.19965>.

## See also

[`powMCT()`](https://openpharma.github.io/DoseFinding/reference/powMCT.md)
[`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md),
[`optContr()`](https://openpharma.github.io/DoseFinding/reference/optContr.md)

## Examples

``` r

# Setup the scenario.
doses <- c(0, 0.5, 1, 2, 4, 8)
mods <- Mods(
  emax = c(0.5, 1, 2, 4),
  sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4, 3)),
  quadratic = -0.1,
  doses = doses
)
w <- c(1, 0.5, 0.5, 0.5, 1, 1)
contMat <- optContr(models = mods, w = w)$contMat
sigma <- 0.3
n_final <- round(531 * w / sum(w))
n <- floor(n_final / 2)
S_0t <- diag(sigma^2 / n)
S_01 <- diag(sigma^2 / n_final)
## assumed interim estimate
mu_0t <- 0.05 * doses / (doses + 1) + rnorm(6, 0, 0.382 / sqrt(n))
## assumed mu (needed for conditional power)
mu_assumed <- 0.135 * doses / (doses + 1)

# Calculate predictive and conditional power.
powMCTInterim(
  contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t,
  type = "predictive"
)
#> [1] 0.0006866183
#> attr(,"error")
#> [1] 2.387509e-06
#> attr(,"msg")
#> [1] "Normal Completion"
powMCTInterim(
  contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t,
  type = "conditional", mu_assumed = mu_assumed
)
#> [1] 0.06110965
#> attr(,"error")
#> [1] 0.0008431539
#> attr(,"msg")
#> [1] "Normal Completion"
powMCTInterim(
  contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t,
  type = "predictive", alternative = "two.sided"
)
#> [1] 0.5590982
#> attr(,"error")
#> [1] 0.0005531084
#> attr(,"msg")
#> [1] "Normal Completion"
powMCTInterim(
  contMat = contMat, S_0t = S_0t, S_01 = S_01, mu_0t = mu_0t,
  type = "predictive", control = mvtnorm.control(maxpts = 1e5)
)
#> [1] 0.0007878812
#> attr(,"error")
#> [1] 0.000365923
#> attr(,"msg")
#> [1] "Normal Completion"
```
