# Calculate multiplicity adjusted p-values for multiple contrast test

Calculate multiplicity adjusted p-values for a maximum contrast test
corresponding to a set of contrasts and given a set of observed test
statistics. This function is exported as it may be a useful building
block and used in more complex testing situations that are not covered
by
[`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md).
Most users probably don't need to use this function.

## Usage

``` r
MCTpval(
  contMat,
  corMat,
  df,
  tStat,
  alternative = c("one.sided", "two.sided"),
  control = mvtnorm.control()
)
```

## Arguments

- contMat:

  Contrast matrix to use. The individual contrasts should be saved in
  the columns of the matrix

- corMat:

  Correlation matrix of contrasts

- df:

  Degrees of freedom to use for calculation.

- tStat:

  Vector of contrast test statistics

- alternative:

  Character determining the alternative for the multiple contrast trend
  test.

- control:

  A list specifying additional control parameters for the `qmvt` and
  `pmvt` calls in the code, see also
  [`mvtnorm.control()`](https://openpharma.github.io/DoseFinding/reference/mvtnorm-control.md)
  for details.

## Value

Numeric containing the calculated p-values.

## References

Pinheiro J, Bornkamp B, Bretz F (2006). “Design and Analysis of Dose
Finding Studies Combining Multiple Comparisons and Modeling Procedures.”
*Journal of Biopharmaceutical Statistics*, **16**, 639-656.
[doi:10.1080/10543400600860428](https://doi.org/10.1080/10543400600860428)
.

## See also

[`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md),
[`optContr()`](https://openpharma.github.io/DoseFinding/reference/optContr.md)

## Author

Bjoern Bornkamp

## Examples

``` r
data(biom)
## define shapes for which to calculate optimal contrasts
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 1, 1, 1), doses = c(0, 0.05, 0.2, 0.6, 1))
contMat <- optContr(modlist, w=1)$contMat
## calculate inputs needed for MCTpval
fit <- lm(resp~factor(dose)-1, data=biom)
est <- coef(fit)
vc <- vcov(fit)
ct <- as.vector(est %*% contMat)
covMat <- t(contMat) %*% vc %*% contMat
den <- sqrt(diag(covMat))
tStat <- ct/den
corMat <- cov2cor(t(contMat) %*% vc %*% contMat)
MCTpval(contMat, corMat, df=100-5, tStat)
#> [1] 0.001566885 0.004235790 0.007480728 0.001206962
## compare to
test <- MCTtest(dose, resp, biom, models=modlist)
attr(test$tStat, "pVal")
#> [1] 0.001572649 0.004888976 0.007579961 0.001275986
```
