# Calculate critical value for multiple contrast test

Calculation of the critical value for a maximum contrast test. This is
based on the equicoordinate quantile function of the multivariate normal
or t distribution as implemented in the `qmvt` function from the mvtnorm
package.

## Usage

``` r
critVal(
  corMat,
  alpha = 0.025,
  df = NULL,
  alternative = c("one.sided", "two.sided"),
  control = mvtnorm.control()
)
```

## Arguments

- corMat:

  Correlation matrix of contrasts

- alpha:

  Significance level for the multiple contrast test

- df:

  Specify the degrees of freedom to use, if this argument is missing
  `df = Inf` is used (which corresponds to the multivariate normal
  distribution).

- alternative:

  Character determining the alternative for the multiple contrast trend
  test.

- control:

  A list specifying additional control parameters for the `qmvt` and
  `pmvt` calls in the code, see also
  [`mvtnorm.control()`](https://openpharma.github.io/DoseFinding/reference/mvtnorm-control.md)
  for details.

## See also

[`powMCT()`](https://openpharma.github.io/DoseFinding/reference/powMCT.md),
[`optContr()`](https://openpharma.github.io/DoseFinding/reference/optContr.md),
[`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md)

## Author

Bjoern Bornkamp

## Examples

``` r
R <- matrix(c(1,0.5,0.5,1), nrow=2)
critVal(R, alpha = 0.05, df = 1)
#> [1] 9.509978
critVal(R, alpha = 0.05, df = 20)
#> [1] 2.027525
critVal(R, alpha = 0.05, df = Inf)
#> [1] 1.916399
```
