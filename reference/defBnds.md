# Calculates default bounds for non-linear parameters in dose-response models

Calculates reasonable bounds for non-linear parameters for the built-in
non-linear regression model based on the dose range under investigation.

## Usage

``` r
defBnds(
  mD,
  emax = c(0.001, 1.5) * mD,
  exponential = c(0.1, 2) * mD,
  logistic = matrix(c(0.001, 0.01, 1.5, 1/2) * mD, 2),
  sigEmax = matrix(c(0.001 * mD, 0.5, 1.5 * mD, 10), 2),
  betaMod = matrix(c(0.05, 0.05, 4, 4), 2)
)
```

## Arguments

- mD:

  Maximum dose in the study.

- emax, exponential, logistic, sigEmax, betaMod:

  values for the nonlinear parameters for these model-functions

## Value

List containing bounds for the model parameters.

## Details

For the logistic model the first row corresponds to the ED50 parameter
and the second row to the delta parameter. For the sigmoid Emax model
the first row corresponds to the ED50 parameter and the second row to
the h parameter, while for the beta model first and second row
correspond to the delta1 and delta2 parameters. See
[`logistic()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md),
[`sigEmax()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
and
[`betaMod()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
for details.

## See also

[`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)

## Author

Bjoern Bornkamp

## Examples

``` r
  defBnds(mD = 1)
#> $emax
#> [1] 0.001 1.500
#> 
#> $logistic
#>       [,1] [,2]
#> [1,] 0.001  1.5
#> [2,] 0.010  0.5
#> 
#> $sigEmax
#>       [,1] [,2]
#> [1,] 0.001  1.5
#> [2,] 0.500 10.0
#> 
#> $exponential
#> [1] 0.1 2.0
#> 
#> $betaMod
#>      [,1] [,2]
#> [1,] 0.05    4
#> [2,] 0.05    4
#> 
  defBnds(mD = 200)
#> $emax
#> [1]   0.2 300.0
#> 
#> $logistic
#>      [,1] [,2]
#> [1,]  0.2  300
#> [2,]  2.0  100
#> 
#> $sigEmax
#>      [,1] [,2]
#> [1,]  0.2  300
#> [2,]  0.5   10
#> 
#> $exponential
#> [1]  20 400
#> 
#> $betaMod
#>      [,1] [,2]
#> [1,] 0.05    4
#> [2,] 0.05    4
#> 
```
