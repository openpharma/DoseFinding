# Control options for pmvt and qmvt functions

Returns a list (an object of class "GenzBretz") with control parameters
for the `pmvt` and `qmvt` functions from the `mvtnorm` package. Note
that the DoseFinding package always uses "GenzBretz" algorithm. See the
mvtnorm documentation for more information.

## Usage

``` r
mvtnorm.control(maxpts = 30000, abseps = 0.001, releps = 0, interval = NULL)
```

## Arguments

- maxpts:

  Maximum number of function values as integer.

- abseps:

  Absolute error tolerance as double.

- releps:

  Relative error tolerance as double.

- interval:

  Interval to be searched, when the quantile is calculated.
