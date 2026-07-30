# Sample size calculations

The `sampSize` function implements a bisection search algorithm for
sample size calculation. The user can hand over a general target
function (via `targFunc`) that is then iterated so that a certain
`target` is achieved. The `sampSizeMCT` is a convenience wrapper of
`sampSize` for multiple contrast tests using the power as target
function.

## Usage

``` r
sampSize(
  upperN,
  lowerN = floor(upperN/2),
  targFunc,
  target,
  tol = 0.001,
  alRatio,
  Ntype = c("arm", "total"),
  verbose = FALSE
)

sampSizeMCT(
  upperN,
  lowerN = floor(upperN/2),
  ...,
  power,
  sumFct = mean,
  tol = 0.001,
  alRatio,
  Ntype = c("arm", "total"),
  verbose = FALSE
)

targN(
  upperN,
  lowerN,
  step,
  targFunc,
  alRatio,
  Ntype = c("arm", "total"),
  sumFct = c("min", "mean", "max")
)

powN(
  upperN,
  lowerN,
  step,
  ...,
  alRatio,
  Ntype = c("arm", "total"),
  sumFct = c("min", "mean", "max")
)

# S3 method for class 'targN'
plot(x, superpose = TRUE, line.at = NULL, xlab = NULL, ylab = NULL, ...)
```

## Arguments

- upperN, lowerN:

  Upper and lower bound for the target sample size. `lowerN` defaults to
  `floor(upperN/2)`.

- targFunc, target:

  The target function needs to take as an input the vector of sample
  sizes in the different dose groups. For `sampSize` it needs to return
  a univariate number. For function `targN` it should return a numerical
  vector.  
    
  Example: `targFunc` could be a function that calculates the power of a
  test, and `target` the desired target power value.  
  For function `sampSize` the bisection search iterates the sample size
  so that a specific target value is achieved (the implicit assumption
  is that targFunc is monotonically increasing in the sample size).  
    
  Function `targN` simply calculates `targFunc` for a given set of
  sample sizes.

- tol:

  A positive numeric value specifying the tolerance level for the
  bisection search algorithm. Bisection is stopped if the `targFunc`
  value is within `tol` of `target`.

- alRatio:

  Vector describing the relative patient allocations to the dose groups
  up to proportionality, e.g. `rep(1, length(doses))` corresponds to
  balanced allocations.

- Ntype:

  One of "arm" or "total". Determines, whether the sample size in the
  smallest arm or the total sample size is iterated in bisection search
  algorithm.

- verbose:

  Logical value indicating if a trace of the iteration progress of the
  bisection search algorithm should be displayed.

- ...:

  Arguments directly passed to the
  [`powMCT()`](https://openpharma.github.io/DoseFinding/reference/powMCT.md)
  function in the `sampSizeMCT` and `powN` function.

- power, sumFct:

  power is a numeric defining the desired summary power to achieve (in
  `sampSizeMCT`).

- step:

  Only needed for functions `targN` and `powN`. Stepsize for the sample
  size at which the target function is calculated. The steps are
  calculated via `seq(lowerN,upperN,by=step)`.

- x, superpose, line.at, xlab, ylab:

  arguments for the plot method of `targN` and `powN`, additional
  arguments are passed down to the low-level lattice plotting routines.

## Details

The `targN` functions calculates a general target function for different
given sample sizes. The `powN` function is a convenience wrapper of
`targN` for multiple contrast tests using the power as target function.

## References

Pinheiro J, Bornkamp B, Bretz F (2006). “Design and Analysis of Dose
Finding Studies Combining Multiple Comparisons and Modeling Procedures.”
*Journal of Biopharmaceutical Statistics*, **16**, 639-656.
[doi:10.1080/10543400600860428](https://doi.org/10.1080/10543400600860428)
.

Pinheiro J, Bornkamp B (2017). “Designing phase II dose-finding studies:
sample size, doses, and dose allocation weights.” In *Handbook of
Methods for Designing, Monitoring, and Analyzing Dose-Finding Trials*,
229–246. Chapman and Hall/CRC.

## See also

[`powMCT()`](https://openpharma.github.io/DoseFinding/reference/powMCT.md)

## Author

Jose Pinheiro, Bjoern Bornkamp

## Examples

``` r

## sampSize examples

## first define the target function
## first calculate the power to detect all of the models in the candidate set
fmodels <- Mods(linear = NULL, emax = c(25),
                logistic = c(50, 10.88111), exponential=c(85),
                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2),
                doses = c(0,10,25,50,100,150), placEff=0, maxEff=0.4,
                addArgs = list(scal=200))
## contrast matrix to use
contMat <- optContr(fmodels, w=1)
## this function calculates the power under each model and then returns
## the average power under all models
tFunc <- function(n){
  powVals <- powMCT(contMat, altModels=fmodels, n=n, sigma = 1,
                    alpha=0.05)
  mean(powVals)
}

## assume we want to achieve 80% average power over the selected shapes
## and want to use a balanced allocations
if (FALSE) { # \dontrun{
sSize <- sampSize(upperN = 80, targFunc = tFunc, target=0.8,
                  alRatio = rep(1,6), verbose = TRUE)
sSize


## Now the same using the convenience sampSizeMCT function
sampSizeMCT(upperN=80, contMat = contMat, sigma = 1, altModels=fmodels,
            power = 0.8, alRatio = rep(1, 6), alpha = 0.05)
## Alternatively one can also specify an S matrix
## covariance matrix in one observation (6 total observation result in a
## variance of 1 in each group)
S <- 6*diag(6)
## this uses df = Inf, hence a slightly smaller sample size results
sampSizeMCT(upperN=500, contMat = contMat, S=S, altModels=fmodels,
            power = 0.8, alRatio = rep(1, 6), alpha = 0.05, Ntype = "total")


## targN examples
## first calculate the power to detect all of the models in the candidate set
fmodels <- Mods(linear = NULL, emax = c(25),
                logistic = c(50, 10.88111), exponential=c(85),
                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2),
                doses = c(0,10,25,50,100,150), placEff=0, maxEff=0.4,
                addArgs = list(scal=200))
## corresponding contrast matrix
contMat <- optContr(fmodels, w=1)
## define target function
tFunc <- function(n){
  powMCT(contMat, altModels=fmodels, n=n, sigma = 1, alpha=0.05)
}
powVsN <- targN(upperN = 100, lowerN = 10, step = 10, tFunc,
                alRatio = rep(1, 6))
plot(powVsN)

## the same can be achieved using the convenience powN function
## without the need to specify a target function
powN(upperN = 100, lowerN=10, step = 10, contMat = contMat,
     sigma = 1, altModels = fmodels, alpha = 0.05, alRatio = rep(1, 6))
} # }
```
