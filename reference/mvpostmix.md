# Prior to posterior updating for a multivariate normal mixture

Calculate conjugate posterior mixture of multivariate normals with known
covariance matrix

## Usage

``` r
mvpostmix(priormix, mu_hat, S_hat)
```

## Arguments

- priormix:

  Prior multivariate normal mixture given as a list of length 3. The
  first list entry contains the mixture weights, the second component
  the mean vectors and the third component of the list the covariance
  matrices.

- mu_hat:

  estimated mean response for each dose

- S_hat:

  estimated covariance matrix

## Value

Returns a posterior multivariate normal mixture as a list of length 3,
containing mixture weights, mean vectors and covariance matrices.

## References

Bernardo JM, Smith AF, Berliner M (1994). *Bayesian theory*, volume 586.
Wiley Online Library.

## Author

Marius Thomas
