# Irritable Bowel Syndrome Dose Response data with covariates

A subset of the data used by (Biesheuvel and Hothorn, 2002). The data
are part of a dose ranging trial on a compound for the treatment of the
irritable bowel syndrome with four active treatment arms, corresponding
to doses 1,2,3,4 and placebo. Note that the original dose levels have
been blinded in this data set for confidentiality. The primary endpoint
was a baseline adjusted abdominal pain score with larger values
corresponding to a better treatment effect. In total 369 patients
completed the study, with nearly balanced allocation across the doses.

## Usage

``` r
data(IBScovars)
```

## Format

A data frame with 369 observations on the following 2 variables.

- `gender`:

  a factor specifying the gender

- `dose`:

  a numeric vector

- `resp`:

  a numeric vector

## Source

Biesheuvel E, Hothorn LA (2002). “Many-to-one comparisons in stratified
designs.” *Biometrical Journal*, **44**(1), 101–116.
[doi:10.1002/1521-4036(200201)44:13.0.CO;2-H](https://doi.org/10.1002/1521-4036%28200201%2944%3A13.0.CO%3B2-H)
.
