# Migraine Dose Response data

Data set obtained from clinicaltrials.gov (NCT00712725). This was
randomized placebo controlled dose-response trial for treatment of acute
migraine. The primary endpoint was "pain freedom at 2 hours postdose" (a
binary measurement).

## Usage

``` r
data(migraine)
```

## Format

A data frame with 517 columns corresponding to the patients that
completed the trial

- `dose`:

  a numeric vector containing the dose values

- `painfree`:

  number of treatment responders

- `ntrt`:

  number of subject per treatment group

## Source

http://clinicaltrials.gov/ct2/show/results/NCT00712725
