# Glycopyrronium Bromide dose-response data

Data from a clinical study evaluating Efficacy and Safety of Four Doses
of Glycopyrronium Bromide in Patients With Stable Chronic Obstructive
Pulmonary Disease (COPD). This data set was obtained from
clinicaltrials.gov (NCT00501852). The study design was a 4 period
incomplete cross-over design. The primary endpoint is the trough forced
expiratory volume in 1 second (FEV1) following 7 days of Treatment.

## Usage

``` r
data(glycobrom)
```

## Format

A data frame with 5 summary estimates (one per dose). Variables: A data
frame with 5 summary estimates (one per dose). Variables:

- `dose`:

  a numeric vector containing the dose values

- `fev1`:

  a numeric vector containing the least square mean per dose

- `sdev`:

  a numeric vector containing the standard errors of the least square
  means per dose

- `n`:

  Number of participants analyzed per treatment group

## Source

http://clinicaltrials.gov/ct2/show/results/NCT00501852

## Details

The data given here are summary estimates (least-square means) for each
dose.

## Examples

``` r
 ## simulate a full data set with given means and sdv (here we ignore
  ## the original study was a cross-over design, and simulate a parallel
  ## group design)
  simData <- function(mn, sd, n, doses, fixed = TRUE){
    ## simulate data with means (mns) and standard deviations (sd), for
    ## fixed = TRUE, the data set will have observed means and standard
    ## deviations as given in mns and sd
    resp <- numeric(sum(n))
    uppind <- cumsum(n)
    lowind <- c(0,uppind)+1
    for(i in 1:length(n)){
      rv <- rnorm(n[i])
      if(fixed)
        rv <- scale(rv)
      resp[lowind[i]:uppind[i]] <- mn[i] + sd[i]*rv
    }
    data.frame(doses=rep(doses, n), resp=resp)
  }
  data(glycobrom)
  fullDat <- simData(glycobrom$fev1, glycobrom$sdev, glycobrom$n,
                     glycobrom$dose)
```
