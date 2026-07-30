# Neurodegenerative disease simulated longitudinal dose-finding data set

This simulated data set is motivated by a real Phase 2 clinical study of
a new drug for a neurodegenerative disease. The state of the disease is
measured through a functional scale, with smaller values corresponding
to more severe neurodeterioration. The goal of the drug is to reduce the
rate of disease progression, which is measured by the linear slope of
the functional scale over time.

## Usage

``` r
data(neurodeg)
```

## Format

A data frame with 100 observations on the following 2 variables.

- `resp`:

  a numeric vector containing the response values

- `dose`:

  a numeric vector containing the dose values

- `id`:

  Patient ID

- `time`:

  time of measurement

## Source

Pinheiro J, Bornkamp B, Glimm E, Bretz F (2014). “Model-based dose
finding under model uncertainty using general parametric models.”
*Statistics in Medicine*, **33**, 1646-1661.
[doi:10.1002/sim.6052](https://doi.org/10.1002/sim.6052) .

## Details

The trial design includes placebo and four doses: 1, 3, 10, and 30 mg,
with balanced allocation of 50 patients per arm. Patients are followed
up for one year, with measurements of the functional scale being taken
at baseline and then every three months.

The functional scale response is assumed to be normally distributed and,
based on historical data, it is believed that the longitudinal
progression of the functional scale over the one year of follow up can
be modeled a simple linear trend. See the example below on how to
analyse this type of data.

This data set was used in Pinheiro et al. (2014) to illustrate the
generalized MCPMod methodology.

## Examples

``` r

if (FALSE) { # \dontrun{
## reproduce analysis from Pinheiro et al. (2014)
data(neurodeg)
## first fit the linear mixed effect model
library(nlme)
fm <- lme(resp ~ as.factor(dose):time, neurodeg, ~time|id, method = "ML")
muH <- fixef(fm)[-1] # extract estimates
covH <- vcov(fm)[-1,-1]

## derive optimal contrasts for candidate shapes
doses <- c(0, 1, 3, 10, 30)
mod <- Mods(emax = 1.11, quadratic= -0.022, exponential = 8.867,
            linear = NULL, doses = doses) # 
contMat <- optContr(mod, S=covH) # calculate optimal contrasts
## multiple contrast test
MCTtest(doses, muH, S=covH, type = "general", critV = TRUE,
        contMat=contMat)
## fit the emax model
fitMod(doses, muH, S=covH, model="emax", type = "general",
       bnds=c(0.1, 10))


## alternatively one can also fit the model using nlme
nlme(resp ~ b0 + (e0 + eM * dose/(ed50 + dose))*time, neurodeg,
     fixed = b0 + e0 + eM + ed50 ~ 1, random = b0 + e0 ~ 1 | id,
     start = c(200, -4.6, 1.6, 3.2))
## both approaches lead to rather similar results
} # }
```
