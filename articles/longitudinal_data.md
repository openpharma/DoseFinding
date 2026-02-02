# Longitudinal Data MCP-Mod

In this vignette we illustrate how to use the `DoseFinding` package with
longitudinal, continuous observations by fitting a mixed model for
repeated measures (MMRM) and applying the generalized MCP-Mod
methodology to the resulting estimates. We also show how to apply
futility analyses during the trial, following Bornkamp et al.
([2024](#ref-bornkamp2024)).

``` r
library(DoseFinding)
```

## Data Simulation Code

In order to illustrate the methodology, we are going to consider data as
simulated in Section 4.1 of Bornkamp et al. ([2024](#ref-bornkamp2024)).
We include the corresponding functions here for completeness. This can
be helpful for users who want to simulate their own data during trial
design.

### Calculate Mean Response from Emax Model

The `calcMean` function calculates the mean response at each dose and
each time point under the Emax model with \\\text{E}\_0 = 0\\ and
\\\text{ED}\_{50} = 1\\ with the \\\text{E}\_{max}\\ parameter gradually
increasing so that a pre-specified maximum effect is achieved at the
last visit. The function takes three arguments: `dose` is a vector of
doses, `times` is a vector of time points, and `maxEff` is the maximum
effect at the final analysis (i.e., for the highest dose at the last
visit). It returns a matrix of mean responses, where the rows correspond
to doses and the columns correspond to time points.

``` r
calcMean <- function(dose, times, maxEff) {
  emaxLastVisit <- Mods(emax = 1, doses = dose, maxEff = maxEff)$emax["eMax"]
  maxTime <- max(times)
  outer(dose, times, function(doses, times) {
    Emax <- emaxLastVisit * (1 - exp(-0.5 * times)) / (1 - exp(-0.5 * maxTime))
    Emax * doses / (doses + 1)
  })
}

# Let's try it out:
calcMean(
  dose = c(0, 0.5, 1, 2, 4),
  times = 0:10,
  maxEff = 0.1
)
```

         [,1]       [,2]       [,3]       [,4]       [,5]       [,6]       [,7]
    [1,]    0 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    [2,]    0 0.01650577 0.02651703 0.03258916 0.03627210 0.03850591 0.03986079
    [3,]    0 0.02475866 0.03977554 0.04888374 0.05440814 0.05775886 0.05979118
    [4,]    0 0.03301154 0.05303405 0.06517832 0.07254419 0.07701182 0.07972157
    [5,]    0 0.03961385 0.06364086 0.07821399 0.08705303 0.09241418 0.09566588
               [,8]       [,9]      [,10]      [,11]
    [1,] 0.00000000 0.00000000 0.00000000 0.00000000
    [2,] 0.04068256 0.04118099 0.04148330 0.04166667
    [3,] 0.06102384 0.06177149 0.06222496 0.06250000
    [4,] 0.08136512 0.08236198 0.08296661 0.08333333
    [5,] 0.09763814 0.09883438 0.09955993 0.10000000

### Generate Normal Mean Vector and Covariance Matrix

The `parGen` function returns the mean and (compound-symmetry)
covariance matrix for the continuous outcomes, given the same input as
for `calcMean` above, plus the baseline mean `bslMean`, and the standard
deviation `errSD` and constant correlation `rho` of the residuals:

``` r
parGen <- function(times, doses, maxEff, bslMean = 0, errSD = 0.55, rho = 0.9) {
  times <- sort(times)
  ndim <- length(times)

  ## matrix of outcome means for each dose and visit
  MeanMat <- calcMean(doses, times, maxEff) + bslMean
  rownames(MeanMat) <- paste0("Dose", doses)
  colnames(MeanMat) <- paste0("Week", times)

  ## cov mat. for err
  sdV <- if (length(errSD) == 1) {
    rep(errSD, ndim)
  } else if (length(errSD) == ndim) {
    errSD
  } else {
    stop("Length of err sd should either be 1 (homogeneous) or equal to the number of visits in times")
  }

  ## CS covariance structure
  R <- diag(1, ndim)
  R[lower.tri(R)] <- R[upper.tri(R)] <- rho
  Sigm <- diag(sdV) %*% R %*% diag(sdV)

  list(visit = times, Sigm = Sigm, MeanMat = MeanMat, bslMean = bslMean)
}

# Let's try it out:
pars <- parGen(
  times = 0:10,
  doses = c(0, 0.5, 1, 2, 4),
  maxEff = 0.1,
  bslMean = 1.5
)
pars
```

    $visit
     [1]  0  1  2  3  4  5  6  7  8  9 10

    $Sigm
             [,1]    [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]    [,9]
     [1,] 0.30250 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225
     [2,] 0.27225 0.30250 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225
     [3,] 0.27225 0.27225 0.30250 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225
     [4,] 0.27225 0.27225 0.27225 0.30250 0.27225 0.27225 0.27225 0.27225 0.27225
     [5,] 0.27225 0.27225 0.27225 0.27225 0.30250 0.27225 0.27225 0.27225 0.27225
     [6,] 0.27225 0.27225 0.27225 0.27225 0.27225 0.30250 0.27225 0.27225 0.27225
     [7,] 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.30250 0.27225 0.27225
     [8,] 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.30250 0.27225
     [9,] 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.30250
    [10,] 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225
    [11,] 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225 0.27225
            [,10]   [,11]
     [1,] 0.27225 0.27225
     [2,] 0.27225 0.27225
     [3,] 0.27225 0.27225
     [4,] 0.27225 0.27225
     [5,] 0.27225 0.27225
     [6,] 0.27225 0.27225
     [7,] 0.27225 0.27225
     [8,] 0.27225 0.27225
     [9,] 0.27225 0.27225
    [10,] 0.30250 0.27225
    [11,] 0.27225 0.30250

    $MeanMat
            Week0    Week1    Week2    Week3    Week4    Week5    Week6    Week7
    Dose0     1.5 1.500000 1.500000 1.500000 1.500000 1.500000 1.500000 1.500000
    Dose0.5   1.5 1.516506 1.526517 1.532589 1.536272 1.538506 1.539861 1.540683
    Dose1     1.5 1.524759 1.539776 1.548884 1.554408 1.557759 1.559791 1.561024
    Dose2     1.5 1.533012 1.553034 1.565178 1.572544 1.577012 1.579722 1.581365
    Dose4     1.5 1.539614 1.563641 1.578214 1.587053 1.592414 1.595666 1.597638
               Week8    Week9   Week10
    Dose0   1.500000 1.500000 1.500000
    Dose0.5 1.541181 1.541483 1.541667
    Dose1   1.561771 1.562225 1.562500
    Dose2   1.582362 1.582967 1.583333
    Dose4   1.598834 1.599560 1.600000

    $bslMean
    [1] 1.5

Note that we assume here that the time is in weeks, but that could of
course easily be adapted as needed for your project.

### Simulate Longitudinal Outcomes

The `datGen` function simulates the longitudinal outcome data, given a
sample size of patients `N`, a list of parameters `pars` (as generated
by `parGen`), a vector of `doses` and corresponding allocation ratios
`alRatio`, the enrollment time of the last patient `LPFV.t`. The
distribution of the recruitment time can be specified by `RecrDist` and
`RandRecr`, with the following options:

- `RecrDist = "Quad"`: quadratic, \\a t^2\\, where \\a = N / LPFV.t^2\\.
  - `RandRecr = TRUE`: \\t \sim \sqrt{N \cdot \text{Unif}(0,1) / a}\\.
  - `RandRecr = FALSE`: \\t = \sqrt{i / a}\\, where \\i = 1, \dotsc,
    N\\.
- `Exp`: exponential, with rate \\\lambda = -\log(0.9 / N) / LPFV.t\\.
  - `RandRecr = TRUE`: \\t \sim \text{Exp}(\lambda)\\.
  - `RandRecr = FALSE`: \\t = F\_{\text{Exp}(\lambda)}^{-1}(i/N)\\ for
    \\i = 1, \dotsc, N-1\\; \\t = LPFV.t\\ for \\i = N\\.
- `Unif`: uniform, \\t \sim \text{Unif}(0, LPFV.t)\\.
  - `RandRecr = TRUE`: \\t \sim \text{Unif}(0, LPFV.t)\\.
  - `RandRecr = FALSE`: \\t = i \cdot LPFV.t / N\\ for \\i = 1, \dotsc,
    N\\.

``` r
datGen <- function(N = 300, pars, doses, alRatio, LPFV.t, RecrDist = "Quad", RandRecr = TRUE) {
  K <- length(doses)
  nt <- length(pars$visit)

  ## generate list of dose arms
  nblk <- ceiling(N / sum(alRatio)) # number of blocks
  d <- unlist(lapply(1:nblk, function(i) sample(rep(1:K, alRatio), sum(alRatio))))
  d <- d[1:N]

  err <- mvtnorm::rmvnorm(N, mean = rep(0, nt), sigma = pars$Sigm) # error term

  Y <- pars$MeanMat[d, ] + err # complete outcomes

  ## enrollment time (week)
  randT <- if (RecrDist == "Quad") {
    a <- N / LPFV.t^2
    if (RandRecr) {
      sqrt(N * runif(N) / a)
    } else {
      sqrt(c(1:N) / a)
    }
  } else if (RecrDist == "Exp") {
    lamb <- -log(0.9 / N) / LPFV.t
    if (RandRecr) {
      rexp(N, rate = lamb)
    } else {
      c(qexp((1:(N - 1)) / N, lamb), LPFV.t)
    }
  } else if (RecrDist == "Unif") {
    if (RandRecr) {
      runif(N, 0, LPFV.t)
    } else {
      (1:N) * LPFV.t / N
    }
  }

  trdat <- cbind(1:N, doses[d], randT, Y)
  colnames(trdat)[1:3] <- c("SUBJID", "Dose", "enroll.time")
  out <- tibble::as_tibble(trdat) |>
    tidyr::gather("Visit", "response", -c("SUBJID", "Dose", "enroll.time")) |>
    dplyr::arrange(SUBJID)
  out$week <- rep(pars$visit, N)
  out$cal.time <- out$enroll.time + out$week

  bsl <- out |>
    dplyr::filter(week <= 0) |>
    dplyr::select("SUBJID", "response") |>
    dplyr::rename(base = "response")
  out <- bsl |>
    dplyr::left_join(out, by = "SUBJID", multiple = "all") |>
    dplyr::mutate(chg = response - base) |>
    dplyr::select(- dplyr::all_of("response")) |>
    dplyr::rename(response = "chg")

  out |>
    dplyr::mutate(USUBJID = paste0("PAT", SUBJID)) |>
    dplyr::mutate_at("Visit", function(x) factor(x, levels = paste0("Week", pars$visit))) |>
    dplyr::mutate_at("Dose", function(x) factor(x, levels = doses)) |>
    dplyr::arrange(SUBJID, week)
}

# Let's try it out:
dat <- datGen(
  N = 300,
  pars = pars,
  doses = c(0, 0.5, 1, 2, 4),
  alRatio = c(1, 1, 1, 1, 1),
  LPFV.t = 10,
  RecrDist = "Quad",
  RandRecr = TRUE
)
head(dat, 15)
```

     [38;5;246m# A tibble: 15 × 9 [39m
       SUBJID  base Dose  enroll.time Visit   week cal.time response USUBJID
         [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<fct> [39m [23m        [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<fct> [39m [23m   [3m [38;5;246m<int> [39m [23m     [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<chr> [39m [23m  
     [38;5;250m 1 [39m      1 1.83  4            9.75 Week0      0     9.75   0      PAT1   
     [38;5;250m 2 [39m      1 1.83  4            9.75 Week1      1    10.8    0.198  PAT1   
     [38;5;250m 3 [39m      1 1.83  4            9.75 Week2      2    11.8   - [31m0 [39m [31m. [39m [31m103 [39m  PAT1   
     [38;5;250m 4 [39m      1 1.83  4            9.75 Week3      3    12.8    0.270  PAT1   
     [38;5;250m 5 [39m      1 1.83  4            9.75 Week4      4    13.8   - [31m0 [39m [31m. [39m [31m150 [39m  PAT1   
     [38;5;250m 6 [39m      1 1.83  4            9.75 Week5      5    14.8    0.183  PAT1   
     [38;5;250m 7 [39m      1 1.83  4            9.75 Week6      6    15.8    0.158  PAT1   
     [38;5;250m 8 [39m      1 1.83  4            9.75 Week7      7    16.8    0.454  PAT1   
     [38;5;250m 9 [39m      1 1.83  4            9.75 Week8      8    17.8   - [31m0 [39m [31m. [39m [31m0 [39m [31m28 [4m1 [24m [39m PAT1   
     [38;5;250m10 [39m      1 1.83  4            9.75 Week9      9    18.8    0.474  PAT1   
     [38;5;250m11 [39m      1 1.83  4            9.75 Week10    10    19.8    0.354  PAT1   
     [38;5;250m12 [39m      2 0.462 1            1.52 Week0      0     1.52   0      PAT2   
     [38;5;250m13 [39m      2 0.462 1            1.52 Week1      1     2.52   0.405  PAT2   
     [38;5;250m14 [39m      2 0.462 1            1.52 Week2      2     3.52   0.271  PAT2   
     [38;5;250m15 [39m      2 0.462 1            1.52 Week3      3     4.52   0.382  PAT2   

So we see that the function generates a `data.frame` with the subject
ID, baseline value, dose, enrollment time, visit and week plus calendar
time for each longitudinal observation.

### Cut Interim Data

The `intDat` function cuts the interim data when a certain proportion
`pct` of patients have their final outcome. It returns a list with the
interim analysis time `int.time`, the cut data `dat` and the
distribution of the last visit `dist` (i.e., the number of patients with
their last visit at each time point).

``` r
intDat <- function(data, pct = 0.5) {
  N <- max(data$SUBJID)
  tmax <- max(data$week)
  sorted_final_cal_time <- data |>
    dplyr::filter(week == tmax) |> 
    dplyr::pull(cal.time) |> 
    sort()
  N_required <- ceiling(pct * N)
  t.cut <- sorted_final_cal_time[N_required]
  out <- list(int.time = t.cut)
  out$dat <- subset(data, cal.time <= t.cut)
  out$dist <- out$dat |>
    dplyr::group_by(SUBJID) |>
    dplyr::summarize(lastvis = max(week)) |>
    dplyr::group_by(lastvis) |>
    dplyr::summarize(
        n = dplyr::n(), 
        pct = (100 * dplyr::n() / N)
    )
  out
}

# Let's try it out:
int_dat <- intDat(dat, pct = 0.5)
int_dat
```

    $int.time
    [1] 17.06187

    $dat
     [38;5;246m# A tibble: 2,998 × 9 [39m
       SUBJID  base Dose  enroll.time Visit  week cal.time response USUBJID
         [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<fct> [39m [23m        [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<fct> [39m [23m  [3m [38;5;246m<int> [39m [23m     [3m [38;5;246m<dbl> [39m [23m     [3m [38;5;246m<dbl> [39m [23m  [3m [38;5;246m<chr> [39m [23m  
     [38;5;250m 1 [39m      1 1.83  4            9.75 Week0     0     9.75    0     PAT1   
     [38;5;250m 2 [39m      1 1.83  4            9.75 Week1     1    10.8     0.198 PAT1   
     [38;5;250m 3 [39m      1 1.83  4            9.75 Week2     2    11.8    - [31m0 [39m [31m. [39m [31m103 [39m PAT1   
     [38;5;250m 4 [39m      1 1.83  4            9.75 Week3     3    12.8     0.270 PAT1   
     [38;5;250m 5 [39m      1 1.83  4            9.75 Week4     4    13.8    - [31m0 [39m [31m. [39m [31m150 [39m PAT1   
     [38;5;250m 6 [39m      1 1.83  4            9.75 Week5     5    14.8     0.183 PAT1   
     [38;5;250m 7 [39m      1 1.83  4            9.75 Week6     6    15.8     0.158 PAT1   
     [38;5;250m 8 [39m      1 1.83  4            9.75 Week7     7    16.8     0.454 PAT1   
     [38;5;250m 9 [39m      2 0.462 1            1.52 Week0     0     1.52    0     PAT2   
     [38;5;250m10 [39m      2 0.462 1            1.52 Week1     1     2.52    0.405 PAT2   
     [38;5;246m# ℹ 2,988 more rows [39m

    $dist
     [38;5;246m# A tibble: 4 × 3 [39m
      lastvis     n   pct
         [3m [38;5;246m<int> [39m [23m  [3m [38;5;246m<int> [39m [23m  [3m [38;5;246m<dbl> [39m [23m
     [38;5;250m1 [39m       7    48  16  
     [38;5;250m2 [39m       8    56  18.7
     [38;5;250m3 [39m       9    46  15.3
     [38;5;250m4 [39m      10   150  50  

## Longitudinal Data Analysis

Given a longitudinal dataset (at an interim analysis), we can perform
two analyses:

1.  Completers analysis: Subset the data to patients who have reached
    the last visit, and then apply a simple linear model.
2.  Repeated measures analysis: Fit a mixed model for repeated measures
    (MMRM) to the data.

In both cases, we then calculate the least square means
\\\widehat{\mu}\_{0t}\\ and the covariance matrix \\S\_{0t}\\ of the
least square means, as well as the residual standard deviation
\\\sigma\\ of the model.

Let’s code these in two corresponding functions:

### Completers Analysis

``` r
analyzeCompleters <- function(dat, eval.time = "Week12") {
  complDat <- dat |> 
    dplyr::filter(Visit == eval.time)
  fit <- lm(response ~ Dose + base, data = complDat)
  lsm <- emmeans::emmeans(fit, ~ Dose)
  summ <- subset(summary(lsm))
  rowIDs <- as.numeric(rownames(summ))
  list(
    mu0t = summ$emmean,
    sigma = sigma(fit),
    S0t = vcov(lsm)
  )
}

# Let's try it out:
resultCompleters <- analyzeCompleters(
  dat = int_dat$dat,
  eval.time = "Week10"
)
resultCompleters
```

    $mu0t
    [1] -0.06011853  0.04853842  0.06463132  0.11936622  0.14242645

    $sigma
    [1] 0.2499916

    $S0t
                    0           0.5             1             2             4
    0    1.838760e-03  3.272396e-06 -4.632441e-06  4.028387e-07  6.215931e-07
    0.5  3.272396e-06  2.248530e-03 -2.341013e-05  2.035753e-06  3.141233e-06
    1   -4.632441e-06 -2.341013e-05  2.049133e-03 -2.881835e-06 -4.446765e-06
    2    4.028387e-07  2.035753e-06 -2.881835e-06  2.500083e-03  3.866923e-07
    4    6.215931e-07  3.141233e-06 -4.446765e-06  3.866923e-07  1.953591e-03

### Repeated Measures Analysis

``` r
analyzeRepeated <- function(dat, eval.time = "Week12") {
  form <- "response ~  Visit + Dose + base + Dose*Visit + base*Visit"
  postBaselineDat <- dat |>
    dplyr::filter(week > 0) |> 
    droplevels()
  fit <- mmrm::mmrm(as.formula(paste0(form, " + us(Visit | USUBJID)")), data = postBaselineDat)
  lsm <- emmeans::emmeans(fit, ~ Dose + Visit)
  summ <- subset(summary(lsm), Visit == eval.time)
  rowIDs <- as.numeric(rownames(summ))
  list(
    mu0t = summ$emmean,
    sigma = sqrt(diag(mmrm::VarCorr(fit)))[eval.time],
    S0t = vcov(lsm)[rowIDs, rowIDs]
  )
}

# Let's try it out:
resultRepeated <- analyzeRepeated(
  dat = int_dat$dat,
  eval.time = "Week10"
)
resultRepeated
```

    $mu0t
    [1] -0.02818037  0.05291721  0.09861362  0.13468919  0.14456095

    $sigma
       Week10 
    0.2513171 

    $S0t
                    0 Week10    0.5 Week10      1 Week10      2 Week10
    0 Week10    1.430501e-03 -1.818752e-06  1.529028e-06 -5.639547e-07
    0.5 Week10 -1.818752e-06  1.626728e-03 -1.358336e-05  1.101611e-06
    1 Week10    1.529028e-06 -1.358336e-05  1.539021e-03 -8.100511e-08
    2 Week10   -5.639547e-07  1.101611e-06 -8.100511e-08  1.743068e-03
    4 Week10    1.596990e-07 -6.294172e-07  7.022021e-07 -1.325705e-07
                    4 Week10
    0 Week10    1.596990e-07
    0.5 Week10 -6.294172e-07
    1 Week10    7.022021e-07
    2 Week10   -1.325705e-07
    4 Week10    1.484844e-03

## Generalized MCP-Mod

### Final Analysis

At the final analysis, we would like to use the generalized MCP-Mod
methodology to test the hypothesis of a dose-response relationship. This
is easily done by providing the least square means via `resp` and their
covariance matrix via `S` to the `MCTtest` function, which will then
perform the MCP-Mod analysis. We just need to also specify the dose
levels and the models we want to test: Here we just use three models, it
could be more in practice.

``` r
doses <- c(0, 0.5, 1, 2, 4)
models <- Mods(
  emax = 2,
  sigEmax = c(0.5, 3),
  quadratic = -0.2,
  doses = doses
)

resultFinal <- analyzeRepeated(
  dat = dat,
  eval.time = "Week10"
)

testFinal <- MCTtest(
  dose = doses,
  resp = resultFinal$mu0t,
  S = resultFinal$S0t,
  models = models,
  alternative = "one.sided",
  type = "general",
  critV = TRUE,
  pVal = TRUE,
  alpha = 0.025
)

testFinal
```

    Multiple Contrast Test

    Contrasts:
          emax sigEmax quadratic
    0   -0.657  -0.788    -0.722
    0.5 -0.271  -0.203    -0.223
    1   -0.013   0.251     0.167
    2    0.309   0.363     0.611
    4    0.632   0.378     0.167

    Contrast Correlation:
               emax sigEmax quadratic
    emax      1.000   0.921     0.827
    sigEmax   0.921   1.000     0.941
    quadratic 0.827   0.941     1.000

    Multiple Contrast Test:
              t-Stat  adj-p
    emax       3.987 <0.001
    sigEmax    3.588 <0.001
    quadratic  3.574 <0.001

    Critical value: 2.18 (alpha = 0.025, one-sided)

``` r
testFinal$critV
```

    [1] 2.180088
    attr(,"Calc")
    [1] TRUE

``` r
testFinal$tStat
```

         emax   sigEmax quadratic 
     3.986915  3.588199  3.574414 
    attr(,"pVal")
    [1] 3.940677e-05 2.383234e-04 2.215375e-04

``` r
attr(testFinal$tStat, "pVal")
```

    [1] 3.940677e-05 2.383234e-04 2.215375e-04

### Futility Interim Analysis

Here we show how to apply the futility analysis during the trial, as
described in Bornkamp et al. ([2024](#ref-bornkamp2024)). The idea is
that based on the available interim data, we calculate the predictive or
conditional power to detect a significant dose-response relationship at
the final analysis. If this value is below a certain threshold, we may
want to stop the trial for futility. We can use the `powMCTInterim`
function to calculate these values. To do so, we need two additional
ingredients:

1.  The contrast matrix for the models we want to test, which is
    provided via the argument `contMat`.
2.  The covariance matrix at the final analysis, denoted by
    \\S\_{\[0,1\]}\\ in Bornkamp et al. ([2024](#ref-bornkamp2024)) and
    provided via the argument `S_01` here. Here we assume that this
    matrix is diagonal, with constant entries \\\sigma^2 / n\\, where
    \\\sigma\\ is the residual standard deviation of the model and \\n\\
    is the number of patients in the final analysis.
3.  For the conditional power, we also need to specify the assumed dose
    group means at the primary analysis time-point. This is provided via
    the argument `mu_assumed`. Typically, this is based on adding the
    observed placebo response rate to the originally planned treatment
    effect (before study start).

Let’s calculate these now for the example interim data from above:

``` r
w <- rep(1, length(doses))
contMat <- optContr(models = models, w = w)$contMat

n_final <- rep(60, length(doses))
S01 <- diag(resultRepeated$sigma^2 / n_final)

# Predictive power:
predPower <- powMCTInterim(
  contMat = contMat,
  mu_0t = resultRepeated$mu0t,
  S_0t = resultRepeated$S0t,
  S_01 = S01,
  alpha = 0.025,
  type = "predictive",
  alternative = "one.sided"
)
predPower
```

    [1] 0.9996943
    attr(,"error")
    [1] 1.504715e-07
    attr(,"msg")
    [1] "Normal Completion"

``` r
# Conditional power:
deltaAssumed <- pars$MeanMat[, "Week10"] - pars$bslMean 
deltaAssumed <- deltaAssumed - deltaAssumed[1]  # assumed treatment difference vs placebo
mu_assumed <- resultRepeated$mu0t[1] + deltaAssumed # to obtain group means add observed placebo response

condPower <- powMCTInterim(
  contMat = contMat,
  mu_0t = resultRepeated$mu0t,
  S_0t = resultRepeated$S0t,
  S_01 = S01,
  mu_assumed = mu_assumed,
  alpha = 0.025,
  type = "conditional",
  alternative = "one.sided"
)
condPower
```

    [1] 0.9978589
    attr(,"error")
    [1] 6.728031e-07
    attr(,"msg")
    [1] "Normal Completion"

## Simulation study

In practice, it will in most cases be important to run a simulation
study to evaluate the performance of the proposed combination of
generalized MCP-Mod with MMRM. In particular, in order to calibrate the
futility threshold for the interim analysis, it is important to
understand the power loss at the final analysis in relation to the
futility threshold. The following code may serve as the starting point
for such a simulation study, and can be adapted as needed. We include
here the argument documentation in `roxygen2` format for simplicity. The
code can also be used to reproduce the simulation study in Section 4.1
of Bornkamp et al. ([2024](#ref-bornkamp2024)).

``` r
##' @param times Vector assessment times
##' @param bslMean Mean at baseline
##' @param errSD Residual standard deviation (on absolute scale)
##' @param rho Correlation of measurements over time (within patient)
##' @param maxEff Maximum effect achieved at last time-point for the highest dose
##' @param RecrDist Recruitment distribution (see function datGen)
##' @param LPFV.t Time for the last patient first visit
##' @param models DoseFinding::Mods object
##' @param ia_timing Information time when IAs are performed (% of patients having last visit)
##' @param N Overall sample size
##' @param doses doses to use
##' @param alRatio allocation ratios to the different doses
##' @param eval.time Name of final
##' @param alpha Type 1 error for test
##' @param nSim Number of simulations
##' @param delta_assumed Vector of treatment effects assumed at the different doses
##' @return Data frame with all results
sim_one_trial <- function(times, bslMean, errSD, rho, maxEff,
                          RecrDist, LPFV.t, models, ia_timing = seq(0.3, 0.8, 0.05), N = 531,
                          doses = c(0, 0.5, 1, 2, 4, 8), alRatio = c(2, 1, 1, 1, 2, 2),
                          eval.time = "Week12", alpha = 0.025, nSim = 1000, delta_assumed) {
  contMat <- optContr(models = models, w = alRatio / sum(alRatio))$contMat
  pars <- parGen(times, doses,
    maxEff = maxEff, bslMean = bslMean, errSD = errSD,
    rho = rho
  )
  mydat <- datGen(N, pars, doses, alRatio, LPFV.t = LPFV.t, RecrDist = RecrDist, RandRecr = TRUE)

  ## fit final data
  fit.final <- analyzeCompleters(mydat)
  test <- MCTtest(
    dose = doses, resp = fit.final$mu0t, S = fit.final$S0t, models = models,
    alternative = "one.sided", type = "general", critV = TRUE, pVal = FALSE,
    alpha = alpha
  )
  maxT <- max(test$tStat)
  cVal <- test$critVal

  n_final <- mydat |>
    dplyr::filter(Visit == eval.time) |>
    dplyr::group_by(Dose) |>
    dplyr::count()

  IAtm <- pp_long <- cp_long <- cp_long2 <- pp_compl <- cp_compl <- cp_compl2 <- inf_long <- inf_compl <- numeric()
  IAdist <- matrix(nrow = length(times), ncol = length(ia_timing))

  ## fit interim data
  for (ia in seq_along(ia_timing)) {
    tp <- ia_timing[ia] # prop pts completes visit of eval.time
    cutDat <- intDat(mydat, tp)
    middat <- cutDat$dat
    IAtm[ia] <- cutDat$int.time
    IAdist[times %in% cutDat$dist$lastvis, ia] <- cutDat$dist$pct

    ## using all data from every patient (longitudinal MMRM)
    fit_int_long <- analyzeRepeated(middat, eval.time = eval.time)
    ## only using the completers by visit eval time (cross-sectional linear model)
    fit_int_compl <- analyzeCompleters(middat, eval.time = eval.time)
    ## The covariance matrix anticipated for the estimates at EoS
    S_end <- diag(fit_int_long$sigma^2 / n_final$n)

    # information fraction
    inf_long[ia] <- det(diag(fit_int_long$sigma^2 / n_final$n))^(1 / length(doses)) /
      det(fit_int_long$S0t)^(1 / length(doses))
    inf_compl[ia] <- det(diag(fit_int_compl$sigma^2 / n_final$n))^(1 / length(doses)) /
      det(fit_int_compl$S0t)^(1 / length(doses))

    mu_assumed <- fit_int_long$mu0t[1] + delta_assumed # for conditional power use planned treatment effect
    pp_long[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_long$S0t,
      S_01 = S_end, mu_0t = fit_int_long$mu0t, type = "predictive", alpha = alpha
    )
    cp_long[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_long$S0t,
      S_01 = S_end, mu_0t = fit_int_long$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = mu_assumed
    )
    cp_long2[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_long$S0t,
      S_01 = S_end, mu_0t = fit_int_long$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = fit_int_long$mu0t
    )
    pp_compl[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_compl$S0t,
      S_01 = S_end, mu_0t = fit_int_compl$mu0t, type = "predictive", alpha = alpha
    )
    cp_compl[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_compl$S0t,
      S_01 = S_end, mu_0t = fit_int_compl$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = mu_assumed
    )
    cp_compl2[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_compl$S0t,
      S_01 = S_end, mu_0t = fit_int_compl$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = fit_int_compl$mu0t
    )
  }
  data.frame(
    final_maxT = maxT, final_cVal = as.numeric(cVal),
    ia_timing, pp_long, cp_long, cp_long2, pp_compl, cp_compl, cp_compl2, inf_long, inf_compl
  )
}

## Function to run one scenario (arguments described in previous function)
run_scen <- function(n_sim, rho, maxEff,
                     LPFV.t, ia_timing, delta_assumed) {
  ## fixed parameters
  doses <- c(0, 0.5, 1, 2, 4, 8) # doses
  alRatio <- c(2, 1, 1, 1, 2, 2) # allocation ratio for doses
  alpha <- 0.025
  times <- c(0, 2, 4, 8, 12)
  bslMean <- 0
  errSD <- 0.56
  RecrDist <- "Quad"

  if (rho == 0.6 & LPFV.t == 50) N <- 820
  if (rho == 0.6 & LPFV.t == 100) N <- 825
  if (rho == 0.9 & LPFV.t == 50) N <- 248
  if (rho == 0.9 & LPFV.t == 100) N <- 236

  models <- Mods(
    emax = c(0.5, 1, 2, 4), sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4, 3)),
    quadratic = -0.1, doses = doses
  )

  lst <- vector("list", n_sim)
  for (i in 1:n_sim) {
    res <- try(sim_one_trial(
      times = times, bslMean = bslMean,
      errSD = errSD, rho = rho, maxEff = maxEff, RecrDist = RecrDist,
      LPFV.t = LPFV.t, models = models, ia_timing = ia_timing, N = N,
      doses = doses, alRatio = alRatio, delta_assumed = delta_assumed,
      alpha = alpha
    ), silent = FALSE)
    if (!is.character(res)) {
      lst[[i]] <- data.frame(sim = i, res)
    }
  }

  out <- do.call("rbind", lst)
  out$rho <- rho
  out$maxEff <- maxEff
  out$RecrDist <- RecrDist
  out$N <- N
  out$LPFV.t <- LPFV.t
  out
}
```

These functions can be used as follows:

``` r
doses <- c(0, 0.5, 1, 2, 4, 8) # doses
alRatio <- c(2, 1, 1, 1, 2, 2) # allocation ratio for doses
alpha <- 0.025 # for MCTtest
times <- c(0, 2, 4, 8, 12) # observation times (weeks)
bslMean <- 0 # mean at baseline
errSD <- 0.56 # SD for outcome on absolute scale
RecrDist <- "Quad" # recruitment distribution

rho <- 0.9 # correlation across time for outcome on absolute scale
LPFV.t <- 100 # time of last patient first visit
maxEff <- 0.12 # effect assumed for the highest dose at the last time point
delta_assumed <- calcMean(doses, times, maxEff)[, length(times)]
N <- 236
ia_timing <- 0.5

models <- Mods(
  emax = c(0.5, 1, 2, 4), sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4, 3)),
  quadratic = -0.1, doses = doses
)

oneSim <- sim_one_trial(
  times = times, bslMean = bslMean,
  errSD = errSD, rho = rho, maxEff = maxEff, RecrDist = RecrDist,
  LPFV.t = LPFV.t, models = models, ia_timing = ia_timing, N = N,
  doses = doses, alRatio = alRatio, delta_assumed = delta_assumed,
  alpha = alpha
)

scenRes <- run_scen(
  n_sim = 10, rho = 0.9, maxEff = 0.12,
  LPFV.t = 50, ia_timing = seq(0.3, 0.8, 0.05), delta_assumed = delta_assumed
)
```

## References

Bornkamp, B., Zhou, J., Xi, D., and Cao, W. (2024), “[Futility analyses
for the MCP-mod methodology based on longitudinal
models](https://arxiv.org/abs/2406.19965).”
