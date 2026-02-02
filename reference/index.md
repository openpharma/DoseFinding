# Package index

## Package

- [`DoseFinding`](https://openpharma.github.io/DoseFinding/reference/DoseFinding-package.md)
  [`DoseFinding-package`](https://openpharma.github.io/DoseFinding/reference/DoseFinding-package.md)
  : DoseFinding: Planning and Analyzing Dose Finding Experiments

## Design dose-reponse trials

- [`Mods()`](https://openpharma.github.io/DoseFinding/reference/Mods.md)
  [`getResp()`](https://openpharma.github.io/DoseFinding/reference/Mods.md)
  [`plotMods()`](https://openpharma.github.io/DoseFinding/reference/Mods.md)
  [`plot(`*`<Mods>`*`)`](https://openpharma.github.io/DoseFinding/reference/Mods.md)
  : Define dose-response models
- [`guesst()`](https://openpharma.github.io/DoseFinding/reference/guesst.md)
  : Calculate guesstimates based on prior knowledge
- [`optContr()`](https://openpharma.github.io/DoseFinding/reference/optContr.md)
  [`plot(`*`<optContr>`*`)`](https://openpharma.github.io/DoseFinding/reference/optContr.md)
  [`plotContr()`](https://openpharma.github.io/DoseFinding/reference/optContr.md)
  : Calculate optimal contrasts
- [`optDesign()`](https://openpharma.github.io/DoseFinding/reference/optDesign.md)
  [`calcCrit()`](https://openpharma.github.io/DoseFinding/reference/optDesign.md)
  [`rndDesign()`](https://openpharma.github.io/DoseFinding/reference/optDesign.md)
  [`plot(`*`<DRdesign>`*`)`](https://openpharma.github.io/DoseFinding/reference/optDesign.md)
  : Function to calculate optimal designs
- [`planMod()`](https://openpharma.github.io/DoseFinding/reference/planMod.md)
  [`summary(`*`<planMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/planMod.md)
  [`plot(`*`<planMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/planMod.md)
  : Evaluate performance metrics for fitting dose-response models
- [`powMCT()`](https://openpharma.github.io/DoseFinding/reference/powMCT.md)
  : Calculate power for multiple contrast test
- [`sampSize()`](https://openpharma.github.io/DoseFinding/reference/sampSize.md)
  [`sampSizeMCT()`](https://openpharma.github.io/DoseFinding/reference/sampSize.md)
  [`targN()`](https://openpharma.github.io/DoseFinding/reference/sampSize.md)
  [`powN()`](https://openpharma.github.io/DoseFinding/reference/sampSize.md)
  [`plot(`*`<targN>`*`)`](https://openpharma.github.io/DoseFinding/reference/sampSize.md)
  : Sample size calculations
- [`DesignMCPModApp()`](https://openpharma.github.io/DoseFinding/reference/DesignMCPModApp.md)
  : Start externally hosted DesignMCPMod Shiny App

## Analyze dose-response trials

- [`MCPMod()`](https://openpharma.github.io/DoseFinding/reference/MCPMod.md)
  [`predict(`*`<MCPMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/MCPMod.md)
  [`plot(`*`<MCPMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/MCPMod.md)
  : MCPMod - Multiple Comparisons and Modeling
- [`MCTpval()`](https://openpharma.github.io/DoseFinding/reference/MCTpval.md)
  : Calculate multiplicity adjusted p-values for multiple contrast test
- [`MCTtest()`](https://openpharma.github.io/DoseFinding/reference/MCTtest.md)
  : Performs multiple contrast test
- [`bMCTtest()`](https://openpharma.github.io/DoseFinding/reference/bMCTtest.md)
  : Performs Bayesian multiple contrast test
- [`powMCTInterim()`](https://openpharma.github.io/DoseFinding/reference/powMCTInterim.md)
  : Calculate Conditional or Predictive Power for Multiple Contrast Test
- [`critVal()`](https://openpharma.github.io/DoseFinding/reference/critVal.md)
  : Calculate critical value for multiple contrast test
- [`mvpostmix()`](https://openpharma.github.io/DoseFinding/reference/mvpostmix.md)
  : Prior to posterior updating for a multivariate normal mixture
- [`mvtnorm.control()`](https://openpharma.github.io/DoseFinding/reference/mvtnorm-control.md)
  : Control options for pmvt and qmvt functions

## Fit dose-response model

- [`bFitMod()`](https://openpharma.github.io/DoseFinding/reference/bFitMod.md)
  [`predict(`*`<bFitMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/bFitMod.md)
  [`plot(`*`<bFitMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/bFitMod.md)
  [`coef(`*`<bFitMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/bFitMod.md)
  : Fit a dose-response model using Bayesian or bootstrap methods.

- [`defBnds()`](https://openpharma.github.io/DoseFinding/reference/defBnds.md)
  : Calculates default bounds for non-linear parameters in dose-response
  models

- [`emax()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`emaxGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`sigEmax()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`sigEmaxGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`exponential()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`exponentialGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`quadratic()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`quadraticGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`betaMod()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`betaModGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`linear()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`linearGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`linlog()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`linlogGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`logistic()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`logisticGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`linInt()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  [`linIntGrad()`](https://openpharma.github.io/DoseFinding/reference/drmodels.md)
  : Built-in dose-response models in DoseFinding

- [`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  [`coef(`*`<DRMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  [`vcov(`*`<DRMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  [`predict(`*`<DRMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  [`plot(`*`<DRMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  [`logLik(`*`<DRMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  [`AIC(`*`<DRMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  [`gAIC(`*`<DRMod>`*`)`](https://openpharma.github.io/DoseFinding/reference/fitMod.md)
  : Fit non-linear dose-response model

- [`maFitMod()`](https://openpharma.github.io/DoseFinding/reference/maFitMod.md)
  [`predict(`*`<maFit>`*`)`](https://openpharma.github.io/DoseFinding/reference/maFitMod.md)
  [`plot(`*`<maFit>`*`)`](https://openpharma.github.io/DoseFinding/reference/maFitMod.md)
  : Fit dose-response models via bootstrap model averaging (bagging)

- [`TD()`](https://openpharma.github.io/DoseFinding/reference/targdose.md)
  [`ED()`](https://openpharma.github.io/DoseFinding/reference/targdose.md)
  :

  Calculate dose estimates for a fitted dose-response model (via
  [`fitMod()`](https://openpharma.github.io/DoseFinding/reference/fitMod.md),
  [`bFitMod()`](https://openpharma.github.io/DoseFinding/reference/bFitMod.md))
  or
  [`maFitMod()`](https://openpharma.github.io/DoseFinding/reference/maFitMod.md))
  or a
  [`Mods()`](https://openpharma.github.io/DoseFinding/reference/Mods.md)
  object

## Datasets

- [`biom`](https://openpharma.github.io/DoseFinding/reference/biom.md) :
  Biometrics Dose Response data
- [`glycobrom`](https://openpharma.github.io/DoseFinding/reference/glycobrom.md)
  : Glycopyrronium Bromide dose-response data
- [`IBScovars`](https://openpharma.github.io/DoseFinding/reference/IBScovars.md)
  : Irritable Bowel Syndrome Dose Response data with covariates
- [`migraine`](https://openpharma.github.io/DoseFinding/reference/migraine.md)
  : Migraine Dose Response data
- [`neurodeg`](https://openpharma.github.io/DoseFinding/reference/neurodeg.md)
  : Neurodegenerative disease simulated longitudinal dose-finding data
  set
