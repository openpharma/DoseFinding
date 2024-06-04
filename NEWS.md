# DoseFinding 1.1-1
* Big thanks to Marius Thomas for adding the bMCTtest function,
implementing a generalized version of the Bayesian MCP-Mod
methodology from Fleischer et al (2022) (https://doi.org/10.1002/pst.2193)
* Thanks to Sebastian Bossert for feedback on bMCTtest
* Function critVal is now exported

# DoseFinding 1.0-5
* Fixed bug in non-exported function powMCTBinCount, for situation
when user defined contrast matrix is handed over.
* Added function DesignMCPModApp which starts the externally hosted
R-Shiny app DesignMCPMod to perform power and sample size
calculations (main authors of the app are Sophie Sun and Danyi
Xiong).

# DoseFinding 1.0-4
* Added non-exported function powMCTBinCount, for power calculation
for binary and count data.

# DoseFinding 1.0-3
* Stop and throw error if calculated df=0 in powMCT (passed to
mvtnorm where df=0 implies use of a normal distribution)
* Added function plotMods to plot candidate models using ggplot2
and function plotContr to plot optimal contrasts using ggplot2
(thanks to Sophie Sun for testing and feedback)
* Added documentation for quadratic model (definition of delta)
* Fixed local options in guesst for logistic and sigEmax models

# DoseFinding 1.0-2
* Define USE_FC_LEN_T and add length of character arguments in
Fortran code called from C, to reflect recent changes in gfortran.
* Fix incorrect error message in fitMod (in case placAdj = TRUE and
data are handed over in a data frame via data argument)

# DoseFinding 1.0-1
* Big thanks to Ludger Sandig, who was instrumental in adding
vignettes for practical MCP-Mod implementation guidance;
introducing tests based on testthat and further bug fixes.
* Thanks to Dong Xi, Hilke Kracker for review of earlier versions of
the draft vignettes
* Thanks to Julia Duda for her helpful comments on the package


# DoseFinding 0.9-17
* Added citation to DESCRIPTION file
* Removed alpha argument for pValues function (not used)
* Propagate error messages from mvtnorm in pValues function
(e.g. cov-matrix not psd), (thx to Daisy Bai)
* Make direction attribute in Mods object unique (thx to Yuhan Li)

# DoseFinding 0.9-16
* Fixed minor bug in print.summary.planMod

# DoseFinding 0.9-15
* Mods Added parameter names for all models in the output list
(thanks to Dong Xi for catching this)

DoseFinding 0.9-14
* planMod.Rd Documentation slightly extended.
* qmvtDF moved back to qmvt function from mvtnorm, as problems
in mvtnorm are fixed.

DoseFinding 0.9-13
* projPrBnds now also covers the case when parameter was exactly
on the bound
* bFitMod doseNam changes
* critVal Added self-written qmvt function qmvtDF	(as mvtnorm::qmvt
got instable on Windows 32bit from release 1.0-3), hopefully
superfluous once mvtnorm fixes this.

DoseFinding 0.9-12
* glycobrom dataset: Included column for number of observations
per treatment.
* calcCrit now takes into account "nold" in determining whether
enough design points were specified to be able to calculate
the design criteria.
* bFitMod documentation for plot.bFitMod and predict.bFitMod
methods added. coef.bFitMod method added. Thanks to Lieven Nils
Kennes for pointing towards the issue.

DoseFinding 0.9-11
* Mods Introduce fullMod argument to allow specification
of full model parameters (again).
* calcTDgrad now calculates the analytical gradient for TD
optimal designs for the beta model. The previous numerical
gradient could get unstable for particular parameter values.
Thanks to Tobias Mielke for the calculations!
* planMod.Rd, powMCT.Rd More description on what "sigma" is
* optDesign, optContr Catch Mods objects with multiple direction
attributes properly in these functions.	

# DoseFinding 0.9-10
* plot.MCPMod In case of no significant model, do not plot anything.
* optContr Bugfix in function constOptC, previous algorithm selected
in some situation an incorrect active set (and hence a suboptimal
solution), the current implementation uses quadratic programming
(hence the new suggested package quadprog).

# DoseFinding 0.9-9
* bFitMod.Bayes Stop if starting values lie outside of bounds
specified for the prior distribution
* predict.bFitMod Remove incorrect "if" statement (use
"effect-curve" not "EffectCurve")
* fitModels.bndnls now uses narrowed bounds for 1d models again
(as in 0.9-5 and earlier), thanks to Tobias Mielke for reporting
the three points above.
* optContr now allows for constrained contrasts, i.e. where the
contrast coefficients in placebo and active treatment groups are
required to have different signs.

# DoseFinding 0.9-8
* MCPMod Major changes needed (also in fitMod and MCTtest) to allow
for dose/response names different from "dose", "resp" when a
data-frame is specified (the problem existed as MCTtest, fitMod 
were called from inside MCPMod).
* bFitMod.Bayes Ensure that the starting values for the parameters
are within the bounds specified by the prior (if no starting
values are specified). Thanks to Tobias Mielke for reporting this.
* bFitMod.bootstrap Remove bug for model = "linear" and placAdj = TRUE.
Thanks to Tobias Mielke for reporting this.

# DoseFinding 0.9-7
* fitMod ensure that the data set returned with DRMod objects
is in the original order (not sorted by dose). Also ensure the right
S matrix is used for fitting for type = "general" and unsorted
dose, resp.
* MCTtest fixed problems for type = "general" and unsorted
dose, resp.
* glycobrom Added glycobrom data set 
* planMod Added planning functions for non-linear modelling
* Coded calculations of compositions to be able to remove dependency
on the partitions package
* man files: added reference to paper on generalized MCPMod
* plot.DRMod Minor changes to ensure raw means are always inside
the plotting region (for plotData = "means")

# DoseFinding 0.9-6
* optDesign Re-named "fmodels" argument to "models".
* optDesign for solnp if lowbnd and uppbnd are specified now use a
feasible starting value (otherwise solnp might get into problems).
* plot.DRMod, plot.MCPMod now use lattice graphics
* powMCT removed bug in case of placAdj = TRUE (thanks to Tobias
Mielke for reporting this)
* ess.mcmc minor change to avoid occasional NA's
* Mods removed class c("Mods", "standMod"), now there is only a
class "Mods", this changes the API of MCTtest, optContr and
MCPMod function (direction argument no longer needed, as this info is now 
contained in the "Mods" object).
* neurodeg added the simulated longitudinal dose-finding data set neurodeg
* targN catch incorrect matrix dimension, when in case of only
one alternative model
* fitModel.bndnls old version used narrowed bnds for 1-dim model, when a
starting value was supplied manually (instead of calculated via
optGrid); fixed.
* MCTtest re-name of p-value column to "adj-p".

# DoseFinding 0.9-5
* targN, powN added function targN to evaluate a target function
(e.g. power) for different sample sizes (similar to the old
powerMM function). powN is a convenience function for 
multiple contrast tests using the power.
* sampSizeMCT added convenience function for sample size calculation
for multiple contrast tests using the power.
* optContr Re-named "weights" argument to "w"

# DoseFinding 0.9-4
* TD, ED Fixed bug for model = linInt and placAdj = TRUE
* powMCT Fixed bug for nr(altModels)=1 in case placAdj=TRUE
* Mods Add requirement that placebo dose needs to be included
* print.bFitMod Do not show n.eff for bootstrap samples
* ess.mcmc return NA, if there is just one unique value in chain
* fitMod, MCTtest catch situations, where type = "normal" and 
placAdj	= TRUE
* bFitMod fixed bug for column names of linear model in case of
placAdj = TRUE
* MCPMod: Fixed sign error in model selection, when critV was specified

# DoseFinding 0.9-3
* fitMod Improvements of efficiency (removed calls to do.call 
in optLoc)
* MCPMod passes direction argument now also to TD
* optDesign solnp is now the default optimizer
* calcCrit default for arg designCrit in calcCrit changed (to harmonize
calcCrit and optDesign)
* bFitMod use fitMod.raw in bFitMod.bootstrap (for efficiency)
* critVal Remove contMat argument (was unused)
* powMCT Allow power calculation for placebo adjusted data

# DoseFinding 0.9-1
* Complete re-structuring and tidying up of the package. 
Main ideas: (i) smoother integration of g-functions (ii) focus on
core functionality (iii) more general code/easier extensibility.
* New features: Bayesian dose-response fitting, nicer plots,
optimal designs for non-normal responses, ...
* Special Thanks to Tobias Mielke for testing of the package and 
numerous bug reports.
* Previous versions of the source are available under 
http://cran.r-project.org/package=DoseFinding under 
"Old sources", a Windows binary of the last version 
before the changes is available under http://goo.gl/p1UZ7.

# DoseFinding 0.6-3
* Added PACKAGE = "DoseFinding" to ".C" calls

# DoseFinding 0.6-2
* calcOptDesign partial rewrite of optDes.c and optDesign.R
to fix segfault bug.

# DoseFinding 0.6-1
* vcov.gDRMod is now functional, predict.gDRMod now allows
calculation of confidence intervals
* gFitDRModel minor changes in underlying optimizer
* explicitly export the gradients of the model functions now

# DoseFinding 0.5-7
* gFitDRModel now always returns an estimate (either the
best value from nlminb or from the grid search if nlminb fails)
* gMCPtest: use sigma = corMat instead of corr = corMat in
p/qmvnorm calls (mvtnorm complained in 1-dimensional case)
* gFitDRModel: Introduced default for bnds argument.
* plot.MCPMod: Plot clinRel in the right place, when direction is
equal to "decreasing" (thanks to Jan Rekowski)
* planMM, critVal: When vCov is specified now right correlation 
matrix is calculated
* calcOptDesign: Additional argument (standDopt) to allow for optional
standardization (division of the log determinant by the number 
of parameters) of the D-optimal design criterion. 

# DoseFinding 0.5-6
* getGrid corrected bug for Ngrd > 75025
* calcOptDesign: For method = "exact" and n2 > 0 the function
did not return the optimal incremental design but the
overall optimal design

# DoseFinding 0.5-5
* gFitDRModel can now fit dose-response models without
intercept
* gMCPtest minor changes to allow for user defined contrast
matrix

# DoseFinding 0.5-4
* MCPtest now uses correct degrees of freedom if addCovars != ~1
* Feedback from Andreas Krause led to a number 
smaller changes in the package (e.g., plot.(g)DRMod or 
fitDRModel). Thanks Andreas!
* Print lattice plots explicitly to increase compability
with Sweave.

# DoseFinding 0.5-3
* Ensure in rndDesign that N is recognized as
an integer by using N <- round(N), to avoid floating point
problems.
* Remove naming bug in gFitDRModel (drFit instead of drEst)

# DoseFinding 0.5-2
* Corrected bug in b-vector for sigEmax model (calcBvec,
only affected MED-type optimal designs)
* Included INDEX file to order the overview help-page better
* predict.DRMod now stops when type = "fullModel" and the
argument newdata does not contain values for all variables
specified in addCovars (thanks to Mouna). 

# DoseFinding 0.5-1
* Restructured calcOptDesign function to allow for user
defined criteria functions.
* The MCPMod object now always contains a estDose and fm
entry (which is NA in case of non-significance or non-convergence)
* Added generalized fitting code, variances and covariances of 
estimates are not available at the moment.
* Added vCov argument to planMM, sampSize, powerMM (so it is possible
to take into account covariances when calculating optimal contrasts)
* Changed order in trellis plots in plotModels (order as specified in
models list instead of alphanumerical order)
* Restructured and summarized help pages 
* Removed dependency on numDeriv package (only suggested now),
this is only needed for calculating optimal designs involving
the beta model.

# DoseFinding 0.4-3
* Minor change in Makevars file (so that DoseFinding works on 
Solaris).

# DoseFinding 0.4-2
* calcBayesEst, getUpdDesign: Minor changes to make functions
more suited for general purpose use.

# DoseFinding 0.4-1
* Introduced new functions calcBayesEst and getUpdDesign,
both were used for simulation purposes in the paper
Bornkamp et al. (2011) "Response Adaptive Dose-Finding
under Model Uncertainty" (to appear in Annals of Applied
Statistics).

# DoseFinding 0.3-1
* calcOptDesign now has an additional optimizer "exact". This
methods calculates all possible designs for a given sample size
and then selects the best.
* Changed order in MakeVars as requested
* calcCrit now checks whether there are not less design points
than parameters.
* Code now checks for positive sigma in powCalc and powerMM

# DoseFinding 0.2-3
* MED function now checks for clinRel > 0 (thanks to Georgina).
* Changed minor bug in output from print.MCPtest (print one-sided
just once)
* Code now outputs a warning, when 'models' argument is missing
(in MCPMod and fullMod function); in fitDRModel it outputs 
a warning if 'model' is missing
* Introduced a default base = 0 and maxEff = 1 for the plotModels
function.
* Added a summary method for DRMod objects.
* Removed superfluous addCovarVals argument from predict.DRMod
* Removed option method = "mult" in calcOptDesign

# DoseFinding 0.2-2
* calcCrit and calcOptDesign now check for NA, NaN or +-Inf
values in the gradient and bvector (and stop execution when
these values occur) before passing these values to the C
code.
* Introduced a logLik method for DRMod objects
* Changed mvtnorm.control default argument for "interval"
to reflect recent changes in mvtnorm package.

# DoseFinding 0.2-1
* Made the getGrad function (gradient for dose-response model),
including documentation available for end-user (was previously
hidden in NAMESPACE)
* Changes in the plot.MCPMod function (col argument for 
panel.superpose was read in different order depending on
lattice options, now there is a manual workaround with
panel.xyplot calls for each group)

# DoseFinding 0.1-3
* Smaller changes in calcCrit functions (the parameter p is
now calculated by the nPars function as in getOptDesign)
* Add further options to powerScenario function (now possible
to for user-specified row and column names for output matrix)

# DoseFinding 0.1-2
* Removed one example from sampSize to reduce check time.
* modelSelect: Use first model when two models have exactly
the same AIC or BIC value.
* predict.DRMod: Return NA as standard deviation if code
cannot calculate Cholesky transformation of covariance
matrix (thanks to Setia Pramana for the hint).
* calcCrit: Code now allows for specifying multiple designs
in a matrix.

# DoseFinding 0.1-1
* fitModel.nls now checks whether nls with plinear option made
a positive number of iterations (as additional convergence check).
In some cases (eg. when number of parameters = number of 
doses) plinear does not do any iteration and does *not* 
put out a warning message that the algorithm failed.
* The calcOptDesign function now allows for upper and lower bounds
on the allocation weights.
* There is no longer the need to specify clinRel, when one wants
to calculate a D-optimal design.
* Output of bootMCPMod function in case of model averaging now
also includes dose estimates under each model & corrected
bug in print.bootMCPMod function (thanks to Setia Pramana)

# DoseFinding 0.1
* 1st Release as version 0.1. Improvements over MCPMod package:
- Extended and improved version of MCPMod (allowing for covariates
and robustified self-developed optimizer)
- Functions for MCP (MCPtest) and Modelling (fitDRModel) part
now available to the user
- New functions (eg. bootMCPMod, powerScenario)
- Functions for calculating optimal designs

