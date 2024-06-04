
#' @description
#' The DoseFinding package provides functions for the design and analysis of dose-finding experiments (for example pharmaceutical Phase II clinical trials). It provides functions for: multiple
#' contrast tests (MCTtest), fitting non-linear dose-response models (fitMod), a combination of testing and dose-response modelling (MCPMod), and calculating optimal designs (optDesign), both for
#' normal and general response variable.
#' 
#' @details
#' The main functions are:\cr
#' \bold{MCTtest}: Implements a multiple contrast tests\cr
#' \bold{powMCT}: Power calculations for multiple contrast tests\cr
#' \bold{fitMod}: Fits non-linear dose-response models\cr
#' \bold{optDesign}: Calculates optimal designs for dose-response models\cr
#' \bold{MCPMod}: Performs MCPMod methodology\cr
#' \bold{sampSize}: General function for sample size calculation\cr
#' 
#' @references Bornkamp, B., Bretz, F., Dette, H. and Pinheiro, J. C. (2011).
#' Response-Adaptive Dose-Finding under model uncertainty, \emph{Annals of
#' Applied Statistics}, \bold{5}, 1611--1631
#' 
#' Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An R Package for
#' the Design and Analysis of Dose-Finding Studies, \emph{Journal of
#' Statistical Software}, \bold{29}(7), 1--23
#' 
#' Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining multiple
#' comparisons and modeling techniques in dose-response studies,
#' \emph{Biometrics}, \bold{61}, 738--748
#' 
#' Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro, J. C. (2008). Optimal
#' Designs for Dose Finding Studies, \emph{Journal of the American Statisical
#' Association}, \bold{103}, 1225--1237
#' 
#' O'Quigley, J., Iasonos, A. and Bornkamp, B. (2017) Handbook of methods for
#' designing, monitoring, and analyzing dose-finding trials, CRC press, Part 3:
#' Dose-Finding Studies in Phase II
#' 
#' Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of
#' dose finding studies combining multiple comparisons and modeling procedures,
#' \emph{Journal of Biopharmaceutical Statistics}, \bold{16}, 639--656
#' 
#' Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014) Model-based
#' dose finding under model uncertainty using general parametric models,
#' \emph{Statistics in Medicine}, \bold{33}, 1646--1661
#' 
#' Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley
#' @keywords internal
#' @examples
#' 
#' data(IBScovars)
#' 
#' ## perform (model based) multiple contrast test
#' ## define candidate dose-response shapes
#' models <- Mods(linear = NULL, emax = 0.2, quadratic = -0.17,
#'                doses = c(0, 1, 2, 3, 4))
#' ## plot models
#' plot(models)
#' ## perform multiple contrast test
#' test <- MCTtest(dose, resp, IBScovars, models=models,
#'                 addCovars = ~ gender)
#' 
#' ## fit non-linear emax dose-response model
#' fitemax <- fitMod(dose, resp, data=IBScovars, model="emax",
#'                   bnds = c(0.01,5))
#' ## display fitted dose-effect curve
#' plot(fitemax, CI=TRUE, plotData="meansCI")
#' 
#' ## Calculate optimal designs for target dose (TD) estimation
#' doses <- c(0, 10, 25, 50, 100, 150)
#' fmodels <- Mods(linear = NULL, emax = 25, exponential = 85,
#'                 logistic = c(50, 10.8811),
#'                 doses = doses, placEff=0, maxEff=0.4)
#' plot(fmodels, plotTD = TRUE, Delta = 0.2)
#' weights <- rep(1/4, 4)
#' desTD <- optDesign(fmodels, weights, Delta=0.2, designCrit="TD")
#' 
"_PACKAGE"


#' Built-in dose-response models in DoseFinding
#' 
#' @description
#' Dose-response model functions and gradients.
#' 
#' Below are the definitions of the model functions:
#' 
#' \bold{Emax model} \deqn{ }{f(d,theta)=E0+Emax d/(ED50 + d).}\deqn{
#' f(d,\theta)=E_0+E_{max}\frac{d}{ED_{50}+d}}{f(d,theta)=E0+Emax d/(ED50 +
#' d).}
#' 
#' \bold{Sigmoid Emax Model} \deqn{ }{f(d,theta)=E0+Emax d^h/(ED50^h +
#' d^h).}\deqn{
#' f(d,\theta)=E_0+E_{max}\frac{d^h}{ED^h_{50}+d^h}}{f(d,theta)=E0+Emax
#' d^h/(ED50^h + d^h).}
#' 
#' \bold{Exponential Model} \deqn{ }{f(d,theta)=E0+E1 (exp(d/delta)-1).}\deqn{
#' f(d,\theta)=E_0+E_1(\exp(d/\delta)-1)}{f(d,theta)=E0+E1 (exp(d/delta)-1).}
#' 
#' \bold{Beta model} \deqn{ }{f(d,theta)=E0+Emax
#' B(delta1,delta2)(d/scal)^delta1(1-d/scal)^delta2}\deqn{
#' f(d,\theta)=E_0+E_{max}B(\delta_1,\delta_2)(d/scal)^{\delta_1}(1-d/scal)^{\delta_2}
#' }{f(d,theta)=E0+Emax B(delta1,delta2)(d/scal)^delta1(1-d/scal)^delta2}\deqn{
#' }{f(d,theta)=E0+Emax B(delta1,delta2)(d/scal)^delta1(1-d/scal)^delta2} here
#' \deqn{B(\delta_1,\delta_2)=(\delta_1+\delta_2)^{\delta_1+\delta_2}/(\delta_1^{\delta_1}
#' }{B(delta1,delta2)=(delta1+delta2)^(delta1+delta2)/(delta1^delta1
#' delta2^delta2).}\deqn{
#' \delta_2^{\delta_2})}{B(delta1,delta2)=(delta1+delta2)^(delta1+delta2)/(delta1^delta1
#' delta2^delta2).} and \eqn{scal}{scal} is a fixed dose scaling parameter.
#' 
#' \bold{Linear Model} \deqn{ }{f(d,theta)=E0+delta d.}\deqn{
#' f(d,\theta)=E_0+\delta d}{f(d,theta)=E0+delta d.}
#' 
#' \bold{Linear in log Model} \deqn{ }{f(d,theta)=E0+delta log(d + off),}\deqn{
#' f(d,\theta)=E_0+\delta \log(d + off)}{f(d,theta)=E0+delta log(d + off),}
#' here \eqn{off}{off} is a fixed offset parameter.
#' 
#' \bold{Logistic Model} \deqn{
#' f(d, \theta) = E_0 + E_{\max}/\left\{1 + \exp\left[ \left(ED_{50} - d
#' \right)/\delta \right] \right\}}{f(d,theta)=E0+Emax/(1 + exp((ED50-d)/delta)).}
#' 
#' \bold{Quadratic Model} \deqn{ }{f(d,theta)=E0+beta1 d+beta2 d^2.}\deqn{
#' f(d,\theta)=E_0+\beta_1d+\beta_2d^2}{f(d,theta)=E0+beta1 d+beta2 d^2.} The
#' standardized model equation for the quadratic model is \eqn{d+\delta
#' d^2}{d+delta d^2}, with \eqn{\delta=\beta_2/\beta_1}{delta=beta2/beta1}.
#' 
#' \bold{Linear Interpolation model}\cr The linInt model provides linear
#' interpolation at the values defined by the nodes vector. In virtually all
#' situations the nodes vector is equal to the doses used in the analysis. For
#' example the \code{\link{Mods}} and the \code{\link{fitMod}} function
#' automatically use the doses that are used in the context of the function
#' call as nodes. The guesstimates specified in the \code{\link{Mods}} function
#' need to be the treatment effects at the active doses standardized to the
#' interval [0,1] (see the examples in the \code{\link{Mods}} function).
#' 
#' @details
#' The \bold{Emax model} is used to represent monotone, concave dose-response
#' shapes.  To distinguish it from the more general sigmoid emax model it is
#' sometimes also called hyperbolic emax model.
#' 
#' The \bold{sigmoid Emax} model is an extension of the (hyperbolic) Emax model
#' by introducing an additional parameter h, that determines the steepness of
#' the curve at the ed50 value. The sigmoid Emax model describes monotonic,
#' sigmoid dose-response relationships. In the toxicology literature this model
#' is also called four-parameter log-logistic (4pLL) model.
#' 
#' The \bold{quadratic} model is intended to capture a possible non-monotonic
#' dose-response relationship.
#' 
#' The \bold{exponential model} is intended to capture a possible sub-linear or
#' a convex dose-response relationship.
#' 
#' The \bold{beta model} is intended to capture non-monotone dose-response
#' relationships and is more flexible than the quadratic model.  The kernel of
#' the beta model function consists of the kernel of the density function of a
#' beta distribution on the interval [0,scal]. The parameter scal is not
#' estimated but needs to be set to a value larger than the maximum dose. It
#' can be set in most functions (\samp{fitMod}, \samp{Mods}) via the
#' \samp{addArgs} argument, when omitted a value of \samp{1.2*(maximum dose)}
#' is used as default, where the maximum dose is inferred from other input to
#' the respective function.
#' 
#' The \bold{linear in log-dose} model is intended to capture concave shapes.
#' The parameter \code{off} is not estimated in the code but set to a
#' pre-specified value. It can be set in most functions (\samp{fitMod},
#' \samp{Mods}) via the \samp{addArgs} argument, when omitted a value of
#' \samp{0.01*(maximum dose)} is used as default, where the maximum dose is
#' inferred from other input to the respective function.
#' 
#' The \bold{logistic model} is intended to capture general monotone, sigmoid
#' dose-response relationships. The logistic model and the sigmoid Emax model
#' are closely related: The sigmoid Emax model is a logistic model in
#' log(dose).
#' 
#' The \bold{linInt model} provids linear interpolation of the means at the
#' doses. This can be used as a "nonparametric" estimate of the dose-response
#' curve, but is probably most interesting for specifying a "nonparametric"
#' truth during planning and assess how well parametric models work under a
#' nonparametric truth. For the function \samp{Mods} and \samp{fitMod} the
#' interpolation \samp{nodes} are selected equal to the dose-levels specified.
#' 
#' @name drmodels
#' @rdname drmodels
#' @aliases drmodels betaMod emax sigEmax exponential logistic linear linlog
#' quadratic linInt betaModGrad emaxGrad sigEmaxGrad exponentialGrad
#' logisticGrad linearGrad linlogGrad quadraticGrad linIntGrad
#' @param dose Dose variable
#' @param e0 For most models placebo effect. For logistic model left-asymptote
#' parameter, corresponding to a basal effect level (not the placebo effect)
#' @param eMax Beta Model: Maximum effect within dose-range\cr Emax, sigmoid
#' Emax, logistic Model: Asymptotic maximum effect
#' @param ed50 Dose giving half of the asymptotic maximum effect
#' @param h Hill parameter, determining the steepness of the model at the ED50
#' @param e1 Slope parameter for exponential model
#' @param delta Exponential model: Parameter, controlling the convexity of the
#' model.\cr Linear and linlog model: Slope parameter\cr Logistic model:
#' Parameter controlling determining the steepness of the curve
#' @param delta1 delta1 parameter for beta model
#' @param delta2 delta2 parameter for beta model
#' @param b1 first parameter of quadratic model
#' @param b2 second parameter of quadratic model (controls, whether model is
#' convex or concave)
#' @param resp Response values at the nodes for the linInt model
#' @param off Offset value to avoid problems with dose=0 (treated as a fixed
#' value, not estimated)
#' @param scal Scale parameter (treated as a fixed value, not estimated)
#' @param nodes Interpolation nodes for the linear interpolation for the linInt
#' model (treated as a fixed value, not estimated)
#' @param ...  Just included for convenience in the gradient functions, so that
#' for example \code{quadratic(dose, e0=0, b1=1, b2=3)} will not throw an error
#' (although the gradient of the quadratic model is independent of e0, b1 and
#' b2).
#' @return Response value for model functions or matrix containing the gradient
#' evaluations.
#' @seealso \code{\link{fitMod}}
#' @references MacDougall, J. (2006). Analysis of dose-response studies - Emax
#' model,\emph{in} N. Ting (ed.), \emph{Dose Finding in Drug Development},
#' Springer, New York, pp. 127--145
#' 
#' Pinheiro, J. C., Bretz, F. and Branson, M. (2006). Analysis of dose-response
#' studies - modeling approaches, \emph{in} N. Ting (ed.). \emph{Dose Finding
#' in Drug Development}, Springer, New York, pp. 146--171
#' @examples
#' 
#' ## some quadratic example shapes
#' quadModList <- Mods(quadratic = c(-0.5, -0.75, -0.85, -1), doses = c(0,1))
#' plotMods(quadModList)
#' 
#' ## some emax example shapes
#' emaxModList <- Mods(emax = c(0.02,0.1,0.5,1), doses = c(0,1))
#' plotMods(emaxModList)
#' ## example for gradient
#' emaxGrad(dose = (0:4)/4, eMax = 1, ed50 = 0.5)
#' 
#' ## some sigmoid emax example shapes
#' sigEmaxModList <- Mods(sigEmax = rbind(c(0.05,1), c(0.15,3), c(0.4,8),
#'                        c(0.7,8)), doses = c(0,1))
#' plotMods(sigEmaxModList)
#' sigEmaxGrad(dose = (0:4)/4, eMax = 1, ed50 = 0.5, h = 8)
#' 
#' ## some exponential example shapes
#' expoModList <- Mods(exponential = c(0.1,0.25,0.5,2), doses=c(0,1))
#' plotMods(expoModList)
#' exponentialGrad(dose = (0:4)/4, e1 = 1, delta = 2)
#' 
#' ## some beta model example shapes
#' betaModList <- Mods(betaMod = rbind(c(1,1), c(1.5,0.75), c(0.8,2.5),
#'                     c(0.4,0.9)), doses=c(0,1), addArgs=list(scal = 1.2))
#' plotMods(betaModList)
#' betaModGrad(dose = (0:4)/4, eMax = 1, delta1 = 1, delta2 = 1, scal = 5)
#' 
#' ## some logistic model example shapes
#' logistModList <- Mods(logistic = rbind(c(0.5,0.05), c(0.5,0.15),
#'                       c(0.2,0.05), c(0.2,0.15)), doses=c(0,1))
#' plotMods(logistModList)
#' logisticGrad(dose = (0:4)/4, eMax = 1, ed50 = 0.5, delta = 0.05)
#' 
#' ## some linInt shapes
#' genModList <- Mods(linInt = rbind(c(0.5,1,1),
#'                       c(0,1,1), c(0,0,1)), doses=c(0,0.5,1,1.5))
#' plotMods(genModList)
#' linIntGrad(dose = (0:4)/4, resp=c(0,0.5,1,1,1), nodes=(0:4)/4)
#' 
#' 
NULL


## Documentation of datasets

#' Biometrics Dose Response data
#' 
#' An example data set for dose response studies. This data set was used in
#' Bretz et al. (2005) to illustrate the MCPMod methodology.
#' 
#' @name biom
#' @docType data
#' @usage data(biom)
#' @format A data frame with 100 observations on the following 2 variables.
#'   \describe{
#'   \item{\code{resp}}{a numeric vector containing the response values}
#'   \item{\code{dose}}{a numeric vector containing the dose values}
#' }
#' @source Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining
#' multiple comparisons and modeling techniques in dose-response studies,
#' \emph{Biometrics}, \bold{61}, 738--748
#' @keywords datasets
NULL


#' Glycopyrronium Bromide dose-response data
#' 
#' Data from a clinical study evaluating Efficacy and Safety of Four Doses of
#' Glycopyrronium Bromide in Patients With Stable Chronic Obstructive Pulmonary
#' Disease (COPD).  This data set was obtained from clinicaltrials.gov
#' (NCT00501852).  The study design was a 4 period incomplete cross-over
#' design. The primary endpoint is the trough forced expiratory volume in 1
#' second (FEV1) following 7 days of Treatment.
#' 
#' The data given here are summary estimates (least-square means) for each
#' dose.
#' 
#' @name glycobrom
#' @docType data
#' @usage data(glycobrom)
#' @format A data frame with 5 summary estimates (one per dose). Variables:
#'   A data frame with 5 summary estimates (one per dose). Variables:
#'  \describe{
#'    \item{\code{dose}}{a numeric vector containing the dose values}
#'    \item{\code{fev1}}{a numeric vector containing the least square
#'      mean per dose}
#'    \item{\code{sdev}}{a numeric vector containing the standard errors
#'      of the least square means per dose}
#'    \item{\code{n}}{Number of participants analyzed per treatment group}
#'  }
#' @source http://clinicaltrials.gov/ct2/show/results/NCT00501852
#' @keywords datasets
#' @examples
#' 
#'  ## simulate a full data set with given means and sdv (here we ignore
#'   ## the original study was a cross-over design, and simulate a parallel
#'   ## group design)
#'   simData <- function(mn, sd, n, doses, fixed = TRUE){
#'     ## simulate data with means (mns) and standard deviations (sd), for
#'     ## fixed = TRUE, the data set will have observed means and standard
#'     ## deviations as given in mns and sd
#'     resp <- numeric(sum(n))
#'     uppind <- cumsum(n)
#'     lowind <- c(0,uppind)+1
#'     for(i in 1:length(n)){
#'       rv <- rnorm(n[i])
#'       if(fixed)
#'         rv <- scale(rv)
#'       resp[lowind[i]:uppind[i]] <- mn[i] + sd[i]*rv
#'     }
#'     data.frame(doses=rep(doses, n), resp=resp)
#'   }
#'   data(glycobrom)
#'   fullDat <- simData(glycobrom$fev1, glycobrom$sdev, glycobrom$n,
#'                      glycobrom$dose)
#' 
NULL

#' Irritable Bowel Syndrome Dose Response data with covariates
#' 
#' A subset of the data used by (Biesheuvel and Hothorn, 2002).  The data are
#' part of a dose ranging trial on a compound for the treatment of the
#' irritable bowel syndrome with four active treatment arms, corresponding to
#' doses 1,2,3,4 and placebo. Note that the original dose levels have been
#' blinded in this data set for confidentiality. The primary endpoint was a
#' baseline adjusted abdominal pain score with larger values corresponding to a
#' better treatment effect. In total 369 patients completed the study, with
#' nearly balanced allocation across the doses.
#' 
#' 
#' @name IBScovars
#' @docType data
#' @usage data(IBScovars)
#' @format 
#'   A data frame with 369 observations on the following 2 variables.
#'   \describe{
#'     \item{\code{gender}}{a factor specifying the gender}
#'     \item{\code{dose}}{a numeric vector}
#'     \item{\code{resp}}{a numeric vector}
#'   }
#' @source Biesheuvel, E. and Hothorn, L. A. (2002). Many-to-one comparisons in
#' stratified designs, \emph{Biometrical Journal}, \bold{44}, 101--116
#' @keywords datasets
NULL


#' Migraine Dose Response data
#' 
#' Data set obtained from clinicaltrials.gov (NCT00712725).  This was
#' randomized placebo controlled dose-response trial for treatment of acute
#' migraine. The primary endpoint was "pain freedom at 2 hours postdose" (a
#' binary measurement).
#' 
#' 
#' @name migraine
#' @docType data
#' @usage data(migraine)
#' @format 
#'  A data frame with 517 columns corresponding to the patients that
#'  completed the trial
#'  \describe{
#'    \item{\code{dose}}{a numeric vector containing the dose values}
#'    \item{\code{painfree}}{number of treatment responders}
#'    \item{\code{ntrt}}{number of subject per treatment group}
#'  }
#' @source http://clinicaltrials.gov/ct2/show/results/NCT00712725
#' @keywords datasets
NULL



#' Neurodegenerative disease simulated longitudinal dose-finding data set
#' 
#' This simulated data set is motivated by a real Phase 2 clinical study of a
#' new drug for a neurodegenerative disease. The state of the disease is
#' measured through a functional scale, with smaller values corresponding to
#' more severe neurodeterioration. The goal of the drug is to reduce the rate
#' of disease progression, which is measured by the linear slope of the
#' functional scale over time.
#' 
#' The trial design includes placebo and four doses: 1, 3, 10, and 30 mg, with
#' balanced allocation of 50 patients per arm. Patients are followed up for one
#' year, with measurements of the functional scale being taken at baseline and
#' then every three months.
#' 
#' The functional scale response is assumed to be normally distributed and,
#' based on historical data, it is believed that the longitudinal progression
#' of the functional scale over the one year of follow up can be modeled a
#' simple linear trend. See the example below on how to analyse this type of
#' data.
#' 
#' This data set was used in Pinheiro et al. (2014) to illustrate the
#' generalized MCPMod methodology.
#' 
#' 
#' @name neurodeg
#' @docType data
#' @usage data(neurodeg)
#' @format 
#'  A data frame with 100 observations on the following 2 variables.
#'  \describe{
#'    \item{\code{resp}}{a numeric vector containing the response values}
#'    \item{\code{dose}}{a numeric vector containing the dose values}
#'    \item{\code{id}}{Patient ID}    
#'    \item{\code{time}}{time of measurement}
#'  }
#' @source Pinheiro, J. C., Bornkamp, B., Glimm, E. and Bretz, F. (2014)
#' Model-based dose finding under model uncertainty using general parametric
#' models, \emph{Statistics in Medicine}, \bold{33}, 1646--1661
#' @keywords datasets
#' @examples
#' 
#' \dontrun{
#' ## reproduce analysis from Pinheiro et al. (2014)
#' data(neurodeg)
#' ## first fit the linear mixed effect model
#' library(nlme)
#' fm <- lme(resp ~ as.factor(dose):time, neurodeg, ~time|id, method = "ML")
#' muH <- fixef(fm)[-1] # extract estimates
#' covH <- vcov(fm)[-1,-1]
#' 
#' ## derive optimal contrasts for candidate shapes
#' doses <- c(0, 1, 3, 10, 30)
#' mod <- Mods(emax = 1.11, quadratic= -0.022, exponential = 8.867,
#'             linear = NULL, doses = doses) # 
#' contMat <- optContr(mod, S=covH) # calculate optimal contrasts
#' ## multiple contrast test
#' MCTtest(doses, muH, S=covH, type = "general", critV = TRUE,
#'         contMat=contMat)
#' ## fit the emax model
#' fitMod(doses, muH, S=covH, model="emax", type = "general",
#'        bnds=c(0.1, 10))
#' 
#' 
#' ## alternatively one can also fit the model using nlme
#' nlme(resp ~ b0 + (e0 + eM * dose/(ed50 + dose))*time, neurodeg,
#'      fixed = b0 + e0 + eM + ed50 ~ 1, random = b0 + e0 ~ 1 | id,
#'      start = c(200, -4.6, 1.6, 3.2))
#' ## both approaches lead to rather similar results
#' }
#' 
NULL


