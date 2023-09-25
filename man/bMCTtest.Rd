\name{bMCTtest}
\alias{bMCTtest}
\title{
  Performs Bayesian multiple contrast test
}
\description{
  This function performs a Bayesian multiple contrast test using normal mixture
  priors for the response on each dose, as proposed in Fleischer et al. (2022). For a general description of the 
  multiple contrast test see \code{\link{MCTtest}}. 
  
  If \samp{type = "normal"}, an ANCOVA model based on a homoscedastic
  normality assumption is fitted and posteriors for dose-response and contrast vectors 
  are obtained assuming a known variance. 
  
  For \samp{type = "general"} it is assumed multivariate normally
  distributed estimates are specified in \samp{resp} with covariance
  given by \samp{S}, which define the likelihood. 
  Posteriors for dose-response and contrast vectors are then obtained assuming a known covariance matrix S
  
  The multiple contrast test decision is based on the maximum posterior probability of a contrast being 
  greater (or smaller) than zero. Thresholds for the posterior probability can either be supplied or
  will be derived from frequentist critical values. In the latter case the Bayesian test will give approximately
  the same results as the frequentist multiple contrast test if uninformative priors are used.
}
\usage{
bMCTtest(dose, resp, data = NULL, models, S = NULL, type = c("normal", "general"), 
         prior, alpha = 0.025, alternative = c("greater", "less"), na.action = na.fail,
         mvtcontrol = mvtnorm.control(), contMat = NULL,
         critV = NULL)
}
\arguments{
    \item{dose, resp}{
    Either vectors of equal length specifying dose and response values,
    or names of variables in the data frame specified in \samp{data}.
  }
  \item{data}{
    Data frame containing the variables referenced in dose and resp if
    \samp{data} is not specified it is assumed that \samp{dose} and
    \samp{resp} are variables referenced from data (and no vectors)
  }
  \item{models}{
    An object of class \samp{Mods}, see \code{\link{Mods}} for
    details
  }
  \item{S}{
    The covariance matrix of \samp{resp} when \samp{type = "general"},
    see Description.
  }
  \item{type}{
    Determines whether the likelihood is based on an ANCOVA model under a
    homoscedastic normality assumption (when \samp{type = "normal"}), or
    estimates at the doses and their covariance matrix are specified 
    directly in \samp{resp} and \samp{S}. See also \code{\link{MCTtest}} and Pinheiro et al. (2014).
  }
  \item{prior}{
    List of length equal to the number of doses with the prior for each arm. 
    Each element needs to be of class \samp{normMix} (See \samp{RBesT} package documentation).
  }
  \item{alpha}{
    Significance level for the frequentist multiple contrast test. 
    If no critical values are supplied via \samp{critV} this is used 
    to derive critical values for Bayesian decision rule.
}
  \item{alternative}{
    Character determining the alternative for the multiple contrast
    trend test. 
  }
  \item{na.action}{
    A function which indicates what should happen when
    the data contain NAs.
  }
  \item{mvtcontrol}{
    A list specifying additional control parameters for the \samp{qmvt}
    and \samp{pmvt} calls in the code used to obtain frequentist critical values, 
    see also \code{\link{mvtnorm.control}} for details.
  }
  \item{contMat}{
    Contrast matrix to apply to the posterior dose-response estimates. The
    contrasts need to be in the columns of the matrix (i.e. the column
    sums need to be 0).   
  }
  \item{critV}{
    Supply a critical value for the maximum posterior probability of the contrasts 
    being greater (or less) than zero that needs to be surpassed to establish a non-flat dose-response. If this argument is NULL, 
    this will be derived from critical values for frequentist MCP-Mod using the provided \samp{alpha}.
  }
}
\details{
  For the default calculation of optimal contrasts the prior information is ignored (i.e. contrasts are caulcated 
  in the same way as in \code{\link{MCTtest}}). 
  Fleischer et al. (2022) discuss using contrasts that take the prior effective sample sizes into account, which can be
  slightly more favourable for the Bayesian MCT test. Such alternative contrasts can be directly handed over via
  the \samp{contMat} argument.
  
  For analysis with covariate adjustment, covariate-adjusted \samp{resp} and \samp{S} can be supplied together
  with using \samp{type = "general"}. See the vignette
  "Design and analysis template MCP-Mod for binary data" for an example on how to obtain covariate adjusted estimates.
  Adequate care nee
}
\value{
  An object of class bMCTtest, a list containing the output.
}
\references{
  Fleischer, F., Bossert, S., Deng, Q., Loley, C. and Gierse, J. (2022).
  Bayesian MCP-Mod, \emph{Pharmaceutical Statistics}, \bold{21},
  654--670
}
\author{
  Marius Thomas
}
\seealso{
  \code{\link{MCTtest}}, \code{\link{optContr}}
}
\examples{

require(RBesT)

###############################
## Normal outcome
###############################

data(biom)
## define shapes for which to calculate optimal contrasts
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 1, 1, 1), doses = c(0, 0.05, 0.2, 0.6, 1))
## specify an informative prior for placebo, weakly informative for other arms
plc_prior <- mixnorm(inf = c(0.8, 0.4, 0.1), rob = c(0.2, 0.4, 10))
vague_prior <- mixnorm(c(1, 0, 10))
prior <- c(list(plc_prior), rep(list(vague_prior), length(unique(biom$dose)) - 1))

m1 <- bMCTtest(dose, resp, biom, models=modlist, prior = prior)
## now supply a critical value
m2 <- bMCTtest(dose, resp, biom, models=modlist, prior = prior, critV = 0.99)

####################################
## Binary outcome with covariates 
####################################

## generate data
logit <- function(p) log(p / (1 - p))
inv_logit <- function(y) 1 / (1 + exp(-y))
doses <- c(0, 0.5, 1.5, 2.5, 4)

## set seed and ensure reproducibility across R versions
set.seed(1, kind = "Mersenne-Twister", sample.kind = "Rejection", normal.kind = "Inversion")
group_size <- 100
dose_vector <- rep(doses, each = group_size)
N <- length(dose_vector)
## generate covariates
x1 <- rnorm(N, 0, 1)
x2 <- factor(sample(c("A", "B"), N, replace = TRUE, prob = c(0.6, 0.4)))
## assume approximately logit(10%) placebo and logit(35%) asymptotic response with ED50=0.5
prob <- inv_logit(emax(dose_vector, -2.2, 1.6, 0.5) + 0.3 * x1 + 0.3 * (x2 == "B"))
dat <- data.frame(y = rbinom(N, 1, prob),
                  dose = dose_vector, x1 = x1, x2 = x2)
                  
## specify an informative prior for placebo (on logit scale), weakly informative for other arms                  
plc_prior <- mixnorm(inf = c(0.8, -2, 0.5), rob = c(0.2, -2, 10))
vague_prior <- mixnorm(c(1, 0, 10))
prior <- c(list(plc_prior), rep(list(vague_prior), length(unique(biom$dose)) - 1))

## candidate models
mods <- Mods(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1),
             placEff = logit(0.1), maxEff = logit(0.35)-logit(0.1),
             doses = doses)

fit_cov <- glm(y~factor(dose) + 0 + x1 + x2, data = dat, family = binomial)

covariate_adjusted_estimates <- function(mu_hat, S_hat, formula_rhs, doses, other_covariates, n_sim) {
  ## predict every patient under *every* dose
  oc_rep <- as.data.frame(lapply(other_covariates, function(col) rep(col, times = length(doses))))
  d_rep <- rep(doses, each = nrow(other_covariates))
  pdat <- cbind(oc_rep, dose = d_rep)
  X <- model.matrix(formula_rhs, pdat)
  ## average on probability scale then backtransform to logit scale
  mu_star <- logit(tapply(inv_logit(X \%*\% mu_hat), pdat$dose, mean))
  ## estimate covariance matrix of mu_star
  pred <- replicate(n_sim, logit(tapply(inv_logit(X \%*\% drop(rmvnorm(1, mu_hat, S_hat))),
                                        pdat$dose, mean)))
  return(list(mu_star = as.numeric(mu_star), S_star = cov(t(pred))))
}

ca <- covariate_adjusted_estimates(coef(fit_cov), vcov(fit_cov), ~factor(dose)+0+x1+x2,
                                   doses, dat[, c("x1", "x2")], 1000)
bMCTtest(doses, ca$mu_star, S = ca$S_star, type = "general", models = mods, prior = prior)

################################################
## example with contrasts handed over
################################################

data(biom)
## define shapes for which to calculate optimal contrasts
modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
                linInt = c(0, 1, 1, 1), doses = c(0, 0.05, 0.2, 0.6, 1))
                
## specify an informative prior for placebo, weakly informative for other arms
plc_prior <- mixnorm(inf = c(0.8, 0.4, 0.1), rob = c(0.2, 0.4, 10), sigma = 0.7)
vague_prior <- mixnorm(c(1, 0, 10), sigma = 0.7)
prior <- c(list(plc_prior), rep(list(vague_prior), length(unique(biom$dose)) - 1))

## use prior effective sample sizes to calculate optimal contrasts
prior_ess <- unlist(lapply(prior, ess))
n_grp <- as.numeric(table(biom$dose))
weights <- n_grp + prior_ess
cmat <- optContr(modlist, w = weights)

bMCTtest(dose, resp, biom, models=modlist, prior = prior, contMat = cmat)
}



