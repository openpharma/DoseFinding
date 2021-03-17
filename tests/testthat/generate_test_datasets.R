# Functions for generating the datasets used in testing

# TODO: unify this mess

# from testsFitting.R ----------------------------------------------------------

genDFdats <- function(model, argsMod, doses, n, sigma, mu = NULL){
  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) n <- rep(n, nD)
  dose <- rep(dose,  n)
  args <- c(list(dose), argsMod)
  mu <- do.call(model, args)
  data.frame(dose = dose, resp = mu + rnorm(sum(n), sd = sigma))
}

getDosSampSiz <- function(){
  # generate dose levels
  mD <- runif(1, 0, 1500)
  nD <- max(rpois(1, 5), 4)
  p <- rgamma(nD, 3)
  p <- cumsum(p/sum(p))
  doses <- signif(c(0, mD*p), 3)

  # sample size allocations
  totSS <- rpois(1, rexp(1, 1/250))
  totSS <- max(totSS, 50)
  p <- rgamma(nD+1, 3);p <- p/sum(p)
  n <- round(p*totSS)
  n[n==0] <- rpois(sum(n==0), 1)+1
  list(doses=doses, n=n)
}

getDFdataSet <- function(doses, n){
  if(missing(doses) & missing(n)){
    ll <- getDosSampSiz()
  } else {
    ll <- list(doses = doses, n=n)
  }

  e0 <- rnorm(1, 0, 10)
  eMax <- rgamma(1, abs(e0)*0.5, 0.5)
  sig <- eMax/runif(1, 0.5, 5)
  if(runif(1)<0.3){
    aa <- genDFdats("betaMod", c(e0 = e0, eMax = eMax, delta1=runif(1, 0.5, 4),
                delta2=runif(1, 0.5, 4), scal=1.2*max(ll$doses)),
                ll$doses, ll$n, sig)
  } else {
    aa <- genDFdats("sigEmax", c(e0 = e0, eMax = eMax,
                                 ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                                 h=runif(1, 0.5, 4)), ll$doses, ll$n, sig)
  }
  N <- sum(ll$n)
  center <- c("blue", "green", "red", "yellow", "silver")
  aa <- data.frame(x= aa$dose, y=aa$resp, center=as.factor(sample(center, N, replace = T)),
                   age=runif(N, 1, 100))
  aa[sample(1:nrow(aa)),]
}

# from testsMCT.R ---------------------------------------------------------------

getDFdataSet_testsMCT <- function(doses, n){
  ll <- getDosSampSiz()
  e0 <- rnorm(1, 0, 10)
  eMax <- rgamma(1, abs(e0)*0.5, 0.5)*I(runif(1)<0.25)
  if(eMax > 0){ sig <- eMax/runif(1, 0.5, 5)}
  else { sig <- rgamma(1, abs(e0)*0.5, 0.5) }
  dosVec <- rep(ll$doses, ll$n)
  if(runif(1)<0.3){
    mnVec <- betaMod(dosVec, e0=e0, eMax=eMax, delta1=runif(1, 0.5, 5),
                     delta2=runif(1, 0.5, 5), scal=1.2*max(ll$doses))
  } else {
    mnVec <- logistic(dosVec, e0 = e0, eMax = eMax,
                      ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                      delta=runif(1, 0.5, max(ll$doses)/2))
  }
  resp <- rnorm(sum(ll$n), mnVec, sig)
  N <- sum(ll$n)
  cov1 <- as.factor(rpois(N, 5))
  cov2 <- runif(N, 1, 100)
  aa <- data.frame(x= dosVec, y=resp, cov1=cov1, cov2=cov2)
  aa[sample(1:nrow(aa)),]
}

getDFdataSet.bin <- function(doses, n){
  ll <- getDosSampSiz()
  ll$n <- ll$n+10
  e0 <- rnorm(1, 0, sqrt(3.28))
  eMax <- rnorm(1, 0, 5)
  dosVec <- rep(ll$doses, ll$n)
  if(runif(1)<0.3){
    mn <- betaMod(dosVec, e0 = e0, eMax = eMax, delta1=runif(1, 0.5, 5),
                  delta2=runif(1, 0.5, 5), scal=1.2*max(ll$doses))
  } else {
    mn <- logistic(dosVec, e0 = e0,
                   eMax = eMax, ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                   delta=runif(1, 0.5, max(ll$doses)/2))
  }
  resp <- rbinom(length(ll$n), ll$n, 1/(1+exp(-mn)))
  aa <- data.frame(dose = ll$doses, resp = resp)
  aa <- data.frame(x= aa$dose, y=aa$resp/ll$n, n=ll$n)
  aa[sample(1:nrow(aa)),]
}
