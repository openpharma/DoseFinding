########################################################################
#### Testing function to generate doses and sample size allocs.
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
