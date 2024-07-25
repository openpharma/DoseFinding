## optimal designs for model-fitting

## calculate gradient of model and gradient of TD
calcGrads <- function(fmodels, doses, weights,
                      Delta, off, scal, direction,
                      designCrit){
  modgrad <- TDgrad <- nPar <- vector("list", modCount(fmodels, fullMod=TRUE))
  z <- 1
  for(nam in names(fmodels)){
    pars <- fmodels[[nam]]
    if(is.matrix(pars)){
      for(i in 1:nrow(pars)){
        modgrad[[z]] <- t(gradCalc(nam, pars[i,], doses, off=off, scal=scal)*sqrt(weights))
        if(designCrit != "Dopt")
          TDgrad[[z]] <- calcTDgrad(nam, pars[i,], Delta, direction, off, scal)
        nPar[[z]] <- nPars(nam)
        z <- z+1
      }
    } else {
      modgrad[[z]] <- t(gradCalc(nam, pars, doses, off=off, scal=scal)*sqrt(weights))
      if(designCrit != "Dopt")
        TDgrad[[z]] <- calcTDgrad(nam, pars, Delta, direction, off, scal)
      nPar[[z]] <- nPars(nam)
      z <- z+1        
    }
  }
  modgrads <- do.call("c", modgrad)
  TDgrad <- do.call("c", TDgrad)
  nPar <- do.call("c", nPar)

  list(modgrads=modgrads, TDgrad=TDgrad, nPar=nPar)
}


## returns the number of parameters (needed for C call)
nPars <- function(mods){
  builtIn <- c("linlog", "linear", "quadratic", 
               "emax", "exponential", "logistic", 
               "betaMod", "sigEmax")
  ind <- match(mods, builtIn)
  if(any(is.na(ind))){
    stop(mods[which(is.na(ind))], " model not allowed in optDesign")
  }
  c(2,2,3,3,3,4,4,4)[ind]
}

## function which calls different optimizers
callOptim <- function(func, method, nD, control, lowbnd, uppbnd){
  ## actual optimizer
  if(method == "nlminb"){ # nlminb and optim run on transformed values
    res <- nlminb(getStart(nD), objective=func, control = control,
                  lower=rep(0, nD), upper=rep(pi, nD))
  } else if(method == "Nelder-Mead"){
    res <- optim(getStart(nD), fn=func, control = control)
  } else if(method == "solnp"){ # no need for transformed values for solnp
    avail <- requireNamespace("Rsolnp", quietly = TRUE)
    if(!avail)
      stop("Need suggested package Rsolnp for this calculation to use solnp optimizer")
    ## get starting value (need feasible starting value for solnp)
    ## try whether equal allocation is feasible
    eq <- rep(1/nD, nD)
    if(all(eq > lowbnd+0.001) & all(eq < uppbnd-0.001)){
      strt <- eq
    } else {
      slb <- sum(lowbnd)
      sub <- sum(uppbnd)
      gam <- (1-slb)/(sub-slb)
      strt <- lowbnd+gam*(uppbnd-lowbnd)
    }
    eqfun <- function(x, ...){
      sum(x)
    }
    con <- list(trace = 0)
    con[(namc <- names(control))] <- control
    res <- Rsolnp::solnp(strt, fun=func, eqfun=eqfun, eqB=1,
                         control = con, LB = lowbnd, UB = uppbnd)
  } 
  res
}

## transforms from unconstrained values R^k into constrained
## values in S^k = {w|sum_i w_i=1 and w_i >= 0}
transTrig <- function(y, k){
  a <- numeric(k)  
  if(k == 2){
    a[1] <- sin(y[1])^2
  } else {
    a[1:(k-1)] <- sin(y)^2
    a[2:(k-1)] <- a[2:(k-1)]*cumprod(cos(y[1:(k-2)])^2)
  }
  a[k] <- prod(cos(y[1:(k-1)])^2)
  a
}

## identity function
idtrans <- function(y, k){
  y
}

## calculate uniform design but on R^k scale
## (inverse of transTrig at uniform design)
getStart <- function(k){
  y <- numeric(k-1)
  eq <- 1/k
  y[1] <- asin(sqrt(eq))
  for(j in 2:(k-1)){
    y[j] <- asin(sqrt(eq/prod(cos(y[(1:j)-1])^2)))
  }
  y
}

## function called in the optimization (design criterion is
## implemented in C and called "critfunc")
optFunc <- function(x, xvec, pvec, nD, probs, M, n, nold, bvec, designCrit,
                    trans, standInt){
  xtrans <- do.call("trans", list(x, nD))
  res <- .C("critfunc", xvec, pvec, nD, probs, M, xtrans, n,
            nold, double(16), as.double(1e-15), bvec, designCrit, standInt,
            double(1), PACKAGE = "DoseFinding")
  res[[14]]
}

## auxiliary function for efficient rounding
which.is.max <- function (x){
    y <- seq_along(x)[x == max(x)]
    if (length(y) > 1L) 
        sample(y, 1L)
    else y
}

getCompositions <- function(N, M){
  nC <- choose(N+M-1, M-1)
  lst <- .C("getcomp", comp=integer(nC*M), integer(M-1),
            as.integer(N), as.integer(M-1), as.integer(nC),
            PACKAGE = "DoseFinding")
  matrix(lst$comp, byrow = TRUE, nrow = nC)
}


## calculate all possible compositions of n patients to nDoses groups
## (assuming a certain block-size) upper and lower bounds on the
## allocations can also be specified
getDesMat <- function(n, nDoses, lowbnd = rep(0, nDoses), 
                      uppbnd = rep(1, nDoses), groupSize,
                      maxvls1, maxvls2){
  if(n %% groupSize)
    stop("n needs to be divisible by groupSize")
  nG <- n/groupSize
  combn <- choose(nG+nDoses-1,nDoses-1)
  if(combn > maxvls1)
    stop(combn, " (unrestricted) combinations, increase maxvls1 in control 
         argument if this calculation should be performed")

  desmat <- getCompositions(nG, nDoses)/nG
 
  if(any(lowbnd > 0) | any(uppbnd < 1)){
    comp <- matrix(lowbnd, byrow = TRUE, ncol = nDoses, nrow=nrow(desmat))
    LindMat <- desmat >= comp
    comp <- matrix(uppbnd, byrow=TRUE, ncol = nDoses, nrow=nrow(desmat))
    UindMat <- desmat <= comp
    ind <- rowSums(LindMat*UindMat) == nDoses
    desmat <- desmat[ind,]
    if(nrow(desmat) == 0)
      stop("no design is compatible with bounds specified in lowbnd and uppbnd")
  }
  if(nrow(desmat) > maxvls2)
    stop(nrow(desmat), " combinations, increase maxvls2 in control argument if
         this calculation should be performed")
  desmat
}
