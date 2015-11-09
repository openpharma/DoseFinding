## calculations of equi-coordinate quantiles for pmvt

wgts <- function(typred, ttarg){
  dists <- abs(typred-ttarg)
  if(sum(dists < 1 & dists > 0.0001) < 3){
    return(rep(1, length(typred)))
  } 
  (1-dists^3)^3*(dists < 1)
}

stop_crit <- function(xtol=NULL, ytol=NULL, fit, xest, ttarg, level=0.05){
  if(!is.null(xtol)){
    xvals <- c(xest-xtol, xest+xtol)
    pred <- predict(fit, data.frame(x=xvals), se.fit=TRUE)
    crit <- qt(1-level/2, df=fit$df.residual)
    cond1 <- pred$fit[1]+crit*pred$se.fit[1] < ttarg
    cond2 <- pred$fit[2]-crit*pred$se.fit[2] > ttarg
    return(list(stop=cond1 & cond2))
  }
  if(!is.null(ytol)){
    pred <- predict(fit, data.frame(x=xest), se.fit=TRUE)
    crit <- qt(1-level/2, df=fit$df.residual)
    cond1 <- pred$fit[1]+crit*pred$se.fit[1] < ttarg+ytol
    cond2 <- pred$fit[1]-crit*pred$se.fit[1] > ttarg-ytol
    return(list(stop=cond1 & cond2))
  }
}

get_new_points <- function(fit, xest, targ, level=0.05){
  pred <- predict(fit, data.frame(x=xest), se.fit=TRUE)
  cf <- coef(fit)
  xLB <- (pred$fit-pred$se.fit-cf[1])/cf[2]
  xUB <- (pred$fit+pred$se.fit-cf[1])/cf[2]
  c(xLB, xUB)
}

get_est <- function(fit, ttarg){
  cf <- coef(fit)
  as.numeric(ttarg-cf[1])/cf[2]
}

fitlr <- function(x, y, weights, method, N, sterr){
  ind <- weights > 0
  weights <- weights[ind]
  x <- x[ind]  
  if(method == "lm"){
    y <- y[ind]
    fit <- glm(y~x, weights = weights,
               family=quasi("probit"))
  }
  if(method == "wlm"){
    y <- y[ind]
    sterr <- sterr[ind]
    sterr <- sterr/min(sterr)
    wgts <- weights/sterr^2
    fit <- glm(y~x, weights = wgts,
               family=quasi("probit"))
  }
  if(method == "glm"){
    ymat <- cbind(y, 1-y)*N
    ymat <- ymat[ind,]
    fit <- glm(ymat~x, weights = weights,
               family=binomial("probit"))
  }
  return(fit)
}

get_quant_loclin <- function(func, targ, interval, 
                             ytol=NULL, xtol=NULL, maxiter = 500, 
                             method = c("lm", "wlm", "glm"),
                             verbose = FALSE, ...){
  ## input argument checks
  if(interval[2] <= interval[1])
    stop("first entry of interval needs to be smaller than second")
  method <- match.arg(method)
  ## set up output vectors
  x <- y <- sterr <- numeric(4+maxiter)
  xest <- numeric(1+maxiter)
  ## first three evaluations
  x[1:3] <- c(interval, mean(interval))
  res <- lapply(x[1:3], function(x) func(x, ...))
  y[1:3] <- sapply(res, function(x) x)
  if(method == "wlm"){
    if(is.null(attr(res[[1]], "sterr")))
      stop("For method == wlm, func neets to have an attribute sterr")
    sterr[1:3] <- sapply(res, function(x) attr(x, "sterr"))
  }
  if(method == "glm"){
    if(is.null(attr(res[[1]], "n")))
      stop("For method == glm, func neets to have an attribute n")
    N <- attr(res[[1]], "n")
  }
  ## check for non-monotone function
  if(y[1] > y[2])
    stop("func does not appear to be monotone")
  ttarg <- qnorm(targ) # transform target as well
  
  ## start iterating
  ycur <- y[1:3];xcur <- x[1:3]
  fit <- fitlr(xcur, ycur, weights=rep(1,3),
               method = method, N=N, sterr=sterr[1:3])
  xest[1] <- get_est(fit, ttarg)
  x[4] <- xest[1]
  res <- func(x[4], ...)
  y[4] <- res
  if(method == "wlm")
    sterr[4] <- attr(res, "sterr")
  count <- 4

  for(i in 1:maxiter){
    ind <- 1:count
    xcur <- x[ind];ycur <- y[ind]
    tycurpred <- predict(fit, newdata=data.frame(x=xcur))
    weights <- wgts(tycurpred, ttarg)
    fit <- fitlr(xcur, ycur, weights=weights,
                 method = method, N=N, sterr=sterr[ind])
    stp <- stop_crit(xtol, ytol, fit, xest[i], ttarg)
    if(verbose & i > 1){
      txt <- sprintf("Iteration %i, Estimate: %f, LB: %f (%f), UB: %f (%f)\n",
                     i, xest[i],
                     x[count-1], y[count-1],
                     x[count], y[count])
      cat(txt)
    }
    if(stp$stop)
      break
    xest[i+1] <- get_est(fit, ttarg)
    new_x <- get_new_points(fit, xest[i+1], targ)
    new_y <- lapply(new_x, function(x) func(x, ...))
    ind <- (count+1):(count+2)
    x[ind] <- new_x
    y[ind] <- c(new_y[[1]], new_y[[2]])
    if(method == "wlm")
      sterr[ind] <- sapply(new_y, function(x) attr(x, "sterr"))
    count <- count+2
  }
  out <- xest[i]
  if(i == maxiter){
    attr(out, "message") <- "Maximum number of iterations reached without sufficient accuracy"
  } else {
    attr(out, "message") <- "Normal Completion"
  }
  attr(out, "iterations") <- list(xest=c(rep(NA,4), xest[1:i]),
                                  x=x[1:count],
                                  y=ycur[1:count])
  out
}

########################################################################
## non-exported objects from mvtnorm (only temporarily needed)
isInf <- function (x) 
  x > 0 & is.infinite(x)
dots2GenzBretz <- function (...) 
{
    addargs <- list(...)
    fm1 <- sapply(names(addargs), function(x) length(grep(x, 
        names(formals(GenzBretz)))) == 1)
    fm2 <- sapply(names(addargs), function(x) length(grep(x, 
        names(formals(uniroot)))) == 1)
    algorithm <- NULL
    uniroot <- NULL
    if (any(fm1)) 
        algorithm <- do.call("GenzBretz", addargs[fm1])
    if (any(fm2)) 
        uniroot <- addargs[fm2]
    list(algorithm = algorithm, uniroot = uniroot)
}
checkmvArgs <- function (lower, upper, mean, corr, sigma) 
{
    UNI <- FALSE
    if (!is.numeric(lower) || !is.vector(lower)) 
        stop(sQuote("lower"), " is not a numeric vector")
    if (!is.numeric(upper) || !is.vector(upper)) 
        stop(sQuote("upper"), " is not a numeric vector")
    if (!is.numeric(mean) || !is.vector(mean)) 
        stop(sQuote("mean"), " is not a numeric vector")
    if (is.null(lower) || any(is.na(lower))) 
        stop(sQuote("lower"), " not specified or contains NA")
    if (is.null(upper) || any(is.na(upper))) 
        stop(sQuote("upper"), " not specified or contains NA")
    rec <- cbind(lower, upper, mean)
    lower <- rec[, "lower"]
    upper <- rec[, "upper"]
    if (!all(lower <= upper)) 
        stop("at least one element of ", sQuote("lower"), " is larger than ", 
            sQuote("upper"))
    mean <- rec[, "mean"]
    if (any(is.na(mean))) 
        stop("mean contains NA")
    if (is.null(corr) && is.null(sigma)) {
        corr <- diag(length(lower))
    }
    if (!is.null(corr) && !is.null(sigma)) {
        sigma <- NULL
        warning("both ", sQuote("corr"), " and ", sQuote("sigma"), 
            " specified: ignoring ", sQuote("sigma"))
    }
    if (!is.null(corr)) {
        if (!is.numeric(corr)) 
            stop(sQuote("corr"), " is not numeric")
        if (!is.matrix(corr)) {
            if (length(corr) == 1) 
                UNI <- TRUE
            if (length(corr) != length(lower)) 
                stop(sQuote("diag(corr)"), " and ", sQuote("lower"), 
                  " are of different length")
        }
        else {
            if (length(corr) == 1) {
                UNI <- TRUE
                corr <- corr[1, 1]
                if (length(lower) != 1) 
                  stop(sQuote("corr"), " and ", sQuote("lower"), 
                    " are of different length")
            }
            else {
                if (length(diag(corr)) != length(lower)) 
                  stop(sQuote("diag(corr)"), " and ", sQuote("lower"), 
                    " are of different length")
                if (!chkcorr(corr)) 
                  stop(sQuote("corr"), " is not a correlation matrix")
            }
        }
    }
    if (!is.null(sigma)) {
        if (!is.numeric(sigma)) 
            stop(sQuote("sigma"), " is not numeric")
        if (!is.matrix(sigma)) {
            if (length(sigma) == 1) 
                UNI <- TRUE
            if (length(sigma) != length(lower)) 
                stop(sQuote("diag(sigma)"), " and ", sQuote("lower"), 
                  " are of different length")
        }
        else {
            if (length(sigma) == 1) {
                UNI <- TRUE
                sigma <- sigma[1, 1]
                if (length(lower) != 1) 
                  stop(sQuote("sigma"), " and ", sQuote("lower"), 
                    " are of different length")
            }
            else {
                if (length(diag(sigma)) != length(lower)) 
                  stop(sQuote("diag(sigma)"), " and ", sQuote("lower"), 
                    " are of different length")
                if (!isTRUE(all.equal(sigma, t(sigma))) || any(diag(sigma) < 
                  0)) 
                  stop(sQuote("sigma"), " is not a covariance matrix")
            }
        }
    }
    list(lower = lower, upper = upper, mean = mean, corr = corr, 
        sigma = sigma, uni = UNI)
}
chkcorr <- function (x) 
{
    if (!is.matrix(x)) 
        return(FALSE)
    rownames(x) <- colnames(x) <- NULL
    storage.mode(x) <- "numeric"
    ONE <- 1 + sqrt(.Machine$double.eps)
    ret <- (min(x) < -ONE || max(x) > ONE) || !isTRUE(all.equal(diag(x), 
        rep(1, nrow(x))))
    !ret
}
########################################################################

## main function
qmvtDF <- function (p, tail = c("lower.tail", "upper.tail", "both.tails"),
                    df = 1, delta = 0, corr = NULL, sigma = NULL, 
                    algorithm = GenzBretz(), type = c("Kshirsagar", "shifted"),
                    ytol = 0.001, maxiter = 500, ...) {
  if (length(p) != 1 || (p <= 0 || p >= 1)) 
    stop(sQuote("p"), " is not a double between zero and one")
  dots <- dots2GenzBretz(...) 
  if (!is.null(dots$algorithm) && !is.null(algorithm)) 
    algorithm <- dots$algorithm
  type <- match.arg(type)
  tail <- match.arg(tail)
  if (tail == "both.tails" && p < 0.5) 
    stop("cannot compute two-sided quantile for p < 0.5")
  dim <- 1
  if (!is.null(corr)) 
    dim <- NROW(corr)
  if (!is.null(sigma)) 
    dim <- NROW(sigma)
  lower <- upper <- rep.int(0, dim)
  args <- checkmvArgs(lower, upper, delta, corr, sigma)
  if (args$uni) {
    if (tail == "both.tails") 
      p <- ifelse(p < 0.5, p/2, 1 - (1 - p)/2)
    if (df == 0 || isInf(df)) {
      q <- qnorm(p, mean = args$mean, lower.tail = (tail != 
                      "upper.tail"))
    } else {
      q <- qt(p, df = df, ncp = args$mean, lower.tail = (tail != 
                            "upper.tail"))
    }
    qroot <- list(quantile = q, f.quantile = 0)
    return(qroot)
  }
  dim <- length(args$mean)
  pfct <- function(q) {
    switch(tail, both.tails = {
      low <- rep(-abs(q), dim)
      upp <- rep(abs(q), dim)
    }, upper.tail = {
      low <- rep(q, dim)
      upp <- rep(Inf, dim)
    }, lower.tail = {
      low <- rep(-Inf, dim)
      upp <- rep(q, dim)
    }, )
    ret <- pmvt(lower = low, upper = upp, df = df, delta = args$mean, 
                corr = args$corr, sigma = args$sigma, algorithm = algorithm, 
                type = type)
    return(ret)
  }
  intp <- switch(tail,
                 both.tails = (1 - (1 - p)/2)^(c(1,1/dim)),
                 upper.tail = 1 - p^c(1,1/dim),
                 lower.tail = p^c(1,1/dim))
  if(is.finite(df) && (df > 0)){
    intq <- qt(intp, df = df)
  } else {
    intq <- qnorm(intp)
  }
  qroot <- get_quant_loclin(pfct, p, intq*c(0.8, 1.25),
                            ytol=ytol, maxiter=maxiter,
                            verbose=FALSE, method="lm")
  attr(qroot, "iterations") <- NULL
  qroot
}
