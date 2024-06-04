## helper functions for calculating optimal contrasts and critical value

optC <- function(mu, Sinv = NULL, placAdj = FALSE){
  ## calculate optimal contrast for given mu and Sinv (Sinv = proportional to inv covariance matrix)
  if(!placAdj){
    aux <- rowSums(Sinv)  # Sinv %*% 1
    mn <- sum(mu * aux)/sum(aux) # formula is: S^(-1)(mu-mu*S^(-1)*1/(1*S^(-1)1)1)
    val <- Sinv %*% (mu - mn)
    ## now center so that sum is 0
    ## and standardize to have norm 1
    val <- val - sum(val)
  } else { # placAdj = TRUE
    val <- Sinv %*% mu     
  }
    val/sqrt(sum(val^2))
}

constOptC <- function(mu, Sinv = NULL, placAdj = FALSE, direction){
  ## calculate optimal contrasts under the additional constraint that
  ## the control and the active treatment groups have a different sign
  ## in the contrast
  S <- solve(Sinv) # ugly fix, we should use S as argument
  if(!placAdj){
    k <- length(mu)
    CC <- cbind(-1,diag(k-1))
    SPa <- CC%*%S%*%t(CC)
    muPa <- as.numeric(CC%*%mu)
  } else {
    k <- length(mu)+1
    SPa <- S
    muPa <- mu
  }
  ## determine direction of effect
  unContr <- solve(SPa)%*%muPa # unconstrained optimal contrast
  mult <- ifelse(direction == "increasing", 1, -1) # 1 increasing, -1 decreasing
  ## prepare call of quadprog::solve.QP
  D <- SPa
  d <- rep(0,k-1)
  tA <- rbind(muPa, 
              mult*diag(k-1))
  A <- t(tA)
  bvec <- c(1,rep(0,k-1))
  contr <- quadprog::solve.QP(D, d, A, bvec, meq=1)$solution
  contr[abs(contr) < 1e-10] <- 0
  if(!placAdj)
    contr <- c(-sum(contr), contr)
  contr/sqrt(sum(contr^2))
}


modContr <- function(means, W = NULL, Sinv = NULL, placAdj = FALSE,
                     type, direction){
  ## call optC on matrix
  ## check whether constant shape was specified and remove (can happen for linInt model)
  if(!placAdj){ 
    ind <- apply(means, 2, function(x){
      length(unique(x)) > 1 
    })
  } else { ## placAdj
    ind <- apply(means, 2, function(x){
      any(x != 0) 
    })
  }
  if(all(!ind))
    stop("All models correspond to a constant shape, no optimal contrasts calculated.")
  if(any(!ind)){
    nam <- colnames(means)[!ind]
    namsC <- paste(nam, collapse = ", ")
    if(length(nam) == 1){
      message("The ", namsC, " model has a constant shape, cannot
calculate optimal contrasts for this shape.")
    } else {
      message("The ", namsC, " models have a constant shape, cannot
calculate optimal contrasts for these shapes.")
    }
    means <- means[,ind, drop=FALSE]
  }

  if(is.null(Sinv))
    Sinv <- solve(W)
  if(type == "unconstrained"){
    out <- apply(means, 2, optC, Sinv = Sinv, placAdj = placAdj)
  } else { # type == "constrained"
    out <- apply(means, 2, constOptC, Sinv = Sinv,
                 placAdj = placAdj, direction = direction)
  }
  if(!is.matrix(out)){ ## can happen for placAdj=T and only 1 act dose
    nam <- names(out)
    out <- matrix(out, nrow = 1)
    colnames(out) <- nam
  }
  out
}


