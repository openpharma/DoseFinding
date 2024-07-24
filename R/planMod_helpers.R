## various functions for assessing the operating characteristics of a design
## for model-based estimation of dose-response functions

## calculates the variance of the estimated curve
getPredVar <- function(model, cf, V,  pDose, off, scal){
  gr <- gradCalc(model, cf, pDose, off, scal)
  gr0 <- gradCalc(model, cf, 0, off, scal)
  grd <- sweep(gr, 2, gr0)
  out <- apply(grd, 1, function(x){
    as.numeric(t(x)%*%V%*%x)
  })
  out
}

## calculates the variance of the EDp estimate
getEDVar <- function(model, cf, V, scale = c("unrestricted", "logit"),
                     p, maxD, off, scal, nodes){
  grd <- calcEDgrad(model, cf, maxD, p, off, scal, nodes)
  if(scale == "logit"){
    tmp <- calcED(model, cf, p, maxD, "continuous",
                                off=off, scal=scal, nodes=nodes)
    grd <- grd*(-maxD/(tmp*(tmp-maxD)))
  }
  grd <- as.numeric(grd)
  return(as.numeric(t(grd)%*%V%*%grd))
}

## calculates the variance of the TD estimate
getTDVar <- function(model, cf, V, scale = c("unrestricted", "log"),
                     Delta, direction = c("increasing", "decreasing"), off, scal, nodes){
  tmp <- calcTD(model, cf, Delta, "continuous", 
                              direction, off=off, scal=scal,
                              nodes = nodes)
  grd <- calcTDgrad(model, cf, Delta, direction, off, scal, nodes)
  if(scale == "log")
    grd <- grd/tmp
  grd <- as.numeric(grd)
  return(as.numeric(t(grd)%*%V%*%grd))
}

## calculates approximate covariance matrix for parameter estimates
aprCov <- function(doses, model, cf, S, off, scal){
  F <- gradCalc(model, cf, doses, off, scal)
  V <- try(solve(t(F)%*%solve(S)%*%F))
  if(inherits(V, "try-error")){
    warning("Could not calculate covariance matrix; Fisher information singular.")
    return(NA)
  }
  V
}


tableMatch <- function(x, match){
  ## like "table", but also returns categories with 0 counts
  out <- numeric(length(match))
  for(i in 1:length(match)){
    out[i] <- sum(x == match[i], na.rm=TRUE)
  }
  names(out) <- match
  out
}

## calculate the predictions for the fitted models
getSimEst <- function(x, type = c("dose-response", "ED", "TD"),
                      doseSeq, direction, p, Delta, placAdj = FALSE){
  modelSel <- attr(x$sim, "modelSel") 
  model <- attr(x, "model")
  coefs <- attr(x$sim, "coefs")
  off <- attr(x, "off")
  scal <- attr(x, "scal")
  nSim <- attr(x$sim, "nSim")
  altModels <- attr(x, "altModels")
  nAlt <- modCount(altModels, fullMod=TRUE)
  doses <- attr(x, "doses")
  maxD <- max(doses)
  type <- match.arg(type)
  if(type == "TD"){
    if(missing(direction))
      stop("need direction for TD calculation")
    if(Delta <= 0)
      stop("\"Delta\" needs to be > 0")
  }
  out <- vector("list", nAlt)
  for(i in 1:nAlt){
    ind <- matrix(ncol = length(model), nrow = nSim)
    if(type == "dose-response"){
      resMat <- matrix(nrow = nSim, ncol = length(doseSeq))
      colnames(resMat) <- doseSeq
      rownames(resMat) <- 1:nSim
      for(j in 1:length(model)){
        ind[,j] <- modelSel[,i] == model[j]
        if(any(ind[,j])){
          cf <- do.call("rbind", (coefs[[i]])[ind[,j]])
          resMat[ind[,j]] <- predSamples(samples=cf, placAdjfullPars = TRUE,
                                         doseSeq=doseSeq,
                                         placAdj=placAdj, model=model[j],
                                         scal=scal, off=off, nodes = NULL)
        }
        out[[i]] <- resMat
      }
    }
    if(is.element(type, c("TD", "ED"))){
      resVec <- numeric(nSim)
      for(j in 1:length(model)){
        ind[,j] <- modelSel[,i] == model[j]
        if(any(ind[,j])){
          cf <- do.call("rbind", (coefs[[i]])[ind[,j]])
          if(type == "TD"){
            resVec[ind[,j]] <- apply(cf, 1, function(z){
              calcTD(model[j], z, Delta, "continuous", direction,
                     off=off, scal=scal)
            })
          }
          if(type == "ED"){
            resVec[ind[,j]] <- apply(cf, 1, function(z){
              calcED(model[j], z, p, maxD, "continuous", off=off, scal=scal)
            })
          }
        } 
      }
      out[[i]] <- resVec
    }
  }
  names(out) <- colnames(getResp(attr(x, "altModels"), doses=0)) ## horrible hack need to improve!
  out
}
  

plotDoseSims <- function(x, type = c("ED", "TD"), p, Delta, xlab){
  altMods <- attr(x, "altModels")
  direction <- attr(altMods, "direction")
  if(type == "ED"){
    out <- getSimEst(x, "ED", p=p)
    trueDoses <- ED(altMods, p=p, EDtype="continuous")
  } else {
    out <- getSimEst(x, "TD", Delta=Delta, direction=direction)
    trueDoses <- TD(altMods, Delta=Delta, TDtype="continuous",
                    direction=direction)
  }
  ## write plotting data frame
  nams <- names(out)
  group <- factor(rep(1:length(nams), each=length(out[[1]])), labels=nams)
  pdat <- data.frame(est = do.call("c", out),
                     group = group)
  ## determine limits for x-axis
  rngQ <- tapply(pdat$est, pdat$group, function(x){
    quantile(x, c(0.025, 0.975), na.rm=TRUE)
  })
  rngQ <- do.call("rbind", rngQ)
  rng <- c(min(rngQ[,1], na.rm = TRUE), max(rngQ[,2], na.rm = TRUE))
  delt <- diff(rng)*0.04
  ## truncate x-axis to 2*maxdose
  maxdose <- max(attr(x, "doses"))
  xlimU <- min(2*maxdose, max(rng[2], maxdose)+delt)
  xlimL <- max(-0.05*maxdose, min(0, rng[1])-delt)
  xlim <- c(xlimL, xlimU)
  parVal <- ifelse(type == "ED", paste("p=", p, sep=""), paste("Delta=", Delta, sep=""))
  maintxt <- paste("95%, 80%, 50% intervals and median of simulated ", type,
                   " estimates (", parVal, ")", sep = "")
  key <- list(text = list(maintxt, cex = 0.9))
  lattice::bwplot(~est|group, data=pdat, xlab = xlab, trueDoses=trueDoses,
                  xlim = xlim,
                  panel = function(...){
                    z <- lattice::panel.number()
                    lattice::panel.grid(v=-1, h=0, lty=2, col = "lightgrey")
                    lattice::panel.abline(v=trueDoses[z], col = "red", lwd=2)
                    lattice::panel.abline(v=c(0, max(attr(x, "doses"))), col = "grey", lwd=2)
                    probs <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
                    simDoseEst <- list(...)$x
                    quants <- quantile(simDoseEst, probs, na.rm = TRUE)
                    lattice::llines(c(quants[1], quants[7]), c(1,1), lwd=2, col=1)
                    lattice::llines(c(quants[2], quants[6]), c(1,1), lwd=5, col=1)
                    lattice::llines(c(quants[3], quants[5]), c(1,1), lwd=10, col=1)
                    lattice::lpoints(quants[4], 1, cex=2, pch="|", col=1)
                    if(type == "TD")
                      lattice::ltext(xlim[2], 1.5, pos = 2, cex = 0.75,
                                     labels = paste("% No TD:", mean(is.na(simDoseEst))*100, "%"))
                  }, layout = c(1,length(out)), as.table = TRUE, key = key)
}

plotDRSims <- function(x, placAdj = FALSE, xlab, ylab){
  altMods <- attr(x, "altModels")
  rng <- range(attr(x, "doses"))
  doseSeq <- seq(rng[1], rng[2], length = 51)
  out <- getSimEst(x, type = "dose-response", doseSeq=doseSeq, placAdj = placAdj)
  trueMn <- getResp(altMods, doses=doseSeq)
  if(placAdj){
    trueMn <- trueMn-trueMn[1,]
  }
  nM <- length(out)
  resp <- vector("list", length=nM)
  for(i in 1:nM){
    qMat <-apply(out[[i]], 2, function(y){
      quantile(y, c(0.025, 0.25, 0.5, 0.75, 0.975))
    })
    resp[[i]] <- c(t(qMat))
  }
  
  resp <- do.call("c", resp)
  quant <- rep(rep(c(0.025, 0.25, 0.5, 0.75, 0.975), each = 51), nM)
  dose <- rep(doseSeq, nM*5)
  model <- factor(rep(1:nM, each = 5*51), labels = names(out))
  key <- list(text =
                list("Pointwise 95%, 50% intervals and median of simulated dose-response estimates", cex = 0.9))
  
  lattice::xyplot(resp~dose|model, groups = quant, xlab=xlab, ylab = ylab,
                  panel = function(...){
                    ## plot grid
                    lattice::panel.grid(v=-1, h=-1, col = "lightgrey", lty=2)
                    ## plot estimates
                    panel.dat <- list(...)
                    ind <- panel.dat$subscripts
                    LB95.x <- panel.dat$x[panel.dat$groups[ind] == 0.025]
                    LB95 <- panel.dat$y[panel.dat$groups[ind] == 0.025]
                    UB95.x <- panel.dat$x[panel.dat$groups[ind] == 0.975]
                    UB95 <- panel.dat$y[panel.dat$groups[ind] == 0.975]
                    lattice::lpolygon(c(LB95.x, rev(UB95.x)), c(LB95, rev(UB95)),
                             col = "lightgrey", border = "lightgrey")
                    LB50.x <- panel.dat$x[panel.dat$groups[ind] == 0.25]
                    LB50 <- panel.dat$y[panel.dat$groups[ind] == 0.25]
                    UB50.x <- panel.dat$x[panel.dat$groups[ind] == 0.75]
                    UB50 <- panel.dat$y[panel.dat$groups[ind] == 0.75]
                    lattice::lpolygon(c(LB50.x, rev(UB50.x)), c(LB50, rev(UB50)),
                                      col = "darkgrey", border = "darkgrey")
                    MED.x <- panel.dat$x[panel.dat$groups[ind] == 0.5]
                    MED <- panel.dat$y[panel.dat$groups[ind] == 0.5]
                    lattice::llines(MED.x, MED, col = 1,lwd = 1.5)
                    ## plot true curve
                    z <- lattice::panel.number()
                    lattice::llines(doseSeq, trueMn[,z], col=2, lwd=1.5)
                  }, as.table = TRUE, key=key)
}
