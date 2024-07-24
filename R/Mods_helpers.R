## functions related to creating, plotting candidate model sets

modCount <- function(models, fullMod = FALSE){
  ## counts the number of models in a candidate model-list
  if(!fullMod){
    nr <- lapply(names(models), function(x){
      xx <- models[[x]]
      if(is.null(xx))
        return(1)
      if(is.element(x, c("emax", "quadratic", "exponential")))
        return(length(xx))
      if(is.element(x, c("sigEmax", "logistic", "betaMod")))
        return(length(xx)/2)
      if(x == "linInt"){
        if(is.vector(xx))
          return(1)
        if(is.matrix(xx))
          return(nrow(xx))
      }
    })
  } else {
    nr <- lapply(models, function(x){
      if(is.vector(x))
        return(1)
      if(is.matrix(x))
        return(nrow(x))
    })
  }
  Reduce("+",nr)
}

getAddArgs <- function(addArgs, doses = NULL){
  if(!is.null(doses)){
    addArgs0 <- list(scal = 1.2*max(doses), off = 0.01*max(doses))
  } else {
    addArgs0 <- list(scal = NULL, off = NULL)
  }
  if(!is.null(addArgs)){
    if(!is.list(addArgs))
      stop("addArgs needs to be of class list")
    namA <- names(addArgs)
    if(!all(namA %in% c("scal", "off")))
      stop("addArgs need to have entries named scal and/or off")
    addArgs0[namA] <- addArgs
    if(length(addArgs0$scal) > 1 | length(addArgs0$off) > 1)
      stop("scal and/or off need to be of length 1")
  }
  list(scal=addArgs0$scal, off=addArgs0$off)
}

checkEntries <- function(modL, doses, fullMod){
  biModels <- c("emax", "linlog", "linear", "quadratic",
                "exponential", "logistic", "betaMod", "sigEmax",
                "linInt")
  checkNam <- function(nam){
    if(is.na(match(nam, biModels)))
      stop("Invalid model specified: ", nam)
  }
  checkStand <- function(nam){
    pars <- modL[[nam]]
    ## checks for as many invalid values as possible
    if(!is.numeric(pars) & !is.null(pars))
      stop("entries in Mods need to be of type: NULL, or numeric.\n",
           " invalid type specified for model ", nam)
    if((nam %in% c("linear", "linlog")) & !is.null(pars))
      stop("For model ", nam, ", model entry needs to be equal to NULL")
    if((nam %in% c("emax", "sigEmax", "betaMod", "logistic", "exponential")) & any(pars <= 0))
      stop("For model ", nam, " model entries needs to be positive")
    if((nam %in% c("emax", "exponential", "quadratic")) & is.matrix(nam))
      stop("For model ", nam, " parameters need to specified in a vector")
    if((nam %in% c("sigEmax", "betaMod", "logistic"))){
      if(is.matrix(pars)){
        if(ncol(pars) != 2)
          stop("Matrix for ", nam, " model needs to have two columns")
      }
      if(length(pars)%%2 > 0)
        stop("Specified parameters need to be a multiple of two for ", nam, " model")
    }
    if(nam == "linInt"){
      if(is.matrix(pars)){
        len <- ncol(pars)
      } else {
        len <- length(pars)
      }
      if(len != (length(doses)-1))
        stop("Need to provide guesstimates for each active dose. ", len,
             " specified, need ", length(doses)-1, ".")
    }
  }
  if(!fullMod){
    lapply(names(modL), function(nam){
      checkNam(nam)
      checkStand(nam)
    })
  } else {
    lapply(names(modL), function(nam){
      checkNam(nam)
    })
  }
}
  
## calculates parameters for all models in the candidate set returns a
## list with all model parameters.
fullMod <-  function(models, doses, placEff, maxEff, scal, off){
  ## check for valid placEff and maxEff arguments
  nM <- modCount(models, fullMod = FALSE)
  if(length(placEff) > 1){
    if(length(placEff) != nM)
      stop("placEff needs to be of length 1 or length equal to the number of models")
  } else {
    placEff <- rep(placEff, nM)
  }
  if(length(maxEff) > 1){
    if(length(maxEff) != nM)
      stop("maxEff needs to be of length 1 or length equal to the number of models")
  } else {
    maxEff <- rep(maxEff, nM)
  }
  nodes <- doses # nodes parameter for linInt
  
  ## calculate linear parameters of models (with standardized
  ## parameters as in models), to achieve the specified placEff and maxEff
  complMod <- vector("list", length=length(models))
  i <- 0;z <- 1
  for(nm in names(models)){
    pars <- models[[nm]]
    if(is.null(pars)){ ## linear and linlog
      Pars <- getLinPars(nm, doses, NULL, placEff[z], maxEff[z], off); i <- i+1; z <- z+1
    } 
    if(is.element(nm,c("emax", "exponential", "quadratic"))){
        nmod <- length(pars)
        if(nmod > 1){
          Pars <- matrix(ncol=3, nrow=nmod)
          for(j in 1:length(pars)){
            tmp <- getLinPars(nm, doses, as.vector(pars[j]), placEff[z], maxEff[z])
            Pars[j,] <- tmp
            z <- z+1
          }
          colnames(Pars) <- names(tmp)
          rownames(Pars) <- 1:length(pars)
          i <- i+1
        } else {
          Pars <-  getLinPars(nm, doses, as.vector(pars), placEff[z], maxEff[z])
          i <- i+1; z <- z+1
        }
      }
    if(is.element(nm,c("logistic", "betaMod", "sigEmax"))){
      if(is.matrix(pars)){
        Pars <- matrix(ncol=4, nrow=nrow(pars))
        for(j in 1:nrow(pars)){
          tmp <- getLinPars(nm, doses, as.vector(pars[j,]), placEff[z], maxEff[z])
          Pars[j,] <- tmp
          z <- z+1
        }
        colnames(Pars) <- names(tmp)
        rownames(Pars) <- 1:nrow(pars)
        i <- i+1
      } else {
        Pars <-  getLinPars(nm, doses, as.vector(pars), placEff[z], maxEff[z]); i <- i+1; z <- z+1
      }
    }
    if(nm == "linInt"){
      if(is.matrix(pars)){
        Pars <- matrix(ncol=length(nodes), nrow=nrow(pars))
        for(j in 1:nrow(pars)){
          Pars[j,] <-  getLinPars(nm, doses, as.vector(pars[j,]), placEff[z], maxEff[z])
          z <- z+1
        }
        colnames(Pars) <- paste("d", doses, sep="")
        rownames(Pars) <- 1:nrow(pars)
        i <- i+1
      } else {
        Pars <- getLinPars(nm, doses, as.vector(pars), placEff[z], maxEff[z]); i <- i+1; z <- z+1
        names(Pars) <- paste("d", doses, sep="")
      }
    }
    complMod[[i]] <- Pars
  }
  names(complMod) <- names(models)
  complMod
}

plotModels <- function(models, nPoints = 200, superpose = FALSE,
                       xlab = "Dose", ylab = "Model means",
                       modNams = NULL, plotTD = FALSE, Delta, ...){
  ## models is always assumed to be of class Mods 
  doses <- nodes <- attr(models, "doses")
  placEff <- attr(models, "placEff")
  maxEff <- attr(models, "maxEff")
  off <- attr(models, "off")
  scal <- attr(models, "scal")
  if(!inherits(models, "Mods"))
    stop("\"models\" needs to be of class Mods")
  nM <- modCount(models, fullMod = TRUE)
  
  if(nM > 50)
    stop("too many models in Mods object to plot (> 50 models).")
  
  doseSeq <- sort(union(seq(min(doses), max(doses), length = nPoints), 
                        doses))
  resp <- calcResp(models, doseSeq, off, scal, nodes)
  pdos <- NULL
  if(plotTD){ # also include TD in plot
    if(missing(Delta))
      stop("need Delta, if \"plotTD = TRUE\"")
    ind <- maxEff > 0
    if(length(unique(ind)) > 1)
      stop("inconsistent directions not possible, when \"plotTD = TRUE\"")
    direction <- ifelse(all(ind), "increasing", "decreasing")
    pdos <- TD(models, Delta, direction = direction)
    yax <- rep(ifelse(direction == "increasing", Delta, -Delta), length(pdos))
  }
  if(length(placEff) == 1)
    placEff <- rep(placEff, nM)
  if(length(maxEff) == 1)
    maxEff <- rep(maxEff, nM)
  if(is.null(modNams)){ # use alternative model names
    nams <- dimnames(resp)[[2]]
  } else {
    if(length(modNams) != nM)
      stop("specified model-names in \"modNams\" of invalid length")
    nams <- modNams
  }
  modelfact <- factor(rep(nams, each = length(doseSeq)),
                      levels = nams)
  if(superpose){
    respdata <- data.frame(response = c(resp),
                           dose = rep(doseSeq, ncol(resp)),
                           model = modelfact)
    spL <- lattice::trellis.par.get("superpose.line")
    spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
    spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
    spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
    ## data for plotting function within panel
    panDat <- list(placEff = placEff, maxEff = maxEff, doses = doses)
    ## number of columns
    nCol <- ifelse(nM < 5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))
    key <- list(lines = spL, transparent = TRUE,
                text = list(nams, cex = 0.9),
                columns = nCol)
    ltplot <- lattice::xyplot(response ~ dose, data = respdata, subscripts = TRUE, 
                              groups = respdata$model, panel.data = panDat, xlab = xlab,
                              ylab = ylab,
                              panel = function(x, y, subscripts, groups, ..., panel.data) {
                                lattice::panel.grid(h=-1, v=-1, col = "lightgrey", lty=2)
                                lattice::panel.abline(h = c(panel.data$placEff, panel.data$placEff + 
                                                              panel.data$maxEff), lty = 2)
                                lattice::panel.superpose(x, y, subscripts, groups, type = "l", ...)
                                ind <- !is.na(match(x, panel.data$doses))
                                lattice::panel.superpose(x[ind], y[ind], subscripts[ind], 
                                                         groups, ...)
                                if(plotTD){
                                  for(z in 1:length(pdos)){
                                    lattice::panel.lines(c(0, pdos[z]), c(yax[z], yax[z]),lty=2, col=2)
                                    lattice::panel.lines(c(pdos[z], pdos[z]), c(0, yax[z]),lty=2, col=2)
                                  }
                                }}, key = key, ...)
  } else {
    respdata <- data.frame(response = c(resp), 
                           dose = rep(doseSeq, ncol(resp)), model = modelfact)
    panDat <- list(placEff = placEff, maxEff = maxEff, doses = doses, pdos=pdos)
    ltplot <- lattice::xyplot(response ~ dose | model, data = respdata,
                              panel.data = panDat, xlab = xlab, ylab = ylab, 
                              panel = function(x, y, ..., panel.data){
                                lattice::panel.grid(h=-1, v=-1, col = "lightgrey", lty=2)
                                z <- panel.number()
                                lattice::panel.abline(h = c(panel.data$placEff[z],
                                                            panel.data$placEff[z] + 
                                                              panel.data$maxEff[z]), lty = 2)
                                lattice::panel.xyplot(x, y, type = "l", ...)
                                ind <- match(panel.data$doses, x)
                                lattice::panel.xyplot(x[ind], y[ind], ...)
                                if(plotTD){
                                  if(direction == "increasing"){
                                    delt <- Delta
                                    base <- panel.data$placEff[z]
                                    delt <- panel.data$placEff[z]+Delta
                                  } else {
                                    delt <- -Delta
                                    base <- panel.data$placEff[z]+panel.data$maxEff[z]
                                    delt <- panel.data$placEff[z]-Delta
                                  }
                                  lattice::panel.lines(c(0, pdos[z]), c(delt, delt), lty=2, col=2)
                                  lattice::panel.lines(c(pdos[z], pdos[z]), c(base, delt),lty=2, col=2)
                                }
                              }, strip = function(...) strip.default(..., style = 1), 
                              as.table = TRUE,...)
  }
  print(ltplot)
}


## calculate target dose
calcTD <- function(model, pars, Delta, TDtype = c("continuous", "discrete"),
                   direction = c("increasing", "decreasing"),
                   doses, off, scal, nodes){
  ## calculate the smallest dose x for which
  ## f(x) > f(0) + Delta (increasing) or f(x) < f(0) - Delta (decreasing)
  ## => f0(x) > Delta (increasing) or f0(x) < - Delta (decreasing) (f0 effect-curve)
  ## need to multiply f0(x) (=slope parameter) with -1 then decreasing case
  ## can be covered equivalent to increasing case
  TDtype <- match.arg(TDtype)
  direction <- match.arg(direction)
  if(direction == "decreasing"){ ## transform problem to "increasing" case
    if(model == "linInt"){
      pars <- -pars
    } else {
      pars[2] <- -pars[2]
      if(model == "quadratic") ## also need to negate pars[3]
        pars[3] <- -pars[3]
    }
  }
  if(model == "betaMod" & missing(scal))
    stop("Need \"scal\" parameter for betaMod model")
  if(model == "linlog" & missing(off))
    stop("Need \"off\" parameter for linlog model")    
  if(model == "linInt"){
    if(missing(nodes))
      stop("Need \"nodes\" parameter for linlog model")
    if(length(nodes) != length(pars))
      stop("nodes and pars of incompatible length")
  }
  
  if(TDtype == "continuous"){ ## calculate target dose analytically
    cf <- pars
    if(model == "linear"){
      td <- Delta/cf[2]
      if(td > 0)
        return(td)
      return(NA)
    }
    if(model == "linlog"){
      td <- off*exp(Delta/cf[2])-off
      if(td > 0)
        return(td)
      return(NA)
    }
    if(model == "quadratic"){
      if(4*cf[3]*Delta+cf[2]^2 < 0)
        return(NA)
      d1 <- -(sqrt(4*cf[3]*Delta+cf[2]^2)+cf[2])/(2*cf[3])
      d2 <- (sqrt(4*cf[3]*Delta+cf[2]^2)-cf[2])/(2*cf[3])       
      ind <- c(d1, d2) > 0
      if(!any(ind))
        return(NA)
      return(min(c(d1, d2)[ind]))
    }
    if(model == "emax"){
      if(Delta > cf[2])
        return(NA)
      return(Delta*cf[3]/(cf[2]-Delta))
    }
    if(model == "logistic"){
      if(Delta > cf[2] * (1 - logistic(0, 0, 1, cf[3], cf[4])))
        return(NA)
      .tmp1 <- exp(cf[3]/cf[4])
      num <- .tmp1*cf[2]-Delta*.tmp1-Delta
      den <- cf[2]+Delta*.tmp1+Delta
      return(cf[3]-cf[4]*log(num/den))
    }
    if(model == "sigEmax"){
      if(Delta > cf[2])
        return(NA)
      return((Delta*cf[3]^cf[4]/(cf[2]-Delta))^(1/cf[4]))
    }
    if(model == "betaMod"){
      if(Delta > cf[2])
        return(NA)
      func <- function(x, Emax, delta1, delta2, scal, Delta){
        betaMod(x, 0, 1, delta1, delta2, scal)-Delta/Emax
      }
      mode <- cf[3]/(cf[3]+cf[4])*scal
      out <- uniroot(func, lower=0, upper=mode, delta1=cf[3],
                     delta2=cf[4], Emax=cf[2], scal=scal,
                     Delta=Delta)$root
      return(out)
    }
    if(model == "exponential"){
      if(Delta/cf[2] < 0) ## wrong direction
        return(NA)
      return(cf[3]*log(Delta/cf[2]+1))
    }
    if(model == "linInt"){
      inds <- cf < cf[1] + Delta
      if(all(inds))
        return(NA)
      ind <- min((1:length(cf))[!inds])-1
      tmp <- (cf[1]+Delta-cf[ind])/(cf[ind+1]-cf[ind])
      td <- nodes[ind] + tmp*(nodes[ind+1]-nodes[ind])
      if(td > 0)
        return(td)
      else
        return(NA)
    }
  }
  if(TDtype == "discrete"){
    if(missing(doses))
      stop("For TDtype = \"discrete\" need the possible doses in doses argument")
    if(!any(doses == 0))
      stop("need placebo dose for TD calculation")
    if(model == "betaMod")
      pars <- c(pars, scal)
    if(model == "linlog")
      pars <- c(pars, off)
    doses <- sort(doses)
    if(model != "linInt"){
      resp <- do.call(model, c(list(doses), as.list(pars)))
    } else {
      resp <- do.call(model, c(list(doses), as.list(list(pars, nodes))))
    }
    ind <- resp >= resp[1] + Delta
    if(any(ind)){ ## TD does exist return smallest dose fulfilling threshold
      return(min(doses[ind]))
    } else {
      return(NA)
    }
  }
}

##  calculate gradient of target dose
calcTDgrad <- function(model, pars, Delta,
                       direction = c("increasing", "decreasing"), off, scal, nodes){
  direction <- match.arg(direction)
  if(direction == "decreasing"){ ## transform problem to "increasing" case
    Delta <- -Delta      ## TD is smallest x so that: 
  }                      ## f(x) = f(0) + Delta (incr), f(x) = f(0) - Delta (decr)
  cf <- pars
  if(model == "linear")
    return(c(0, -Delta/cf[2]^2))
  if(model == "linlog"){
    ## version assuming off unknown
    ##c(0, -Delta*off*exp(Delta/cf[2])/cf[2]^2, exp(Delta/cf[2])-1)
    return(c(0, -Delta*off*exp(Delta/cf[2])/cf[2]^2))
  }
  if(model == "quadratic"){
    squrt <- sqrt(4*Delta*cf[3]+cf[2]^2)
    .p1 <- -(squrt-cf[2])/(2*cf[3]*squrt)
    .p2 <- cf[2]*squrt-2*Delta*cf[3]-cf[2]^2
    .p2 <- .p2/(2*cf[3]^2*squrt)
    return(c(0, .p1, .p2))
  }
  if(model == "emax"){
    .p1 <- -Delta*cf[3]/(cf[2]-Delta)^2
    .p2 <- -Delta/((Delta/cf[2]-1)*cf[2])
    return(c(0, .p1, .p2))
  }
  if(model == "logistic"){
    et2t3 <- exp(cf[3]/cf[4])
    t1 <- (1/(1+et2t3)+Delta/cf[2])
    t2 <- (1/t1-1)
    .p1 <- -Delta*cf[4]/(cf[2]^2*t1^2*t2)
    .p2 <- 1-et2t3/((et2t3+1)^2*t1^2*t2)
    .p3 <- cf[3]*et2t3/(cf[4]*(et2t3+1)^2*t1^2*t2)-log(t2)
    return(c(0, .p1, .p2, .p3))
  }
  if(model == "sigEmax"){
    brack <- (-Delta*cf[3]^cf[4]/(Delta-cf[2]))^(1/cf[4])
    .p1 <- brack/((Delta-cf[2])*cf[4])
    .p2 <- brack/cf[3]
    .p3 <- brack*(log(cf[3])/cf[4]-log((-Delta*cf[3]^cf[4])/(Delta-cf[2]))/cf[4]^2)
    return(c(0, .p1, .p2, .p3))
  }
  if(model == "betaMod"){
    h0 <- function(cf, scal, Delta){
      func <- function(x, delta1, delta2, Emax, scal, Delta){
        betaMod(x, 0, 1, delta1, delta2, scal)-Delta/Emax
      }
      mode <- cf[3]/(cf[3]+cf[4])*scal
      uniroot(func, lower=0, upper=mode, delta1=cf[3], delta2=cf[4],
              Emax=cf[2], scal=scal, Delta=Delta)$root
    }
    td <- h0(cf, scal, Delta) ## calculate target dose
    .p1 <- -td*(scal-td)/(cf[2]*(cf[3]*(scal-td)-cf[4]*td))
    .p2 <- .p1*cf[2]*(log(td/scal)+log(cf[3]+cf[4])-log(cf[3]))
    .p3 <- .p1*cf[2]*(log(1-td/scal)+log(cf[3]+cf[4])-log(cf[4]))
    return(c(0, .p1, .p2, .p3))
  }
  if(model == "exponential"){
    .p1 <- -Delta*cf[3]/(cf[2]*Delta+cf[2]^2)
    .p2 <- log(Delta/cf[2] + 1)
    return(c(0, .p1, .p2))
  }
  if(model == "linInt"){
    stop("linInt model not implemented")
    ## ## the below should be correct
    ## out <- numeric(length(cf))
    ## indx <- 1:max(which(cf==max(cf)))
    ## ind <- max(indx[cf[indx] < cf[1] + Delta])
    ## out[1] <- 1/(cf[ind+1]-cf[ind])
    ## out[ind] <- -1/(cf[ind+1]-cf[ind])
    ## out[ind+1] <- -(cf[1]+Delta-cf[ind])/(cf[ind+1]-cf[ind])^2
    ## return(out*(nodes[ind+1]-nodes[ind]))
  }
}

calcED <- function(model, pars, p, maxD, EDtype = c("continuous", "discrete"),
                   doses, off, scal, nodes){
  ## calculate the smallest dose x for which
  ## f(x) > f(0) + p*(f(xmax)-f(0))
  ## e.g. the EDp within the observed dose-range
  EDtype <- match.arg(EDtype)
  if(model == "betaMod" & missing(scal))
    stop("Need \"scal\" parameter for betaMod model")
  if(model == "linlog" & missing(off))
    stop("Need \"off\" parameter for linlog model")    
  if(model == "linInt"){
    if(missing(nodes))
      stop("Need \"nodes\" parameter for linlog model")
    if(length(nodes) != length(pars))
      stop("nodes and pars of incompatible length")
  }
  
  if(EDtype == "continuous"){ ## calculate target dose analytically
    cf <- pars
    if(cf[2] == 0){
      return(NA)
    }
    if(model == "linear"){
      return(p*maxD)
    }
    if(model == "linlog"){
      return(off*(exp(p*(log(maxD+off)-log(off)))-1))
    }
    if(model == "exponential"){
      return(cf[3]*log(p*exp(maxD/cf[3])-p+1))
    }
    if(model == "emax"){
      return(p*cf[3]*maxD/((1-p)*maxD+cf[3]))
    }
    if(model == "logistic"){
      res1 <- ((p-1)*exp(maxD/cf[4]+cf[3]/cf[4])-exp(2*cf[3]/cf[4])-p*exp(cf[3]/cf[4]))
      res2 <- ((p*exp(cf[3]/cf[4])+1)*exp(maxD/cf[4])+(1-p)*exp(cf[3]/cf[4]))
      return(cf[3]-cf[4]*log(-res1/res2))
    }
    if(model == "sigEmax"){
      out <-  p*cf[3]^cf[4]*maxD^cf[4]/((1-p)*maxD^cf[4]+cf[3]^cf[4])
      return(out^(1/cf[4]))
    }
    if(model == "quadratic"){
      mode <- -pars[2]/(2*pars[3])
      if(mode > maxD | mode < 0) ## maximum outside dose range
        mode <- maxD
      const <- pars[2]*mode+pars[3]*mode^2
      d1 <- -(sqrt(4*pars[3]*const*p+pars[2]^2)+pars[2])/pars[3]/2.0
      d2 <- (sqrt(4*pars[3]*const*p+pars[2]^2)-pars[2])/pars[3]/2.0
      ind <- c(d1, d2) > 0
      if(!any(ind))
        return(NA)
      return(min(c(d1, d2)[ind]))
    }
    if(model == "betaMod"){
      func <- function(x, Emax, delta1, delta2, scal, p, mode){
        p - betaMod(x, 0, 1, delta1, delta2, scal)/betaMod(mode, 0, 1, delta1, delta2, scal)
      }
      mode <- cf[3]/(cf[3]+cf[4])*scal
      out <- uniroot(func, lower=0, upper=mode, delta1=cf[3],
                     delta2=cf[4], Emax=cf[2], scal=scal,
                     p=p, mode = mode)$root
      return(out)
    }
    if(model == "linInt"){
      dif <- cf-cf[1]
      ind <- which.max(abs(dif))
      maxEff <- abs(dif)[ind]
      if(dif[ind] > 0){
        direc <- "increasing"
      } else {
        direc <- "decreasing"
      }
      out <- calcTD("linInt", cf, Delta=p*maxEff, TDtype="continuous",
                    direction = direc, off=off, scal=scal, nodes=nodes)
      return(out)
    }
  }
  if(EDtype == "discrete"){
    ## use calcTD function
    if(missing(doses))
      stop("For EDtype = \"discrete\" need the possible doses in doses argument")
    if(!any(doses == 0))
      stop("need placebo dose for ED calculation")
    doses <- sort(doses)
    if(model != "linInt"){
      if(model == "betaMod")
        pars <- c(pars, scal)
      if(model == "linlog")
        pars <- c(pars, off)
      resp0 <- do.call(model, c(list(0), as.list(pars)))
      resp <- abs(do.call(model, c(list(doses), as.list(pars)))-resp0)
    } else {
      resp0 <- do.call(model, c(list(0), as.list(list(pars, nodes))))
      resp <- abs(do.call(model, c(list(doses), as.list(list(pars, nodes))))-resp0)
    }
    ## calculate maximum response
    if(model %in% c("betaMod", "quadratic")){
      func2 <- function(x){
        resp0 <- do.call(model, c(list(0), as.list(pars)))
        abs(do.call(model, c(list(x), as.list(pars)))-resp0)
      }
      opt <- optimize(func2, range(doses), maximum=TRUE)
      maxResp <- opt$objective
    } else {
      maxResp <- max(resp)
    }
  }
  ind <- resp >= p*maxResp
  if(any(ind)){ ## TD does exist return smallest dose fulfilling threshold
    return(min(doses[ind]))
  } else {
    return(NA)
  }
}


calcEDgrad <- function(model, pars, maxD, p, off, scal, nodes){
  cf <- pars
  if(model == "linear")
    return(c(0,0))
  if(model == "linlog"){
    return(c(0,0))
  }
  if(model == "emax"){
    p <- (1-p)*p*maxD^2/(p*maxD-maxD-cf[3])^2
    return(c(0, 0, p))
  }
  if(model == "exponential"){
    p <- log(p*exp(maxD/cf[3])-p+1)-p*maxD*exp(maxD/cf[3])/(cf[3]*(p*exp(maxD/cf[3])-p+1))
    return(c(0, 0, p))
  }
  ## for other models calculate gradient numerically (formulas more complicated)
  if(model == "linInt"){
    stop("linInt model not implemented")
  }
  avail <- requireNamespace("numDeriv", quietly = TRUE)
  if(!avail)
    stop("Need suggested package numDeriv for this calculation")
  func0 <- function(pars, model, p, maxD, off, scal){
    calcED(model, pars, p, maxD, EDtype = "continuous", off=off, scal=scal)
  }
  scal0 <- off0 <- NULL
  if(model == "betaMod")
    scal0 <- scal
  if(model == "linlog")
    off0 <- off
  numDeriv::grad(func0, pars, model=model, p=p, maxD=maxD, off=off, scal=scal)
}


calcResp <- function(models, doses, off, scal, nodes){
  ## generate response vectors for models and guesstimates in "models"
  ## models - candidate model list of class Mods
  nModels <- length(models)             # number of model elements
  parList <- val <- vector("list", modCount(models, fullMod = TRUE))
  k <- 1
  nams <- character()
  for(nm in names(models)) {
    pars <- models[[nm]]
    if (!is.null(pars) && !is.numeric(pars)) {
      stop("elements of \"models\" must be NULL or numeric")
    }
    if (is.matrix(pars)) {            # multiple models
      nmod <- nrow(pars)              # number of models
      if(nm == "linlog")
        pars <- cbind(pars, off)
      if(nm == "betaMod")
        pars <- cbind(pars, scal)
      ind <- 1:nmod
      nams <- c(nams, paste(nm, ind, sep = ""))
      for(j in 1:nmod) {
        if(nm != "linInt"){
          val[[k]] <- do.call(nm, c(list(doses), as.list(pars[j,])))
        } else {
          val[[k]] <- linInt(doses, pars[j,], nodes)
        }
        parList[[k]] <- pars[j,]
        k <- k + 1
      }
    } else {                      # single model
      if(nm == "linlog")
        pars <- c(pars, off)
      if(nm == "betaMod")
        pars <- c(pars, scal)
      nams <- c(nams, nm)
      if(nm != "linInt"){
        val[[k]] <- do.call(nm, c(list(doses), as.list(pars)))
      } else {
        val[[k]] <- linInt(doses, pars, nodes)
      }
      parList[[k]] <- pars
      k <- k + 1
    }       
  }
  muMat <- do.call("cbind", val)
  dimnames(muMat) <- list(doses, nams)
  names(parList) <- nams
  attr(muMat, "parList") <- parList
  muMat
}

## calculates the location and scale parameters corresponding to
## given placEff, maxEff, and guesstimates
getLinPars <- function(model, doses, guesstim, placEff, maxEff, off, scal){
  if(model == "linear"){
    e1 <- maxEff/max(doses)
    return(c(e0=placEff, delta=e1))
  }
  if(model == "linlog"){
    e1 <- maxEff/(log(max(doses) + off) - log(off))
    return(c(e0=(placEff-e1*log(off)), delta=e1))
  }
  if(model == "quadratic"){
    dMax <- 1/(-2*guesstim)
    b1 <- maxEff/(dMax + guesstim*dMax^2)
    b2 <- guesstim * b1
    return(c(e0=placEff, b1=b1, b2=b2))
  }
  if(model == "emax"){
    emax.p <- maxEff * (guesstim + max(doses))/max(doses)
    return(c(e0=placEff, eMax=emax.p, ed50=guesstim))
  }
  if(model == "exponential"){
    e1 <- maxEff/(exp(max(doses)/guesstim) - 1)
    e0 <- placEff
    return(c(e0=e0, e1=e1, delta=guesstim))
  }
  if(model == "logistic"){
    emax.p <- maxEff/
      (logistic(max(doses),0,1, guesstim[1], guesstim[2]) -
       logistic(0, 0, 1, guesstim[1], guesstim[2]))
    e0 <- placEff-emax.p*logistic(0,0,1,guesstim[1], guesstim[2])
    return(c(e0=e0, eMax=emax.p, ed50=guesstim[1], delta=guesstim[2]))
  }
  if(model == "betaMod"){
    return(c(e0=placEff, eMax=maxEff, delta1=guesstim[1], delta2=guesstim[2]))
  }
  if(model == "sigEmax"){
    ed50 <- guesstim[1]
    h <- guesstim[2]
    dmax <- max(doses)
    eMax <- maxEff*(ed50^h+dmax^h)/dmax^h
    return(c(e0 = placEff, eMax = eMax, ed50 = ed50, h = h))
  }
  if(model == "linInt"){
    ind <- which.max(abs(guesstim))
    return(c(placEff, placEff+maxEff*guesstim/guesstim[ind]))
  }
}

getModNams <- function(parList){
  ## extract model names with parameter values
  nM <- length(parList)
  mod_nams <- names(parList)
  for(i in 1:nM){
    if(startsWith(mod_nams[i], "linlog"))
      mod_nams[i] <- sprintf("linlog (off=%s)", parList[[i]][3])
    if(startsWith(mod_nams[i], "emax"))
      mod_nams[i] <- sprintf("emax (ED50=%s)", parList[[i]][3])
    if(startsWith(mod_nams[i], "exponential"))
      mod_nams[i] <- sprintf("exponential (delta=%s)", parList[[i]][3])
    if(startsWith(mod_nams[i], "quadratic"))
      mod_nams[i] <- sprintf("quadratic (delta=%s)", parList[[i]][3]/parList[[i]][2])
    if(startsWith(mod_nams[i], "sigEmax"))
      mod_nams[i] <- sprintf("sigEmax (ED50=%s,h=%s)", parList[[i]][3], parList[[i]][4])
    if(startsWith(mod_nams[i], "logistic"))
      mod_nams[i] <- sprintf("logistic (ED50=%s,delta=%s)",
                             parList[[i]][3], parList[[i]][4])
    if(startsWith(mod_nams[i], "betaMod"))
      mod_nams[i] <- sprintf("betaMod (delta1=%s,delta2=%s,scal=%s)",
                             parList[[i]][3], parList[[i]][4], parList[[i]][5])
      if(startsWith(mod_nams[i], "linInt"))
        mod_nams[i] <- sprintf("linInt (%s)", paste0(parList[[i]], collapse=","))
  }
  mod_nams
}

