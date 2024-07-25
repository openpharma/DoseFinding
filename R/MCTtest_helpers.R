checkAnalyArgs <- function(dose, resp, data, S, type,
                           addCovars, placAdj, na.action, cal){
  if(!inherits(addCovars, "formula"))
    stop("addCovars argument needs to be of class formula")
  if(!inherits(placAdj, "logical"))
    stop("placAdj argument needs to be of class logical")
  if(placAdj){
    if(type == "normal")
      stop("\"placAdj == TRUE\" only allowed for type = \"general\"")
  }
  if(!is.null(data)){ # data handed over in data frame
    if(!is.data.frame(data))
      stop("data argument needs to be a data frame")
    nams <- c(cal[2], cal[3], all.vars(addCovars))
    ind <- match(nams, names(data))
    if (any(is.na(ind)))
      stop("variable(s): ", paste(nams[is.na(ind)], collapse= ", "), " not found in ", cal[4])
    dd <- na.action(data[,nams])
  } else { # data handed over via vectors
    if(addCovars != ~1)
      stop("need to hand over data and covariates in data frame")
    if(!(is.numeric(resp) && is.null(dim(resp)))) {
      warning(cal[3], " is not a numeric but a ", class(resp)[1],
              ", converting with as.numeric()")
      resp <- as.numeric(resp)
    }
    if(length(dose) != length(resp))
      stop(cal[2], " and ", cal[3], " not of equal length")
    dd <- na.action(data.frame(dose, resp))
    cal[2:3] <- gsub("\\$", "", cal[2:3])
    cal[2:3] <- gsub("\\[|\\]", "", cal[2:3])
    colnames(dd) <- cal[2:3]
  }
  doseNam <- cal[2];respNam <- cal[3]
  if(placAdj){
    if(any(dd[[doseNam]] == 0))
      stop("If placAdj == TRUE there should be no placebo group")
  }
  if(any(dd[[doseNam]] < -.Machine$double.eps))
    stop("dose values need to be non-negative")
  if(!is.numeric(dd[[doseNam]]))
    stop("dose variable needs to be numeric")
  if(!is.numeric(dd[[respNam]]))
    stop("response variable needs to be numeric")
  ## check type related arguments
  if(type == "general" & is.null(S))
    stop("S argument missing")
  if(type == "normal" & !is.null(S))
    message("Message: S argument ignored for type == \"normal\"\n")
  if(type == "general" & addCovars != ~1)
    message("Message: addCovars argument ignored for type == \"general\"")
  if(!is.null(S)){
    if(!is.matrix(S))
      stop("S needs to be of class matrix")
    nD <- length(dd[[doseNam]])
    if(nrow(S) != nD | ncol(S) != nD)
      stop("S and dose have non-conforming size")
  }
  ord <- order(dd[[doseNam]])
  dd <- dd[ord, ]
  Sout <- NULL
  if(type == "general")
    Sout <- S[ord, ord]
  return(list(dd=dd, type = type, S = Sout, ord=ord,
              doseNam=doseNam, respNam=respNam))
}
