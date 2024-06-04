

checkAnalyArgs_bMCP <- function (dose, resp, data, S, type, prior, na.action, cal) 
{
  
  if (!is.null(data)) {
    if (!is.data.frame(data)) 
      stop("data argument needs to be a data frame")
    nams <- c(cal[2], cal[3])
    ind <- match(nams, names(data))
    if (any(is.na(ind))) 
      stop("variable(s): ", paste(nams[is.na(ind)], collapse = ", "), 
           " not found in ", cal[4])
    dd <- na.action(data[, nams])
  }
  else {
    if (!(is.numeric(resp) && is.null(dim(resp)))) {
      warning(cal[3], " is not a numeric but a ", class(resp)[1], 
              ", converting with as.numeric()")
      resp <- as.numeric(resp)
    }
    if (length(dose) != length(resp)) 
      stop(cal[2], " and ", cal[3], " not of equal length")
    dd <- na.action(data.frame(dose, resp))
    cal[2:3] <- gsub("\\$", "", cal[2:3])
    cal[2:3] <- gsub("\\[|\\]", "", cal[2:3])
    colnames(dd) <- cal[2:3]
  }
  doseNam <- cal[2]
  respNam <- cal[3]
  if (any(dd[[doseNam]] < -.Machine$double.eps)) 
    stop("dose values need to be non-negative")
  if (!is.numeric(dd[[doseNam]])) 
    stop("dose variable needs to be numeric")
  if (!is.numeric(dd[[respNam]])) 
    stop("response variable needs to be numeric")
  if (type == "general" & is.null(S)) 
    stop("S argument missing")
  if (type == "normal" & !is.null(S)) 
    message("Message: S argument ignored for type == \"normal\"\n")
  if (!is.null(S)) {
    if (!is.matrix(S)) 
      stop("S needs to be of class matrix")
    nD <- length(dd[[doseNam]])
    if (nrow(S) != nD | ncol(S) != nD) 
      stop("S and dose have non-conforming size")
  }
  ord <- order(dd[[doseNam]])
  dd <- dd[ord, ]
  Sout <- NULL
  if (type == "general") 
    Sout <- S[ord, ord]

  if (length(unique(dd[[doseNam]])) != length(prior)) 
    stop("Dose and prior have non-conforming size")
  if (!all(unlist(lapply(prior, function(x) "normMix" %in% class(x))))) 
    stop("priors need to be of class normMix")
  
  return(list(dd = dd, type = type, S = Sout, ord = ord, doseNam = doseNam, 
              respNam = respNam))
}
