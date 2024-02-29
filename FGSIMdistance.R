######################################################
### Fine-Grid Spatial Interaction Matrix           ###
### Based on distances only                        ###
### Wrapper function that returns three functions: ###
### Weights, and its first two derivatives         ###
### and initial value for the parameter            ###
######################################################

# Settings:
# maxlag:     To truncate the maximum lag considered
# pars:       A list that is in the global environment and contains
#             dist & sigma2 & M (if CT) & pop (if PT)
# normalize:  Shall weights be row-normalised?
# log:        If true, parameters are sign restricted
# initial:    Initial values for parameters
# from0:      If FALSE, intra-district weights are zero
# CT:         Contact transformation?
# PT:         Population transformation?

W_pdist = function(pars,
                   maxlag = Inf,
                   normalize = TRUE,
                   log = FALSE,
                   initial = ifelse(log, 0, 1),
                   from0 = TRUE,
                   areaScale = FALSE,
                   contactScale = FALSE,
                   popScale = FALSE){
  if(contactScale){
    if(!areaScale){
      areaScale = T
      warning("areaScale set to TRUE for contact scaling")
    }
    if(!popScale){
      popScale = T
      warning("popScale set to TRUE for contact scaling")
    }
  }
  
  pn = names(pars)
  if(!("dist" %in% pn)) stop("List pars must contain an element dist")
  if(!("s11" %in% pn)) stop("List pars must contain an element s11")
  if(areaScale & (!("area" %in% pn))) stop("pars must contain an element area for areaScale")
  if(popScale & (!("pop" %in% pn))) stop("pars must contain an element pop for popScale")
  if(contactScale & (!("M" %in% pn))) stop("pars must contain an element M for contactScale")
  
  # Initialise
  weights.call = call("WrawDistance", quote(nbmat), quote(distPars), quote(d), ifelse(from0, 0, 1), maxlag, 0)
  weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                        envir=getNamespace("surveillance"))
  
  dweights.call = call("WrawDistance", quote(nbmat), quote(distPars), quote(d), ifelse(from0, 0, 1), maxlag, 1)
  dweights = as.function(c(alist(d=, nbmat=, ...=), call("{", dweights.call)),
                         envir=getNamespace("surveillance"))
  
  d2weights.call = call("WrawDistance", quote(nbmat), quote(distPars), quote(d), ifelse(from0, 0, 1), maxlag, 2)
  d2weights = as.function(c(alist(d=, nbmat=, ...=), call("{", d2weights.call)),
                          envir=getNamespace("surveillance"))
  
  # Assign weights to an object Wraw
  body(weights) = as.call(c(as.name("{"),
                            substitute(out <- weights.call, list(weights.call=weights.call))
  ))
  
  body(dweights) = as.call(c(as.name("{"),
                             substitute(out <- dweights.call, list(dweights.call=dweights.call))
  ))
  
  body(d2weights) = as.call(c(as.name("{"),
                              substitute(out <- d2weights.call, list(d2weights.call=d2weights.call))
  ))
  
  # Read in distPars from global environment
  # Slower than forwarding it from hhh4(), but feasible
  nm = deparse(substitute(pars))
  body(weights) <- as.call(append(as.list(body(weights)),
                                  substitute(distPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                  substitute(distPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                  substitute(distPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  
  # Transform decay parameter if it is in logs
  if (log){
    body(weights) <- as.call(append(as.list(body(weights)),
                                   quote(d <- exp(d)), after = 1))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(d <- exp(d)), after = 1))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(d <- exp(d)), after = 1))
  }
  
  if(areaScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperSC(distPars$area, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperSC(distPars$area, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperSC(distPars$area, out$W, out$W1, out$W2))))
  }
  
  if(contactScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperRN(out$W, out$W1, out$W2))))
    
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperC(distPars$M, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperC(distPars$M, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperC(distPars$M, out$W, out$W1, out$W2))))
  }
  
  if(popScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperSC(distPars$pop, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperSC(distPars$pop, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperSC(distPars$pop, out$W, out$W1, out$W2))))
  }
  
  if(normalize){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperRN(out$W, out$W1, out$W2))))
  }
  
  if(log){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperLOG(d, out$w, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperLOG(d, out$w, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperLOG(d, out$w, out$W1, out$W2))))
  }
  
  body(weights) <- as.call(append(as.list(body(weights)),
                                  quote(out$W)))
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                   quote(out$W1[[1]])))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(out$W2[[1]])))
  
  return(list(w = weights, dw = dweights, d2w = d2weights, initial = initial))
}
