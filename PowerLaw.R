###################################################
### Methods to compute the weights of the       ###
### Fine-Grid Spatial Interaction Model (FGSIM) ###
### For gravity model weights                   ###
###################################################

# Settings:
# pars:          A list that is in the global environment and may contain
#                area, M and/or pop
# maxlag:        To truncate the maximum lag considered
# normalize:     Shall weights be row-normalised?
# log:           If true, parameters are sign restricted
# initial:       Initial values for parameters
# from0:         If FALSE, intra-district weights are zero
# areaScale:     Apply area scaling?
# contactScale:  Apply area scaling?
# popScale:      Apply area scaling?

W_pdist = function(pars = NA,
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
  if(any(is.na(pars)) & any(c(areaScale, popScale, contactScale))) stop("pars must be provided")
  if(areaScale & (!("area" %in% pn))) stop("pars must contain an element area for areaScale")
  if(popScale & (!("pop" %in% pn))) stop("pars must contain an element pop for popScale")
  if(contactScale & (!("M" %in% pn))) stop("pars must contain an element M for contactScale")
  
  # Initialise
  weights.call = call("WrawPowerlaw", quote(nbmat), quote(d), ifelse(from0, 0, 1), maxlag, 0)
  weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                        envir=getNamespace("surveillance"))
  
  dweights.call = call("WrawPowerlaw", quote(nbmat), quote(d), ifelse(from0, 0, 1), maxlag, 1)
  dweights = as.function(c(alist(d=, nbmat=, ...=), call("{", dweights.call)),
                         envir=getNamespace("surveillance"))
  
  d2weights.call = call("WrawPowerlaw", quote(nbmat), quote(d), ifelse(from0, 0, 1), maxlag, 2)
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
                                  substitute(powerlawPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                   substitute(powerlawPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    substitute(powerlawPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  
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
                                    quote(out <- helperSC(powerlawPars$area, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperSC(powerlawPars$area, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperSC(powerlawPars$area, out$W, out$W1, out$W2))))
  }
  
  if(contactScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperRN(out$W, out$W1, out$W2))))
    
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperC(powerlawPars$M, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperC(powerlawPars$M, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperC(powerlawPars$M, out$W, out$W1, out$W2))))
  }
  
  if(popScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperSC(powerlawPars$pop, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperSC(powerlawPars$pop, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperSC(powerlawPars$pop, out$W, out$W1, out$W2))))
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
