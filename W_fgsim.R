# Wrapper functions to be used with hhh4:
# Takes settings and parameters and returns functions that compute weights
# and their derivatives

# Misc functions
removeNAmat = function(x){
  x[is.na(x)] = 0
  x[is.infinite((x))] = 0
  return(x)
}

removeNAlist = function(x){
  x = lapply(x, removeNAmat)
  return(x)
}

# Input:
# pars:      List of parameters with elements mu, sigma, w, l, u
# truncPL:   Truncated power law Y/N
# maxlag:    Maximum spatial lag considered
# normalize: Row normalisation Y/N
# log:       Log-transformation of parameters Y/N
# initial:   Initial parameters
# from0:     Consider intra-district weights Y/N
# popScale:  Scale by population Y/N (if yes, pars should contain "pop")
#            Should have determinant 1 nor numerical stability

W_fgsim = function(pars,
                   truncPL = FALSE,
                   gravity = FALSE,
                   maxlag = Inf,
                   normalize = TRUE,
                   log = FALSE,
                   initial = rep(ifelse(log, 0, 1), 1 + (truncPL | gravity)),
                   from0 = TRUE,
                   popScale = FALSE){
  
  if(gravity & truncPL){
    warning("Truncated power law not supported for gravity model, set to FALSE.")
    truncPL = FALSE
  }
  if(gravity){
    target = c("mu1", "mu2", "sigma1", "sigma12", "sigma2", "w")
  }else{
    target = c("mu", "sigma", "w", "l", "u")
  }
  
  pn = names(pars)
  for(i in 1:length(target)){
    if(!(target[i] %in% pn)) stop(paste("List pars must contain an element", target[i]))
  }
  if(popScale){
    if(!("pop" %in% pn)) stop(paste("List pars must contain a matrix pop"))
  }
  
  if(truncPL & log){
    warning("Truncation threshold for power law is already defined on log distances.")
  }
  
  # Initialise by computing raw weights
  if(truncPL){
    weights.call = call("WrawMTLtrunc", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 0)
    weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                          envir=getNamespace("surveillance"))
    
    dweights.call = call("WrawMTLtrunc", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 1)
    dweights = as.function(c(alist(d=, nbmat=, ...=), call("{", dweights.call)),
                           envir=getNamespace("surveillance"))
    
    d2weights.call = call("WrawMTLtrunc", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 2)
    d2weights = as.function(c(alist(d=, nbmat=, ...=), call("{", d2weights.call)),
                            envir=getNamespace("surveillance"))
  }else{
    if(gravity){
      weights.call = call("WrawMN", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 0)
      weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                            envir=getNamespace("surveillance"))
      
      dweights.call = call("WrawMN", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 1)
      dweights = as.function(c(alist(d=, nbmat=, ...=), call("{", dweights.call)),
                             envir=getNamespace("surveillance"))
      
      d2weights.call = call("WrawMN", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 2)
      d2weights = as.function(c(alist(d=, nbmat=, ...=), call("{", d2weights.call)),
                              envir=getNamespace("surveillance"))
    }else{
      weights.call = call("WrawMTL", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 0)
      weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                            envir=getNamespace("surveillance"))
      
      dweights.call = call("WrawMTL", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 1)
      dweights = as.function(c(alist(d=, nbmat=, ...=), call("{", dweights.call)),
                             envir=getNamespace("surveillance"))
      
      d2weights.call = call("WrawMTL", quote(nbmat), quote(pars), quote(d), ifelse(from0, 0, 1), maxlag, 2)
      d2weights = as.function(c(alist(d=, nbmat=, ...=), call("{", d2weights.call)),
                              envir=getNamespace("surveillance"))
    }
  }
  
  
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
  
  # Read in pars from global environment
  # Slower than forwarding it from hhh4(), but feasible
  nm = deparse(substitute(pars))
  body(weights) <- as.call(append(as.list(body(weights)),
                                  substitute(pars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                   substitute(pars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    substitute(pars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  
  # Transform decay parameter if it is in logs
  if (log){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(d <- exp(d)), after = 1))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(d <- exp(d)), after = 1))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(d <- exp(d)), after = 1))
  }
  
  # Scale by population
  if(popScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperSC(pars$pop, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperSC(pars$pop, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperSC(pars$pop, out$W, out$W1, out$W2))))
  }
  
  # Row normalisation
  if(normalize){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperRN(out$W, out$W1, out$W2))))
  }
  
  # Transform weights/derivatives if log-transformation of parameters
  if(log){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperLOG(d, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperLOG(d, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperLOG(d, out$W, out$W1, out$W2))))
  }
  
  # Remove NA and infinite values (just to be sure)
  body(weights) <- as.call(append(as.list(body(weights)),
                                  quote(removeNAmat(out$W))))
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                   quote(removeNAlist(out$W1))))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(removeNAlist(out$W2))))
  
  return(list(w = weights, dw = dweights, d2w = d2weights, initial = initial))
}
