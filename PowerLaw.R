###################################################
### Methods to compute the weights of the       ###
### Fine-Grid Spatial Interaction Model (FGSIM) ###
### For gravity model weights                   ###
###################################################

W_powerlaw2 = function(maxlag,
                       powerlawPars,
                       normalize = TRUE,
                       log = FALSE,
                       initial = if(log) 0 else 1,
                       from0 = TRUE,
                       CT = FALSE,
                       PT = FALSE){
  
  # Initialise
  if(from0) maxlag = maxlag + 1
  
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
  
  nm = deparse(substitute(powerlawPars))
  body(weights) <- as.call(append(as.list(body(weights)), 
                                  substitute(powerlawPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(dweights) <- as.call(append(as.list(body(dweights)), 
                                   substitute(powerlawPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                    substitute(powerlawPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))

  # Transform decay parameter if it is in logs
  if(log){
    body(weights) = as.call(append(as.list(body(weights)), 
                                   quote(d <- exp(d)), after = 1))
    body(dweights) = as.call(append(as.list(body(dweights)), 
                                     quote(d <- exp(d)), after = 1))
    body(d2weights) = as.call(append(as.list(body(d2weights)), 
                                      quote(d <- exp(d)), after = 1))
  }
  
  # Transform neighbouring order if from0
  if(from0){
    body(weights) = as.call(append(as.list(body(weights)), 
                                   quote(nbmat <- nbmat + 1), after = 1))
    body(dweights) = as.call(append(as.list(body(dweights)), 
                                   quote(nbmat <- nbmat + 1), after = 1))
    body(d2weights) = as.call(append(as.list(body(d2weights)), 
                                   quote(nbmat <- nbmat + 1), after = 1))
  }
  
  if(PT){
    body(weights) <- as.call(append(as.list(body(weights)), 
                                    quote(out <- helperPOPdist(powerlawPars$pop, out, depth = 0))))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(out <- helperPOPdist(powerlawPars$pop, out$W, out$W1, depth = 1))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(out <- helperPOPdist(powerlawPars$pop, out$W, out$W1, out$W2, depth = 2))))
  }
  
  # Normalize weights if selected
  if(normalize | CT){
    body(weights) <- as.call(append(as.list(body(weights)), 
                                    quote(out <- helperRNdist(out, depth = 0))))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(out <- helperRNdist(out$W, out$W1, depth = 1))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(out <- helperRNdist(out$W, out$W1, out$W2, depth = 2))))
  }
  
  if(CT){
    body(weights) <- as.call(append(as.list(body(weights)), 
                                    quote(out <- helperCdist(powerlawPars$M, out, depth = 0))))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(out <- helperCdist(powerlawPars$M, out$W, out$W1, depth = 1))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(out <- helperCdist(powerlawPars$M, out$W, out$W1, out$W2, depth = 2))))
    
    if(normalize){
      body(weights) <- as.call(append(as.list(body(weights)), 
                                      quote(out <- helperRNdist(out, depth = 0))))
      body(dweights) <- as.call(append(as.list(body(dweights)), 
                                       quote(out <- helperRNdist(out$W, out$W1, depth = 1))))
      body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                        quote(out <- helperRNdist(out$W, out$W1, out$W2, depth = 2))))
    }
  }
  
  if(log){
    body(weights) <- as.call(append(as.list(body(weights)), 
                                    quote(out <- helperLOGdist(d, out, depth = 0))))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(out <- helperLOGdist(d, out$W, out$W1, depth = 1))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(out <- helperLOGdist(d, out$W, out$W1, out$W2, depth = 2))))
  }
  
  body(weights) <- as.call(append(as.list(body(weights)), 
                                  quote(out)))
  body(dweights) <- as.call(append(as.list(body(dweights)), 
                                   quote(out$W1)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                    quote(out$W2)))
  
  return(list(w = weights, dw = dweights, d2w = d2weights, initial = initial))
}
