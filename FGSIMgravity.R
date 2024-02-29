###################################################
### Methods to compute the weights of the       ###
### Fine-Grid Spatial Interaction Model (FGSIM) ###
### For gravity model weights                   ###
###################################################

W_gravity = function(pars,
                   maxlag = Inf,
                   normalize = TRUE,
                   log = FALSE,
                   initial = if(log) c(0, 0, 0) else c(1, 1, 1),
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
  
  target = c("dist", "pO", "pD", "s11", "s12", "s13", "s22", "s23", "s33")
  pn = names(pars)
  for(i in 1:9){
    if(!(target[i] %in% pn)) stop(paste("List pars must contain an element", target[i]))
  }
  if(areaScale & (!("area" %in% pn))) stop("pars must contain an element area for areaScale")
  if(popScale & (!("pop" %in% pn))) stop("pars must contain an element pop for popScale")
  if(contactScale & (!("M" %in% pn))) stop("pars must contain an element M for contactScale")
  
  # Initialise
  weights.call = call("WrawGravity", quote(nbmat), quote(gravityPars), quote(d), ifelse(from0, 0, 1), maxlag, 0)
  weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                        envir=getNamespace("surveillance"))
  
  dweights.call = call("WrawGravity", quote(nbmat), quote(gravityPars), quote(d), ifelse(from0, 0, 1), maxlag, 1)
  dweights = as.function(c(alist(d=, nbmat=, ...=), call("{", dweights.call)),
                         envir=getNamespace("surveillance"))
  
  d2weights.call = call("WrawGravity", quote(nbmat), quote(gravityPars), quote(d), ifelse(from0, 0, 1), maxlag, 2)
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
                                  substitute(gravityPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                   substitute(gravityPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    substitute(gravityPars <- get(nameHere, globalenv()), list(nameHere = nm)), after = 1))
  
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
                                    quote(out <- helperSC(gravityPars$area, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperSC(gravityPars$area, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperSC(gravityPars$area, out$W, out$W1, out$W2))))
  }
  
  if(contactScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperRN(out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperRN(out$W, out$W1, out$W2))))
    
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperC(gravityPars$M, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperC(gravityPars$M, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperC(gravityPars$M, out$W, out$W1, out$W2))))
  }
  
  if(popScale){
    body(weights) <- as.call(append(as.list(body(weights)),
                                    quote(out <- helperSC(gravityPars$pop, out$W, out$W1, out$W2))))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(out <- helperSC(gravityPars$pop, out$W, out$W1, out$W2))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(out <- helperSC(gravityPars$pop, out$W, out$W1, out$W2))))
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
                                   quote(out$W1)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(out$W2)))
  
  return(list(w = weights, dw = dweights, d2w = d2weights, initial = initial))
}
