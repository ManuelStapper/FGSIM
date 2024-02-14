###################################################
### Methods to compute the weights of the       ###
### Fine-Grid Spatial Interaction Model (FGSIM) ###
### Only distances                              ###
###################################################


### Function to compute raw weights based on distances
# wRaw = exp(-d*mu + d^2sigma2/2)
# nbmat:      Matrix of neighbourhood orders
# distPars:   Parameters of the normal distribution
# d:          Decay parameter
# minlag:     Minimum ngb order considered
# maxlag:     Maximum ngb order considered
# normalize:  Scale weights by row sums?

WrawDistance = function(nbmat,
                        distPars,
                        d = 1,
                        minlag = 0,
                        maxlag = Inf,
                        normalize = FALSE){
  out = exp(-d*distPars$mu + d^2*distPars$sigma2/2)
  out[nbmat > maxlag] = 0
  out[nbmat < minlag] = 0
  out[is.na(out)] = 0
  if(normalize){
    out = out / rowSums(out)
  }
  dimnames(out) = dimnames(nbmat)
  return(out)
}

######################################################
### Fine-Grid Spatial Interaction Matrix           ###
### Based on distances only                        ###
### Wrapper function that returns three functions: ###
### Weights, and its first two derivatives         ###
### and initial value for the parameter            ###
######################################################

W_pdist = function(maxlag,
                   distPars,
                   normalize = TRUE,
                   log = FALSE,
                   initial = ifelse(log, 0, 1),
                   from0 = TRUE){
  
  ###############
  ### Weights ###
  ###############
  
  # Call to compute raw weights
  weights.call = call("WrawDistance", quote(nbmat), quote(distPars), quote(d), ifelse(from0, 0, 1), maxlag, quote(FALSE))
  # Set up function
  weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                        envir=getNamespace("surveillance"))
  # Assign weights to an object Wraw
  body(weights) <- as.call(c(as.name("{"),
                             substitute(Wraw <- weights.call, list(weights.call=weights.call))
  ))
  
  # Transform decay parameter if it is in logs
  if (log){
    body(weights) <- as.call(append(as.list(body(weights)), 
                                    quote(d <- exp(d)), after = 1))
  }
  
  # Normalize weights if selected
  if(normalize){
    body(weights) <- as.call(append(as.list(body(weights)), 
                                    quote(Wraw / rowSums(Wraw))))
  }else{
    body(weights) <- as.call(append(as.list(body(weights)), 
                                    quote(Wraw)))
  }
  
  #######################
  ### Fist derivative ###
  #######################
  
  # Set up function
  dweights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                         envir=getNamespace("surveillance"))
  # Assign weights to an object Wraw
  body(dweights) <- as.call(c(as.name("{"),
                             substitute(Wraw <- weights.call, list(weights.call=weights.call))
  ))
  
  # Transform decay parameter if it is in logs
  if (log){
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                    quote(d <- exp(d)), after = 1))
  }
  
  # Normalize weights if selected
  if(normalize){
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(rs <- rowSums(Wraw))))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                    quote(W <- Wraw / rs)))
  }else{
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                    quote(W <- Wraw)))
  }
  
  # Compute scaling constants
  body(dweights) <- as.call(append(as.list(body(dweights)), 
                                   quote(Cmat <- d*distPars$sigma2 - distPars$mu)))
  
  
  # Compute derivative of row sums
  if(normalize){
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(drs <- rowSums(Wraw * Cmat))))
  }
  
  if(normalize){
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(Cmat <- Cmat - drs/rs)))
  }
  
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                   quote(deriv <- W*Cmat)))
  
  # Scale by derivative of parameter transformation
  if(log){
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(deriv <- deriv * d)))
  }
  # Return results
  body(dweights) <- as.call(append(as.list(body(dweights)), 
                                   quote(deriv)))
  
  #########################
  ### Second derivative ###
  #########################
  
  # Set up function
  d2weights = as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                         envir=getNamespace("surveillance"))
  # Assign weights to an object Wraw
  body(d2weights) <- as.call(c(as.name("{"),
                              substitute(Wraw <- weights.call, list(weights.call=weights.call))
  ))
  
  # Transform decay parameter if it is in logs
  if (log){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(d <- exp(d)), after = 1))
  }
  
  # Normalize weights if selected
  if(normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(rs <- rowSums(Wraw))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(W <- Wraw / rs)))
  }else{
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(W <- Wraw)))
  }
  
  # Compute scaling constants
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(Cmat <- d*distPars$sigma2 - distPars$mu)))
  
  
  # Compute derivatives of row sums
  if(normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(drs <- rowSums(Wraw * Cmat))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(d2rs <- rowSums(Wraw * (Cmat^2 + distPars$sigma2)))))
  }
  
  # If we use log transformation, save scaling matrix and first derivative
  if(log & normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(Cmat1 <- Cmat - drs/rs)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                     quote(deriv1 <- W*Cmat1)))
  }
  if(log & !normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                     quote(deriv1 <- W*Cmat)))
  }
  
  if(normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(Cmat <- Cmat + distPars$sigma2 - 2*Cmat*drs/rs + 2*drs^2/rs^2 - d2rs/rs)))
  }
  
  # Compute derivative (Holds true if not log transformed)
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                   quote(deriv <- W * Cmat)))
  
  # Scale by derivative of parameter transformation
  if(log){
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                     quote(deriv <- deriv * d^2 + deriv1*d)))
  }
  
  # Return results
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(deriv)))
   
  return(list(w = weights, dw = dweights, d2w = d2weights, initial = ifelse(log, 0, 1)))
}
