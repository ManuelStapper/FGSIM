###################################################
### Methods to compute the weights of the       ###
### Fine-Grid Spatial Interaction Model (FGSIM) ###
### For gravity model weights                   ###
###################################################


### Function to compute raw weights including gravity component
# nbmat:         Matrix of neighbourhood orders
# gravityPars:   Parameters of the normal distribution
# d:             Decay parameters
# minlag:        Minimum ngb order considered
# maxlag:        Maximum ngb order considered
# normalize:     Scale weights by row sums?

WrawGravity = function(nbmat,
                       gravityPars,
                       d = c(1, 1, 1),
                       minlag = 0,
                       maxlag = Inf,
                       normalize = FALSE){
  out = -d[1]*gravityPars$dist
  out = out - d[2]*gravityPars$pO
  out = out + d[3]*gravityPars$pD
  
  out = out + d[1]^2*gravityPars$s11/2
  out = out + d[2]^2*gravityPars$s22/2
  out = out + d[3]^2*gravityPars$s33/2
  
  out = out + d[1]*d[2]*gravityPars$s12
  out = out - d[1]*d[3]*gravityPars$s13
  out = out - d[2]*d[3]*gravityPars$s23
  
  out = exp(out)
  
  out[is.na(out)] = 0
  
  out[nbmat > maxlag] = 0
  out[nbmat < minlag] = 0
  if(normalize){
    out = out / rowSums(out)
  }
  dimnames(out) = dimnames(nbmat)
  return(out)
}

######################################################
### Fine-Grid Spatial Interaction Matrix           ###
### Including gravity model component              ###
### Wrapper function that returns three functions: ###
### Weights, and its first two derivatives         ###
### and initial values for parameters              ###
######################################################

W_gravity = function(maxlag,
                     gravityPars,
                     normalize = TRUE,
                     log = FALSE,
                     initial = ifelse(log, c(0, 0, 0), c(1, 1, 1)),
                     from0 = TRUE){
  
  ###############
  ### Weights ###
  ###############
  
  # Call to compute raw weights
  weights.call = call("WrawGravity", quote(nbmat), quote(gravityPars), quote(d), ifelse(from0, 0, 1), maxlag, quote(FALSE))
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
                                   quote(CmatA <- -gravityPars$dist + d[1]*gravityPars$s11 + d[2]*gravityPars$s12 - d[3]*gravityPars$s13)))
  body(dweights) <- as.call(append(as.list(body(dweights)), 
                                   quote(CmatB <- -gravityPars$pO + d[2]*gravityPars$s22 + d[1]*gravityPars$s12 - d[3]*gravityPars$s23)))
  body(dweights) <- as.call(append(as.list(body(dweights)), 
                                   quote(CmatC <- gravityPars$pD + d[3]*gravityPars$s33 - d[1]*gravityPars$s13 - d[2]*gravityPars$s23)))
  
  
  # Compute derivative of row sums
  if(normalize){
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(drsA <- rowSums(Wraw * CmatA))))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(drsB <- rowSums(Wraw * CmatB))))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(drsC <- rowSums(Wraw * CmatC))))
  }
  
  if(normalize){
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(CmatA <- CmatA - drsA/rs)))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(CmatB <- CmatB - drsB/rs)))
    body(dweights) <- as.call(append(as.list(body(dweights)), 
                                     quote(CmatC <- CmatC - drsC/rs)))
    
  }
  
  body(dweights) <- as.call(append(as.list(body(dweights)),
                                   quote(deriv <- list(W*CmatA, W*CmatB, W*CmatC))))
  
  # Scale by derivative of parameter transformation
  if(log){
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(deriv[[1]] <- deriv[[1]] * d[1])))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(deriv[[2]] <- deriv[[2]] * d[2])))
    body(dweights) <- as.call(append(as.list(body(dweights)),
                                     quote(deriv[[3]] <- deriv[[3]] * d[3])))
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
                                   quote(CmatA <- -gravityPars$dist + d[1]*gravityPars$s11 + d[2]*gravityPars$s12 - d[3]*gravityPars$s13)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatB <- -gravityPars$pO + d[2]*gravityPars$s22 + d[1]*gravityPars$s12 - d[3]*gravityPars$s23)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatC <- gravityPars$pD + d[3]*gravityPars$s33 - d[1]*gravityPars$s13 - d[2]*gravityPars$s23)))
  
  # Replaced here
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatAA <- gravityPars$s11)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatBB <- gravityPars$s22)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatCC <- gravityPars$s33)))
  
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatAB <- gravityPars$s12)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatAC <- - gravityPars$s13)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                   quote(CmatBC <- - gravityPars$s23)))
  
  
  # Compute derivatives of row sums
  if(normalize){
    # First derivatives
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(drsA <- rowSums(Wraw * CmatA))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(drsB <- rowSums(Wraw * CmatB))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(drsC <- rowSums(Wraw * CmatC))))
    
    # Second and cross derivatives
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(d2rsAA <- rowSums(Wraw * (CmatA^2 + CmatAA)))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(d2rsBB <- rowSums(Wraw * (CmatB^2 + CmatBB)))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(d2rsCC <- rowSums(Wraw * (CmatC^2 + CmatCC)))))
    
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(d2rsAB <- rowSums(Wraw * (CmatA * CmatB + CmatAB)))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(d2rsAC <- rowSums(Wraw * (CmatA * CmatC + CmatAC)))))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(d2rsBC <- rowSums(Wraw * (CmatB * CmatC + CmatBC)))))
    
  }
  
  # If we use log transformation, save scaling matrix and first derivative
  if(log & normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(CmatA1 <- CmatA - drsA/rs)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(CmatB1 <- CmatB - drsB/rs)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(CmatC1 <- CmatC - drsC/rs)))
    
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(deriv1A <- W*CmatA1)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(deriv1B <- W*CmatB1)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                     quote(deriv1C <- W*CmatC1)))
    
  }
  
  if(log & !normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(deriv1A <- W*CmatA)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(deriv1B <- W*CmatB)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(deriv1C <- W*CmatC)))
  }
  
  if(normalize){
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(CmatAA <- (CmatA - drsA/rs)^2 + CmatAA - d2rsAA/rs + drsA^2/rs^2)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(CmatBB <- (CmatB - drsB/rs)^2 + CmatBB - d2rsBB/rs + drsB^2/rs^2)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(CmatCC <- (CmatC - drsC/rs)^2 + CmatCC - d2rsCC/rs + drsC^2/rs^2)))
    
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(CmatAB <- (CmatA - drsA/rs)*(CmatB - drsB/rs) + CmatAB - d2rsAB/rs + drsA*drsB/rs^2)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(CmatAC <- (CmatA - drsA/rs)*(CmatC - drsC/rs) + CmatAC - d2rsAC/rs + drsA*drsC/rs^2)))
    body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                      quote(CmatBC <- (CmatB - drsB/rs)*(CmatC - drsC/rs) + CmatBC - d2rsBC/rs + drsB*drsC/rs^2)))
  }
  
  # Compute derivative (Holds true if not log transformed)
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(derivAA <- W * CmatAA)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(derivBB <- W * CmatBB)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(derivCC <- W * CmatCC)))
  
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(derivAB <- W * CmatAB)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(derivAC <- W * CmatAC)))
  body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                    quote(derivBC <- W * CmatBC)))
  
  # Scale by derivative of parameter transformation
  if(log){
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(derivAA <- derivAA * d[1]^2 + deriv1A*d[1])))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(derivBB <- derivBB * d[2]^2 + deriv1B*d[2])))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(derivCC <- derivCC * d[3]^2 + deriv1C*d[3])))
    
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(derivAB <- derivAB * d[1]*d[2])))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(derivAC <- derivAC * d[1]*d[3])))
    body(d2weights) <- as.call(append(as.list(body(d2weights)),
                                      quote(derivBC <- derivBC * d[2]*d[3])))
  }
  
  # Return results
  body(d2weights) <- as.call(append(as.list(body(d2weights)), 
                                    quote(list(derivAA, derivAB, derivAC, derivBB, derivBC, derivCC))))
  
  return(list(w = weights, dw = dweights, d2w = d2weights, initial = if(log) c(0, 0, 0) else c(1, 1, 1)))
}
