### Functions that compute raw weights and its derivatives
### No normalisation or log-transformation (yet)
### Will be used as a starting point for other functions, thus
### matches the notation and structure of the functions in
### "derivativeHelper.R"

### Approach based only on distances

WrawDistance = function(nbmat,
                        distPars,
                        d = 1,
                        minlag = 0,
                        maxlag = Inf,
                        depth = 0){
  # Compute weights
  W = exp(-d*distPars$dist + d^2*distPars$s11/2)
  W[nbmat > maxlag] = 0
  W[nbmat < minlag] = 0
  W[is.na(W)] = 0
  dimnames(W) = dimnames(nbmat)
  if(depth == 0){
    return(list(W = W, W1 = NA, W2 = NA))  
  }
  W1 = W * (d*distPars$s11 - distPars$dist)
  if(depth == 1){
    return(list(W = W, W1 = list(W1), W2 = NA))
  }
  W2 = W * ((d*distPars$s11 - distPars$dist)^2 + distPars$s11)
  return(list(W = W, W1 = list(W1), W2 = list(W2)))
}

# Approach including gravity component

WrawGravity = function(nbmat,
                       gravityPars,
                       d = c(1, 1, 1),
                       minlag = 0,
                       maxlag = Inf,
                       depth = 0){
  W = -d[1]*gravityPars$dist
  W = W - d[2]*gravityPars$pO
  W = W + d[3]*gravityPars$pD
  
  W = W + d[1]^2*gravityPars$s11/2
  W = W + d[2]^2*gravityPars$s22/2
  W = W + d[3]^2*gravityPars$s33/2
  
  W = W + d[1]*d[2]*gravityPars$s12
  W = W - d[1]*d[3]*gravityPars$s13
  W = W - d[2]*d[3]*gravityPars$s23
  
  W = exp(W)
  
  W[is.na(W)] = 0
  W[nbmat > maxlag] = 0
  W[nbmat < minlag] = 0
  dimnames(W) = dimnames(nbmat)
  if(depth == 0){
    return(list(W = W, W1 = NA, W2 = NA))
  }
  Ca = -gravityPars$dist + d[1]*gravityPars$s11 + d[2] * gravityPars$s12 - d[3]*gravityPars$s13
  Cb = -gravityPars$pO + d[2]*gravityPars$s22 + d[1] * gravityPars$s12 - d[3]*gravityPars$s23
  Cc = gravityPars$pD + d[3]*gravityPars$s33 - d[1] * gravityPars$s13 - d[2]*gravityPars$s23
  
  W1a = W * Ca
  W1b = W * Cb
  W1c = W * Cc
  
  if(depth == 1){
    return(list(W = W, W1 = list(W1a, W1b, W1c), W2 = NA))
  }
  
  W2aa = W * (Ca^2 + gravityPars$s11)
  W2bb = W * (Cb^2 + gravityPars$s22)
  W2cc = W * (Cc^2 + gravityPars$s33)
  
  W2ab = W * (Ca * Cb + gravityPars$s12)
  W2ac = W * (Ca * Cc - gravityPars$s13)
  W2bc = W * (Cb * Cc - gravityPars$s23)
  
  W1 = list(W1a, W1b, W1c)
  W2 = list(W2aa, W2ab, W2ac, W2bb, W2bc, W2cc)
  
  return(list(W = W, W1 = W1, W2 = W2))
}

# Powerlaw 

WrawPowerlaw = function(nbmat,
                        d = 1,
                        minlag = 0,
                        maxlag = Inf,
                        depth = 0){
  W = nbmat^(-d)
  
  W[is.na(W)] = 0
  W[nbmat > maxlag] = 0
  W[nbmat < minlag] = 0
  dimnames(W) = dimnames(nbmat)
  if(depth == 0){
    return(list(W = W, W1 = NA, W2 = NA))
  }
  logo = log(nbmat)
  W1 = W*(-logo)
  
  if(depth == 1){
    return(list(W = W, W1 = W1, W2 = NA))
  }
  
  W2 = W*logo^2
  
  return(list(W = W, W1 = list(W1), W2 = list(W2)))
}
