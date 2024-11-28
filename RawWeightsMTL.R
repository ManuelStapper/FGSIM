# Functions to compute raw weights (no transformations)
# Three versions depending on which distribution is assumed for random distance
# measure:

# 1) Mixture of truncated LogNormals and Power Law
# 2) Mixture of truncated LogNormals and truncated Power Law
# 3) Mixture of bivariate LogNormals (for gravity model only)

# Input:
# nbmat:  Neighbourhood matrix
# pars:   Parameters of truncated normal
#         Should contain origin, destination, mu1, ..., muK,
#         sigma1, ..., sigmaK, w1, ..., wK, l, u
# d:      Decay parameter
# minlag: Minimum spatial lag considered
# maxlag: Maximum spatial lag considered
# depth:  Order of derivatives returned

# Output:
# List of weights, first and second derivatives
# Derivatives are lists even for univariate parameters
WrawMTL = function(nbmat, pars, d = 1, minlag = 0, maxlag = Inf, depth = 0){
  n = nrow(nbmat)
  K = dim(pars$mu)[3]
  
  W = matrix(0, nrow = n, ncol = n)
  if(depth > 0) W1 = matrix(0, nrow = n, ncol = n)
  if(depth > 1) W2 = matrix(0, nrow = n, ncol = n)
  
  mu = pars$mu
  sigma = pars$sigma
  w = pars$w
  l = pars$l
  u = pars$u
  
  sigma2 = sigma^2
  mu_d_sigma2 = -mu + d * sigma2
  
  for(k in 1:K){
    Z1 = (u + mu_d_sigma2[, , k]) / sigma[, , k]
    Z2 = (l + mu_d_sigma2[, , k]) / sigma[, , k]
    Z3 = (u - mu[, , k]) / sigma[, , k]
    Z4 = (l - mu[, , k]) / sigma[, , k]
    
    P1 = pnorm(Z1)
    P2 = pnorm(Z2)
    P3 = pnorm(Z3)
    P4 = pnorm(Z4)
    D1 = dnorm(Z1)
    D2 = dnorm(Z2)
    
    hInv = 1 / (P3 - P4)
    f = exp(-d * mu[, , k] + d^2 * sigma2[, , k] / 2)
    q = P1 - P2
    
    W = W + w[, , k] * f * hInv * q
    
    if(depth > 0){
      q1 = (D1 - D2) * sigma[, , k]
      f1 = f * mu_d_sigma2[, , k]
      W1 = W1 + w[, , k] * (f1 * q + f * q1) * hInv
    }
    
    if(depth > 1){
      q2 = (D2 * (l + mu_d_sigma2[, , k]) - D1 * (u + mu_d_sigma2[, , k])) * sigma[, , k]
      f2 = f * (mu_d_sigma2[, , k]^2 + sigma2[, , k])
      W2 = W2 + w[, , k] * (f2 * q + 2 * f1 * q1 + f * q2) * hInv
    }
  }
  
  W[nbmat > maxlag] = 0
  W[nbmat < minlag] = 0
  W[is.na(W)] = 0
  dimnames(W) = dimnames(nbmat)
  
  if(depth == 0){
    return(list(W = W, W1 = NA, W2 = NA))
  }
  
  if(depth == 1){
    W1[nbmat > maxlag] = 0
    W1[nbmat < minlag] = 0
    W1[is.na(W1)] = 0
    dimnames(W1) = dimnames(nbmat)
    
    return(list(W = W, W1 = list(W1), W2 = NA))
  }
  
  if(depth == 2){
    W2[nbmat > maxlag] = 0
    W2[nbmat < minlag] = 0
    W2[is.na(W2)] = 0
    dimnames(W2) = dimnames(nbmat)
    
    return(list(W = W, W1 = list(W1), W2 = list(W2)))
  }
}

# Function for truncated power law
WrawMTLtrunc = function(nbmat, pars, d = c(1, 1), minlag = 0, maxlag = Inf, depth = 0){
  n = nrow(nbmat)
  K = dim(pars$mu)[3]
  
  W = matrix(0, nrow = n, ncol = n)
  if(depth > 0) W1_1 = W1_2 = matrix(0, nrow = n, ncol = n)
  if(depth > 1) W2_11 = W2_12 = W2_22 = matrix(0, nrow = n, ncol = n)
  
  delta = d[2] # Caution: This is in logs
  exp_delta = exp(delta)
  d = d[1]
  
  mu = pars$mu
  sigma = pars$sigma
  w = pars$w
  l = pars$l
  u = pars$u
  
  sigma2 = sigma^2
  sigma3inv = 1/(sigma^3)
  mu_d_sigma2 = -mu + d * sigma2
  mu_d2_sigma2 = -mu + d^2 * sigma2
  
  
  U1 = pmin(pmax(l, delta), u)
  U1_2 = (delta >= l)*(delta <= u)
  
  for(k in 1:K) {
    Z1 = (U1 - mu[, , k])/sigma[, , k]
    Z2 = (l - mu[, , k])/sigma[, , k]
    Z3 = (u + mu_d_sigma2[, , k])/sigma[, , k]
    Z4 = (U1 + mu_d_sigma2[, , k])/sigma[, , k]
    Z5 = (u - mu[, , k])/sigma[, , k]
    
    P1 = pnorm(Z1)
    P2 = pnorm(Z2)
    P3 = pnorm(Z3)
    P4 = pnorm(Z4)
    P5 = pnorm(Z5)
    
    p = P1 - P2
    q = P3 - P4
    
    f = exp(-d*delta)
    g = exp(-d*mu[, , k] + d^2*sigma2[, , k]/2)
    
    hInv = 1 / (P5 - P2)
    
    W = W + w[, , k] * (p * f + q * g) * hInv
    
    if(depth > 0) {
      D1 = dnorm(Z1)
      D2 = dnorm(Z2)
      D3 = dnorm(Z3)
      D4 = dnorm(Z4)
      
      f1 = f*(-delta)
      q1 = (D3 - D4)*sigma[, , k]
      g1 = g * mu_d_sigma2[, , k]
      
      f2 = f*(-d)
      p2 = D1 * U1_2 / sigma[, , k]
      q2 = -D4 * U1_2 / sigma[, , k]
      
      W1_1 = W1_1 + w[, , k] * (p * f1 + q1 * g + q*g1) * hInv
      W1_2 = W1_2 + w[, , k] * (p2 * f + p * f2 + q2 * g) * hInv
    }
    
    if(depth > 1) {
      f11 = delta^2 * f
      q11 = (D4 * (U1 + mu_d_sigma2[, , k]) - D3 * (u + mu_d_sigma2[, , k])) * sigma[, , k]
      g11 = g * (mu_d_sigma2[, , k]^2 + sigma2[, , k])
      
      p22 = D1 * (mu[, , k] - U1) * U1_2 * sigma3inv[, , k]
      f22 = d^2 * f
      q22 = D4 * ((U1 + mu_d_sigma2[, , k])*U1_2) * sigma3inv[, , k]
      
      f12 = (delta*d - 1) * f
      q12 = D4 * U1_2 * (U1 + mu_d_sigma2[, , k])/sigma[, , k]
      
      W2_11 = W2_11 + w[, , k] * (p * f11 + q11 * g + 2 * q1 * g1 + q * g11) * hInv
      W2_12 = W2_12 + w[, , k] * (p2 * f1 + p * f12 + q12 * g + q2 * g1) * hInv
      W2_22 = W2_22 + w[, , k] * (p22 * f + 2 * p2 * f2 + p * f22 + q22 * g) * hInv
    }
  }
  
  W[nbmat > maxlag] = 0
  W[nbmat < minlag] = 0
  W[is.na(W)] = 0
  dimnames(W) = dimnames(nbmat)
  
  if(depth == 0) {
    return(list(W = W, W1 = NA, W2 = NA))
  }
  
  if(depth == 1) {
    W1_1[nbmat > maxlag] = 0
    W1_1[nbmat < minlag] = 0
    W1_1[is.na(W1_1)] = 0
    dimnames(W1_1) = dimnames(nbmat)
    W1_2[nbmat > maxlag] = 0
    W1_2[nbmat < minlag] = 0
    W1_2[is.na(W1_2)] = 0
    dimnames(W1_2) = dimnames(nbmat)
    
    return(list(W = W, W1 = list(W1_1, W1_2), W2 = NA))
  }
  
  if(depth == 2) {
    W2_11[nbmat > maxlag] = 0
    W2_11[nbmat < minlag] = 0
    W2_11[is.na(W2_11)] = 0
    dimnames(W2_11) = dimnames(nbmat)
    W2_12[nbmat > maxlag] = 0
    W2_12[nbmat < minlag] = 0
    W2_12[is.na(W2_12)] = 0
    dimnames(W2_22) = dimnames(nbmat)
    W2_22[nbmat > maxlag] = 0
    W2_22[nbmat < minlag] = 0
    W2_22[is.na(W2_22)] = 0
    dimnames(W2_22) = dimnames(nbmat)
    
    return(list(W = W, W1 = list(W1_1, W1_2), W2 = list(W2_11, W2_12, W2_22)))
  }
}


# Function for mixture of bivariate normal distribution

# Structure of the input list "pars":
# $mu1, $mu2, $sigma1, $sigma12, $sigma2, $w (all 3D arrays)
# mu1 corresponds to distance, mu2 to population ratio

WrawMN = function(nbmat, pars, d = c(1, 1), minlag = 0, maxlag = Inf, depth = 0){
  n = nrow(nbmat)
  K = dim(pars$mu1)[3]
  
  W = matrix(0, nrow = n, ncol = n)
  if(depth > 0) W1_1 = W1_2 = matrix(0, nrow = n, ncol = n)
  if(depth > 1) W2_11 = W2_12 = W2_22 = matrix(0, nrow = n, ncol = n)
  
  alpha = d[2]
  d = d[1]
  
  mu1 = pars$mu1
  mu2 = pars$mu2
  sigma1 = pars$sigma1
  sigma12 = pars$sigma12
  sigma2 = pars$sigma2
  w = pars$w
  
  for(k in 1:K){
    Wtemp = w[, , k] * exp(alpha * mu2[, , k] - d*mu1[, , k] + alpha^2*sigma2[, , k]/2 - alpha*d*sigma12[, , k] + d^2*sigma1[, , k]/2)
    W = W + Wtemp
    
    if(depth > 0) {
      c1 = (-mu1[, , k] - alpha*sigma12[, , k] + d*sigma1[, , k])
      c2 = (mu2[, , k] - d*sigma12[, , k] + alpha*sigma2[, , k])
      W1_1 = W1_1 + Wtemp * c1
      W1_2 = W1_2 + Wtemp * c2
    }
    
    if(depth > 1) {
      W2_11 = W2_11 + Wtemp*(c1^2 + sigma1[, , k])
      W2_12 = W2_12 + Wtemp*(c1*c2 - sigma12[, , k])
      W2_22 = W2_22 + Wtemp*(c2^2 + sigma2[, , k])
    }
  }
  
  W[nbmat > maxlag] = 0
  W[nbmat < minlag] = 0
  W[is.na(W)] = 0
  dimnames(W) = dimnames(nbmat)
  
  if(depth == 0) {
    return(list(W = W, W1 = NA, W2 = NA))
  }
  
  if(depth == 1) {
    W1_1[nbmat > maxlag] = 0
    W1_1[nbmat < minlag] = 0
    W1_1[is.na(W1_1)] = 0
    dimnames(W1_1) = dimnames(nbmat)
    W1_2[nbmat > maxlag] = 0
    W1_2[nbmat < minlag] = 0
    W1_2[is.na(W1_2)] = 0
    dimnames(W1_2) = dimnames(nbmat)
    
    return(list(W = W, W1 = list(W1_1, W1_2), W2 = NA))
  }
  
  if(depth == 2) {
    W2_11[nbmat > maxlag] = 0
    W2_11[nbmat < minlag] = 0
    W2_11[is.na(W2_11)] = 0
    dimnames(W2_11) = dimnames(nbmat)
    W2_12[nbmat > maxlag] = 0
    W2_12[nbmat < minlag] = 0
    W2_12[is.na(W2_12)] = 0
    dimnames(W2_22) = dimnames(nbmat)
    W2_22[nbmat > maxlag] = 0
    W2_22[nbmat < minlag] = 0
    W2_22[is.na(W2_22)] = 0
    dimnames(W2_22) = dimnames(nbmat)
    
    return(list(W = W, W1 = list(W1_1, W1_2), W2 = list(W2_11, W2_12, W2_22)))
  }
}










