### Functions to help with derivatives
### Starting with weights and their derivatives, we transform the weights in
### different ways. The derivatives need to be computed for those transformations.


# RN  = Row normalisation: W_ji --> W_ji / sum_k(W_jk)
# LOG = Log-transformations of parameters: d --> exp(d)
# C   = Quadratic form: W * M * W' for matrix M
# T   = Transpose: W --> W'
# SYM = Ensure symmetry: W --> W + W'
# SR  = Scale rows: W_ji = W_ji * S_j 
# SC  = Scale columns: W_ji = W_ji * S_i


### For each function, the input is similar to the output:
### Original weights, first derivative(s) and second derivative(s)

### The goal is to use these functions inside a weight function that
### is passed to the hhh4 function.

### To avoid unnecessary computations, a "depth" argument specifies which
### derivatives are needed (0 = only weights, 1 = 1st, 2 = 2nd)
### Replaced by checking if W1 and/or W2 is NA

# Row normalisation for one parameter
helperRN = function(W, W1 = NA, W2 = NA){
  n = dim(W)[1]
  RSmat = matrix(1, nrow = n, ncol = n)
  RS = W %*% RSmat
  X = W / RS
  
  if(any(is.na(W1))){
    return(list(W = X, W1 = NA, W2 = NA))
  }
  
  k = length(W1)
  RS1 = lapply(W1, function(xx) xx %*% RSmat)
  X1 = lapply(1:k, function(i) W1[[i]] / RS - W * RS1[[i]] / RS^2)
  
  if(any(is.na(W2))){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  RS2 = lapply(W2, function(xx) xx %*% RSmat)
  RSsq = RS^2
  
  counter = 1
  X2 = list()
  for(i in 1:k){
    for(j in i:k){
      if(i == j){
        X2[[counter]] = W2[[counter]] / RS - 2 * W1[[i]] * RS1[[i]] / RSsq - W * RS2[[counter]] / RSsq + 2 * W * RS1[[i]]^2 / RS^3
      }else{
        X2[[counter]] = W2[[counter]] / RS - W1[[i]] * RS1[[j]] / RSsq - W1[[j]] * RS1[[i]] / RSsq - W * RS2[[counter]] / RSsq + 2 * W * RS1[[i]] * RS1[[j]] / RS^3
      }
      counter = counter + 1
    }
  }
  
  return(list(W = X, W1 = X1, W2 = X2))
}

# Input also includes d, which is treated as exp(d)
# In the outer functions, we transform d -> exp(d)
helperLOG = function(d, W, W1 = NA, W2 = NA){
  n = dim(W)[1]
  X = W
  if(any(is.na(W1))){
    return(list(W = X, W1 = NA, W2 = NA))
  }
  
  k = length(W1)
  X1 = lapply(1:k, function(i) W1[[i]] * d[i])
  if(any(is.na(W2))){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  X2 = list()
  counter = 1
  for(i in 1:k){
    for(j in i:k){
      if(i == j){
        X2[[counter]] = W2[[counter]] * d[i] * d[j] + W1[[i]] * d[i]
      }else{
        X2[[counter]] = W2[[counter]] * d[j] * d[i]
      }
      
      counter = counter + 1
    }
  }

  return(list(W = X, W1 = X1, W2 = X2))
}

# Function for contact transformation
# Now also includes a diagonal matrix M that gives the probability
# that two individuals who are in the same district have contact
# M is computed before estimation
helperC = function(M, W, W1 = NA, W2 = NA){
  n = dim(W)[1]
  Wt = t(W)
  X = W %*% M %*% Wt
  if(any(is.na(W1))){
    return(list(W = X, W1 = NA, W2 = NA))
  }
  
  k = length(W1)
  X1 = lapply(1:k, function(i) W1[[i]] %*% M %*% Wt)
  X1 = lapply(X1, function(xx) xx + t(xx))
  
  if(any(is.na(W2))){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  X2 = list()
  counter = 1
  for(i in 1:k){
    for(j in i:k){
      X2[[counter]] = W2[[counter]] %*% M %*% Wt + W1[[i]] %*% M %*% t(W1[[j]])
      X2[[counter]] = X2[[counter]] + t(X2[[counter]])
      counter = counter + 1
    }
  }
  
  return(list(W = X, W1 = X1, W2 = X2))
}

helperT = function(W, W1 = NA, W2 = NA){
  X = t(W)
  if(any(is.na(W1))){
    return(list(W = X, W1 = NA, W2 = NA))
  }
  
  X1 = lapply(W1, function(xx) t(xx))
  if(any(is.na(W2))){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  X2 = lapply(W2, function(xx) t(xx))
  return(list(W = X, W1 = X1, W2 = X2))
}

helperSYM = function(W, W1 = NA, W2 = NA){
  X = W + t(W)
  if(any(is.na(W1))){
    return(list(W = X, W1 = NA, W2 = NA))
  }
  
  X1 = lapply(W1, function(xx) xx + t(xx))
  if(any(is.na(W2))){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  X2 = lapply(W2, function(xx) xx + t(xx))
  return(list(W = X, W1 = X1, W2 = X2))
}

helperSC = function(S, W, W1 = NA, W2 = NA){
  X = W %*% S
  if(any(is.na(W1))){
    return(list(W = X, W1 = NA, W2 = NA))
  }
  
  X1 = lapply(W1, function(xx) xx %*% S)
  if(any(is.na(W2))){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  X2 = lapply(W2, function(xx) xx %*% S)
  return(list(W = X, W1 = X1, W2 = X2))
}

helperSR = function(S, W, W1 = NA, W2 = NA){
  X = S %*% W
  if(any(is.na(W1))){
    return(list(W = X, W1 = NA, W2 = NA))
  }
  
  X1 = lapply(W1, function(xx) S %*% xx)
  if(any(is.na(W2))){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  X2 = lapply(W2, function(xx) S %*% xx)
  return(list(W = X, W1 = X1, W2 = X2))
}