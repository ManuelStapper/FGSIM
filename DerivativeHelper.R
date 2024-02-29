### Functions to help with derivatives
### Starting with weights and its derivatives, we transform the weights in
### four different ways. The derivatives need to be computed for those
### transformations. We have:
### 1) Row normalisation          W_ij = W_ij/sum(W_i*)
### 2) Parameter transformation   d = exp(d)
### 3) Contact transformation     W = W * M * W^T
### 4) Scaling columns            W = W * P

### For each function, the input is similar to the output:
### Original weights, first derivative(s) and second derivative(s)

### The goal is to use these functions inside a weight function that
### is passed to the hhh4 function.

### To avoid unnecessary computations, a "depth" argument specifies which
### derivatives are needed (0 = only weights, 1 = 1st, 2 = 2nd)
### Replaced by checking if W1 and/or W2 is NA

any(is.na(list(matrix(0, nrow = 2, ncol = 2))))

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
  if(depth == 1){
    return(list(W = X, W1 = X1, W2 = NA))
  }
  
  X2 = list()
  counter = 1
  for(i in 1:k){
    for(j in i:k){
      X2[[counter]] = W2[[counter]] * d[i] * d[j]
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

##################################################################
##################################################################
### Old versions split by number of parameters below this line ###
##################################################################
##################################################################

# # Row normalisation for one parameter
# helperRN1 = function(W, W1 = NA, W2 = NA, depth = 0){
#   n = dim(W)[1]
#   RSmat = matrix(1, nrow = n, ncol = n)
#   RS = W %*% RSmat
#   X = W / RS
#   if(depth == 0){
#     return(X)
#   }
#   RS1 = W1 %*% RSmat
#   X1 = W1 / RS - W * RS1 / RS^2
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   RS2 = W2 %*% RSmat
#   RSsq = RS^2
#   X2 = W2 / RS - 2 * W1 * RS1 / RSsq - W * RS2 / RSsq + 2 * W * RS1^2 / RS^3
#   return(list(W = X, W1 = X1, W2 = X2))
# }
# 
# # Input also includes d, which is treated as exp(d)
# # In the outer functions, we transform d -> exp(d)
# helperLOG1 = function(d, W, W1 = NA, W2 = NA, depth = 0){
#   n = dim(W)[1]
#   X = W
#   if(depth == 0){
#     return(X)
#   }
#   X1 = W1 * d
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   X2 = W2 * d^2 + W1 * d
#   return(list(W = X, W1 = X1, W2 = X2))
# }
# 
# # Function for contact transformation
# # Now also includes a diagonal matrix M that gives the probability
# # that two individuals who are in the same district have contact
# # M is computes before estimation
# helperC1 = function(M, W, W1 = NA, W2 = NA, depth = 0){
#   n = dim(W)[1]
#   Wt = t(W)
#   X = W %*% M %*% Wt
#   if(depth == 0){
#     return(X)
#   }
#   
#   X1 = W1 %*% M %*% Wt
#   # Making use of the symmetry of M
#   X1 = X1 + t(X1)
#   
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   X2 = W2 %*% M %*% Wt + W1 %*% M %*% t(W1)
#   X2 = X2 + t(X2)
#   return(list(W = X, W1 = X1, W2 = X2))
# }
# 
# # Scale columns by diagonal matrix S
# helperSC1 = function(S, W, W1 = NA, W2 = NA, depth = 0){
#   n = dim(W)[1]
#   X = W %*% S
#   if(depth == 0){
#     return(X)
#   }
#   X1 = W1 %*% S
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   X2 = W2 %*% S
#   return(list(W = X, W1 = X1, W2 = X2))
# }
# 
# ### Functions for the gravity model
# 
# # Row normalisation for one parameter
# helperRN3 = function(W, W1 = NA, W2 = NA, depth = 0){
#   n = dim(W)[1]
#   RSmat = matrix(1, nrow = n, ncol = n)
#   RS = W %*% RSmat
#   X = W / RS
#   if(depth == 0){
#     return(X)
#   }
#   RS1a = W1[[1]] %*% RSmat
#   RS1b = W1[[2]] %*% RSmat
#   RS1c = W1[[3]] %*% RSmat
#   
#   X1a = W1[[1]] / RS - W * RS1a / RS^2
#   X1b = W1[[2]] / RS - W * RS1b / RS^2
#   X1c = W1[[3]] / RS - W * RS1c / RS^2
#   
#   X1 = list(X1a, X1b, X1c)
#   
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   RS2aa = W2[[1]] %*% RSmat
#   RS2ab = W2[[2]] %*% RSmat
#   RS2ac = W2[[3]] %*% RSmat
#   RS2bb = W2[[4]] %*% RSmat
#   RS2bc = W2[[5]] %*% RSmat
#   RS2cc = W2[[6]] %*% RSmat
#   
#   RSsq = RS^2
#   
#   X2aa = W2[[1]] / RS - 2 * W1[[1]] * RS1a / RSsq - W * RS2aa / RSsq + 2 * W * RS1a^2 / RS^3
#   X2ab = W2[[2]] / RS - W1[[1]] * RS1b / RSsq - W1[[2]] * RS1a / RSsq - W * RS2ab / RSsq + 2 * W * RS1a * RS1b / RS^3
#   X2ac = W2[[3]] / RS - W1[[1]] * RS1c / RSsq - W1[[3]] * RS1a / RSsq - W * RS2ac / RSsq + 2 * W * RS1a * RS1c / RS^3
#   X2bb = W2[[4]] / RS - 2 * W1[[2]] * RS1b / RSsq - W * RS2bb / RSsq + 2 * W * RS1b^2 / RS^3
#   X2bc = W2[[5]] / RS - W1[[2]] * RS1c / RSsq - W1[[3]] * RS1b / RSsq - W * RS2bc / RSsq + 2 * W * RS1b * RS1c / RS^3
#   X2cc = W2[[6]] / RS - 2 * W1[[3]] * RS1c / RSsq - W * RS2cc / RSsq + 2 * W * RS1c^2 / RS^3
#   
#   X2 = list(X2aa, X2ab, X2ac, X2bb, X2bc, X2cc)
#   
#   return(list(W = X, W1 = X1, W2 = X2))
# }
# 
# # Input also includes d, which is treated as exp(d)
# # In the outer functions, we transform d -> exp(d)
# helperLOG3 = function(d, W, W1 = NA, W2 = NA, depth = 0){
#   n = dim(W)[1]
#   X = W
#   if(depth == 0){
#     return(X)
#   }
#   X1 = W1
#   X1[[1]] = X1[[1]] * d[1]
#   X1[[2]] = X1[[2]] * d[2]
#   X1[[3]] = X1[[3]] * d[3]
#   
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   X2 = W2
#   X2[[1]] = X2[[1]] * d[1] * d[1]
#   X2[[2]] = X2[[2]] * d[1] * d[2]
#   X2[[3]] = X2[[3]] * d[1] * d[3]
#   X2[[4]] = X2[[4]] * d[2] * d[2]
#   X2[[5]] = X2[[5]] * d[2] * d[3]
#   X2[[6]] = X2[[6]] * d[3] * d[3]
#   
#   return(list(W = X, W1 = X1, W2 = X2))
# }
# 
# # Function for contact transformation
# # Now also includes a diagonal matrix M that gives the probability
# # that two individuals who are in the same district have contact
# # M is computed before estimation
# helperC3 = function(M, W, W1 = NA, W2 = NA, depth = 0){
#   n = dim(W)[1]
#   Wt = t(W)
#   X = W %*% M %*% Wt
#   if(depth == 0){
#     return(X)
#   }
#   
#   X1 = W1
#   X1[[1]] = W1[[1]] %*% M %*% Wt
#   X1[[1]] = X1[[1]] + t(X1[[1]])
#   X1[[2]] = W1[[2]] %*% M %*% Wt
#   X1[[2]] = X1[[2]] + t(X1[[2]])
#   X1[[3]] = W1[[3]] %*% M %*% Wt
#   X1[[3]] = X1[[3]] + t(X1[[3]])
#   
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   X2 = W2
#   X2[[1]] = W2[[1]] %*% M %*% Wt + W1[[1]] %*% M %*% t(W1[[1]])
#   X2[[1]] = X2[[1]] + t(X2[[1]])
#   
#   X2[[2]] = W2[[2]] %*% M %*% Wt + W1[[1]] %*% M %*% t(W1[[2]])
#   X2[[2]] = X2[[2]] + t(X2[[2]])
#   
#   X2[[3]] = W2[[3]] %*% M %*% Wt + W1[[1]] %*% M %*% t(W1[[3]])
#   X2[[3]] = X2[[3]] + t(X2[[3]])
#   
#   X2[[4]] = W2[[4]] %*% M %*% Wt + W1[[2]] %*% M %*% t(W1[[2]])
#   X2[[4]] = X2[[4]] + t(X2[[4]])
#   
#   X2[[5]] = W2[[5]] %*% M %*% Wt + W1[[2]] %*% M %*% t(W1[[3]])
#   X2[[5]] = X2[[5]] + t(X2[[5]])
#   
#   X2[[6]] = W2[[6]] %*% M %*% Wt + W1[[3]] %*% M %*% t(W1[[3]])
#   X2[[6]] = X2[[6]] + t(X2[[6]])
#   
#   return(list(W = X, W1 = X1, W2 = X2))
# }
# 
# helperSC3 = function(S, W, W1 = NA, W2 = NA, depth = 0){
#   X = W %*% S
#   if(depth == 0){
#     return(X)
#   }
#   
#   X1 = lapply(W1, function(xx) xx %*% S)
#   if(depth == 1){
#     return(list(W = X, W1 = X1))
#   }
#   
#   X2 = lapply(W2, function(xx) xx %*% S)
#   return(list(W = X, W1 = X1, W2 = X2))
# }
