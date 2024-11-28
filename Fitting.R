# Fitting (truncated) Normal distributions or mixtures
library(truncnorm)
library(mvtnorm)

# Four settings:
# - Normal distribution
# - Truncated Normal distribution
# - Mixture of Normal distributions
# - Mixture of truncated Normal distributions
# - Mixture of truncated Normals where sample is weighted
# - Mixture of bivariate Normals (not truncated)

# Estimates parameters and returns estimates + BIC
# Estimation functions work on logs of distance measures

# Simple Normal distribution
estN = function(x){
  ll = sum(dnorm(x, mean(x), s = sd(x), log = T))
  return(list(m = mean(x), s = sd(x), aic = 4 - 2*ll, bic = 2*log(length(x)) - 2*ll))
}

# Truncated Normal distribution (between l und u)
estTN = function(x, l, u){
  tf = function(theta, x, l, u){
    m = theta[1]
    s = theta[2]
    if(s <= 0) return(Inf)
    -sum(log(dtruncnorm(x = x, a = l, b = u, mean = m, sd = s)))
  }
  optRaw = optim(c(mean(x), sd(x)), function(vars) tf(vars, x, l, u))
  est = optRaw$par
  ll = -optRaw$value
  
  return(list(m = est[1], s = est[2], aic = 4 - 2*ll, bic = 2*log(length(x)) - 2*ll))
}

# Mixture of Normal distributions
# k:     number of mixture components
# nIter: number of iterations for EM algorithm

estMN = function(x, k, nIter = 100){
  # Initialisation
  sx = sort(x)
  n = length(x)
  temp = round(seq(1, length(x), length.out = k+1))
  m = sapply(1:k, function(i) mean(sx[temp[i]:temp[i+1]]))
  s = sapply(1:k, function(i) sd(sx[temp[i]:temp[i+1]]))
  w = rep(1/k, k)
  
  # Iteration
  for(iter in 1:nIter){
    f = sapply(1:k, function(i) w[i]*dnorm(x, m[i], s[i]))
    f = f / rowSums(f)
    w = colMeans(f)
    
    m = c((t(f)%*%x)/colSums(f))
    s = sqrt(sapply(1:k, function(i) sum(f[, i]*(x - m[i])^2)/sum(f[, i])))
  }
  
  ll = sum(log(rowSums(sapply(1:k, function(i) w[i]*dnorm(x, m[i], s[i])))))
  
  return(list(m = m, s = s, w = w, aic = 2*k*3 - 2*ll, bic = k*3*log(length(x)) - 2*ll))
}


# Mixture of truncated Normal distribution
# k:     Number of components
# l:     Lower bound (same for all components)
# u:     Upper bound (same for all components)
# nIter: Number of EM iterations
estMTN = function(x, k, l, u, nIter = 100){
  # Initialisation
  sx = sort(x)
  n = length(x)
  temp = round(seq(1, length(x), length.out = k+1))
  m = sapply(1:k, function(i) mean(sx[temp[i]:temp[i+1]]))
  s = sapply(1:k, function(i) sd(sx[temp[i]:temp[i+1]]))
  w = rep(1/k, k)
  
  # Iteration
  for(iter in 1:nIter){
    scl = w/(pnorm(u, m, s) - pnorm(l, m, s))
    # f = sapply(1:k, function(i) w[i]*dtruncnorm(x, a = l, b = u, mean = m[i], sd = s[i]))
    f = sapply(1:k, function(i) scl[i]*dnorm(x, mean = m[i], sd = s[i]))
    f = f / rowSums(f)
    w = colMeans(f)
    
    m = c((t(f)%*%x)/colSums(f))
    s = sqrt(sapply(1:k, function(i) sum(f[, i]*(x - m[i])^2)/sum(f[, i])))
  }
  
  scl = w/(pnorm(u, m, s) - pnorm(l, m, s))
  
  ll = sum(log(rowSums(sapply(1:k, function(i){
    scl[i]*dnorm(x, mean = m[i], sd = s[i])
  }))))
  
  return(list(m = m, s = s, w = w, aic = 2*k*3 - 2*ll, bic = k*3*log(length(x)) - 2*ll))
}


# Adding a weighted estimation here to be used for travel time data
# Still a bit slowish

# Input: x:  distinct observations
#        nr: number of observations with certain value
estMTNweighted = function(xUni, nr, k, l, u, nIter = 100) {
  x = rep(xUni, nr)
  # Initialisation
  sx = sort(x)
  n = length(x)
  temp = round(seq(1, length(x), length.out = k+1))
  m = sapply(1:k, function(i) mean(sx[temp[i]:temp[i+1]]))
  s = sapply(1:k, function(i) sd(sx[temp[i]:temp[i+1]]))
  w = rep(1/k, k)
  
  # Iteration
  for(iter in 1:nIter){
    scl = w/(pnorm(u, m, s) - pnorm(l, m, s))
    # f = sapply(1:k, function(i) w[i]*dtruncnorm(x, a = l, b = u, mean = m[i], sd = s[i]))
    f = sapply(1:k, function(i) scl[i]*dnorm(x, mean = m[i], sd = s[i]))
    f = f / rowSums(f)
    w = colMeans(f)
    
    m = c((t(f)%*%x)/colSums(f))
    s = sqrt(sapply(1:k, function(i) sum(f[, i]*(x - m[i])^2)/sum(f[, i])))
  }
  
  scl = w/(pnorm(u, m, s) - pnorm(l, m, s))
  
  ll = sum(log(rowSums(sapply(1:k, function(i){
    scl[i]*dnorm(x, mean = m[i], sd = s[i])
  }))))
  
  return(list(m = m, s = s, w = w, aic = 2*k*3 - 2*ll, bic = k*3*log(length(x)) - 2*ll))
}

# Fitting a mixture of bivariate Normals
# x: Matrix with 2 columns
# k: Number of components
estMN2 = function(x, k, nIter = 100){
  x1 = x[, 1]
  x2 = x[, 2]
  
  if(k == 1){
    m = matrix(c(mean(x1), mean(x2)), nrow = 1)
    S = matrix(c(var(x1), cov(x1, x2), cov(x1, x2), var(x2)), nrow = 1)
    w = c(1)
    ll = sum(log(w[1]*dmvnorm(x, m[1, ], matrix(S[1, ], nrow = 2))))
    return(list(m = m, S = S, w = w, aic = 2*5 - 2*ll, bic = 5*log(2*length(x1)) - 2*ll))
  }
  
  # Initialisation
  sx1 = sort(x1)
  sx2 = sort(x2)
  n = length(x1)
  temp = round(seq(1, n, length.out = k+1))
  m1 = sapply(1:k, function(i) mean(sx1[temp[i]:temp[i+1]]))
  s1 = sapply(1:k, function(i) var(sx1[temp[i]:temp[i+1]]))
  s12 = sapply(1:k, function(i) cov(sx1[temp[i]:temp[i+1]], sx2[temp[i]:temp[i+1]]))
  m2 = sapply(1:k, function(i) mean(sx2[temp[i]:temp[i+1]]))
  s2 = sapply(1:k, function(i) var(sx2[temp[i]:temp[i+1]]))
  s12 = sapply(1:k, function(i) cov(sx1[temp[i]:temp[i+1]], sx2[temp[i]:temp[i+1]]))
  w = rep(1/k, k)
  
  m = cbind(m1, m2)
  S = cbind(s1, s12, s12, s2)
  
  # Iteration
  for(iter in 1:nIter){
    f = sapply(1:k, function(i) w[i]*dmvnorm(x, m[i, ], matrix(S[i, ], nrow = 2)))
    f = f / rowSums(f)
    f = t(apply(f, 1, function(xx) if(any(is.na(xx))) rep(1/k, k) else xx))
    w = colMeans(f)
    
    m[, 1] = c((t(f)%*%x1)/colSums(f))
    m[, 2] = c((t(f)%*%x2)/colSums(f))
    
    S[, 1] = sapply(1:k, function(i) sum(f[, i]*(x1 - m[i, 1])^2)/sum(f[, i]))
    S[, 2] = sapply(1:k, function(i) sum(f[, i]*(x1 - m[i, 1])*(x2 - m[i, 2]))/sum(f[, i]))
    S[, 3] = S[, 2]
    S[, 4] = sapply(1:k, function(i) sum(f[, i]*(x2 - m[i, 2])^2)/sum(f[, i]))
  }
  
  ll = sum(log(rowSums(sapply(1:k, function(i) w[i]*dmvnorm(x, m[i, ], matrix(S[i, ], nrow = 2))))))
  
  return(list(m = m, S = S, w = w, aic = 2*k*5 - 2*ll, bic = k*5*log(2*n) - 2*ll))
}
