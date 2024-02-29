# With samples, estimate (Log)Normal parameters
# Returns parameters of the triple
# (distance, population origin, population destination)
# Or only for distance if gravity not included

fitLN = function(sam1,                     # Sample origin district
                 sam2,                     # Sample destination district
                 i,                        # District origin
                 j,                        # District destination
                 distr,                    # Spatial object of districts
                 pop,                      # JSON file of populations
                 gravity = FALSE){         # Returns trivariate parameters
  # Sample size
  n = nrow(sam1)
  # Crop population data to districts
  pop1 <- crop(pop, distr[i])
  pop2 <- crop(pop, distr[j])
  
  # Compute population at points
  p1 = extract(pop1, sam1)
  p2 = extract(pop2, sam2)
  
  # Compute log distance (in km)
  Draw = distGeo(sam1, sam2)
  # Add jittering (needed?) and convert to km
  D = log((Draw + runif(n)*50) / 1000)
  
  # Compute weighted means (of logs)
  md <- mean(D)
  p1l <- log(p1)
  p2l <- log(p2)
  
  if(gravity){
    Pi <- mean(p1l)
    Pj <- mean(p2l)
    sigma2 = cov(cbind(D, p1l, p2l))
    return(list(mu = c(md, Pi, Pj), sigma2 = sigma2))
  } else {
    sigma2 = var(D)
    return(list(mu = md, sigma2 = sigma2))
  }
}
