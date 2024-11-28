### Function to sample coordinates from a district
### Samples a total of n observations, n at a time uniformly in the district
### then uses population density to either accept a sample or reject
### until a total of n observations are accepted.

# poly:     SpatialPolygon of district
# n:        Sample size
# pop:      Population raster
# weighted: Use population as weights?

wspsample = function(poly, n, pop, weighted = F){
  # population data cropped to the MBR of the district
  popi    = mask(crop(pop, poly), poly)
  # Obtain all population values in the district's MBR
  popAll  = values(popi)
  # Remove NAs
  indKeep = which(!is.na(popAll))
  # Population without any NA or zero
  popAll = popAll[indKeep]
  # Respective coordinates
  co = coordinates(popi)[indKeep, ]
  
  if(weighted){
    w = popAll/sum(popAll)
  }else{
    w = as.numeric(popAll > 0)
  }
  
  # Sample indices
  ind = sample(x = 1:length(popAll), size = n, replace = T, prob = w)
  out = cbind(co[ind, ], 0)
  
  # Jittering to avoid duplicates in sample 
  delta1 = min(diff(sort(unique(co[, 1]))))
  delta2 = min(diff(sort(unique(co[, 2]))))
  out[, 1] = out[, 1] + runif(n)*delta1 - delta1/2
  out[, 2] = out[, 2] + runif(n)*delta2 - delta2/2
  out[, 3] = popAll[ind]
  
  return(out)
}

# Notes:

# For weighted sampling according to population density, we could:
# 1) Sample uniformly and use population as weights afterwards
# 2) Sample a cell according to population inside MBR and then check if the
#    cell is inside the district
# 3) Sample a cell uniformly inside district and reject the cell with probability
#    p_i/p_max, where p_i is the population in that cell and p_max is the maximum
#    population in all cells of the district

# 1) could be misleading for small sample sizes
# 2) is slow because of repeated checking and potentially huge MBRs
# 3) is a good trade-off and implemented



