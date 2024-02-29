### Function to sample coordinated from a district
### Samples a total of n observations, n at a time uniformly in the district
### then uses population density to either accept a sample or reject
### until a total of n observations are accepted.

# poly:     SpatialPolygon of district
# n:        Sample size
# pop:      Population raster
# weighted: Use population as weights?
wspsample = function(poly, n, pop, weighted = F){
  # population data cropped to the MBR of the district
  popi    = crop(pop, poly)
  # Obtain all population values in the district's MBR
  popAll  = values(popi)
  # Remove NAs
  popAll[is.na(popAll)] = 0
  # Find maximum population
  M       = max(popAll)
  # Sample uniformly
  sam1_df = spsample(poly, n, "random")@coords
  # Get population for samples locations
  p = extract(popi, sam1_df)
  # Calculate acceptance probabilities
  # Excludes zero-population locations also for unweighted approach
  if(weighted){
    acc = p / M
  }else{
    acc = p > 0
  }
  # Check which locations are accepted
  indKeep = which(acc >= runif(n))
  # Initiate output data
  out = sam1_df[indKeep, ]
  
  # Repeat above steps until we have enough points
  while(nrow(out) < n){
    sam1_df = spsample(poly, n, "random")@coords
    p = extract(popi, sam1_df)
    p[is.na(p)] = 0
    if(weighted){
      acc = p / M
    }else{
      acc = p > 0
    }
    indKeep = which(acc >= runif(n))
    out = rbind(out, sam1_df[indKeep, ])
  }
  out[1:n, ]
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
# 3) is a good trade-off and implemented below