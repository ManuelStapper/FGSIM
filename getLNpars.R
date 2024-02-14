library("raster")
library("sf")
library("geosphere")

############################################
### Goal: Prepare data for FGSIM model   ###
### We need: District boundary data      ###
###          Population density data     ###
############################################

# Different options:
# pdist without population density
# pdist with population density
# gravity without population density
# gravity with population density

# First part: pdist and gravity
# pdist
# We use E[D_ij^(-d)] as raw weights and assume
# log(D_ij) ~ Normal(mu_ij, sigma2_ij)
# such that D_ij^(-d) ~ LogNormal(-d*mu_ij, d^2 * sigma2_ij)
# --> E[D_ij^(-d)] = exp{-d*mu_ij + d^2 * sigma2_ij / 2}
# Distribution of D_ij stems from random selection of individuals/locations
# gravity
# We include population counts in origin (P_i) and destination (P_j)
# They can easily be saved during sampling (if available)
# Assumption: distance and the two populations are trivariate normal (in logs)
# Raw weights: E[D_ij^(-d1) * P_i^(-d2) * P_j^(d3)]
# They have a closed form by assumptions


# Second part refers to the sampling.
# Either locations are samples uniformly across the district
# or weighted by population density in cells

# Set working directory here

# Population counts downloaded from worldpop

# Function objective: Download population data as raster object
# Checks if the data is already in the directory to avoid double downloads
# Disclaimer: Takes a bit to download, timeout should be set higher

# Input:
# countryCode:  ISO-Code for country, see
#               https://en.wikipedia.org/wiki/List_of_ISO_3166_country_codes
# year:         The year for which data is downloaded. Check availability first
# timeout:      Maximum seconds to download the file
getPOP = function(countryCode, year, timeout = 100){
  options(timeout = timeout)
  filename = paste0(countryCode, year, ".tif")
  
  if(!(filename %in% list.files())){
    url = paste0("https://data.worldpop.org/GIS/Population/Global_2000_2020/", year, "/", countryCode , "/", tolower(countryCode), "_ppp_", year,".tif")
    download.file(url, destfile = paste0(countryCode, year, ".tif"), mode = "wb")
  }
  options(timeout = 100)
  return(raster(filename))
}

# Function to download a country's district boundaries as
# geojson file

# Input:
# countryCode:  ISO-3166 country code, same as above
# level:        Level administrative districts to be downloaded
#               1: Country, 2: states, 3: counties, ...
# timeout:      Maximum number of seconds for download
getJSON = function(countryCode, level = 1, timeout = 100){
  options(timeout = timeout)
  filename = paste0(countryCode, level, ".geojson")
  if(!(filename %in% list.files())){
    lvl = paste0("ADM", level)
    url = paste0("https://github.com/wmgeolab/geoBoundaries/raw/9469f09/releaseData/gbOpen/", countryCode, "/", lvl, "/geoBoundaries-", countryCode, "-", lvl,".geojson")
    download.file(url, destfile = paste0(countryCode, level, ".geojson"), mode = "wb")
  }
  options(timeout = 100)
  return(st_read(filename))
}
# An alternative would be:
# install_github("wmgeolab/rgeoboundaries")
# library(rgeoboundaries)
# ger0 = gb_adm0("germany")


# To sample points in a certain district, we have two options:
# - Ignore population density
# - Weighted sampling by population density

# For weighted sampling according to population density, we could:
# 1) Sample uniformly and use population as weights afterwards
# 2) Sample a cell according to population inside MBR and then check if the
#    cell is inside the district
# 3) Sample a cell uniformly inside district and reject the cell with probability
#    p_i/p_max, where p_i is the population in that cell and p_max is the maximum
#    population in all cells of the district

# 2) is slow because of repeated checking and potentially huge MBRs
# 1) could be misleading for small sample sizes
# 3) is a good trade-off and implemented below

# Function to sample coordinated from a district
# distri:   SpatialPolygon of district
# n:        sample size
# pop:      population raster
# weighted: Use population as weights?
wspsample = function(distri, n, pop, weighted = F){
  # population data cropped to the MBR of the district
  popi    = crop(pop, distri)
  # Obtain all population values in the district's MBR
  popAll  = values(popi)
  # Remove NAs
  popAll[is.na(popAll)] = 0
  # Find maximum population
  M       = max(popAll)
  # Sample uniformly
  sam1_df = spsample(distri, n, "random")@coords
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
    sam1_df = spsample(distri, n, "random")@coords
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


# With samples, estimate (Log)Normal parameters
# Returns parameters of the triple
# (distance, population origin, population destination)
# Or only for distance if gravity not included

fitLN = function(sam1_df,                  # Sample origin district
                 sam2_df,                  # Sample destination district
                 i,                        # District origin
                 j,                        # District destination
                 distrS,                   # Spatial object of districts
                 pop,                      # JSON file of populations
                 gravity = FALSE){         # Returns trivariate parameters
  # Sample size
  n = nrow(sam1_df)
  # Crop population data to districts
  pop1 <- crop(pop, distrS[i])
  pop2 <- crop(pop, distrS[j])
  
  # Compute population at points
  p1 = extract(pop1, sam1_df)
  p2 = extract(pop2, sam2_df)
  
  # Compute log distance (in km)
  Draw = apply(cbind(sam1_df, sam2_df), 1, function(x) distGeo(x[1:2], x[3:4]))
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