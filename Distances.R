# Function to compute distance measures between points
library("geosphere")
library("osrm")
library("Rcpp")
library("raster")

sourceCpp("Radiation.cpp")

beeline = function(origin,
                   destination,
                   pairwise = F,
                   log = F){
  n = c(nrow(origin), nrow(destination))
  nMin = min(n)
  if(!pairwise){
    out = distGeo(p1 = origin[1:nMin, 1:2], p2 = destination[1:nMin, 1:2])/1000
    if(log) out = log(out)
    return(out)
  }else{
    out = matrix(0, nrow = n[1], ncol = n[2])
    for(i in 1:n[1]){
      out[i, ] = distGeo(origin[i, 1:2], destination[, 1:2])/1000
    }
    if(log) out = log(out)
    return(c(out))
  }
}

# Returns values and counts
travelTime = function(origin,
                      destination,
                      pairwise = F,
                      batchsize = 100,
                      log = F){
  n = c(nrow(origin), nrow(destination))
  nMin = min(n)
  
  if(!pairwise){
    out = rep(0, nMin)
    for(i in 1:ceiling(nMin/batchsize)){
      ind = ((i-1)*batchsize):min((i*batchsize), nMin)
      out[ind] = diag(osrmTable(src = sam1[ind, 1:2], dst = sam2[ind, 1:2])$durations)
    }
    if(log) out = log(out)
    outTab = table(out)
    outNum = as.numeric(names(outTab)) + 0.1
    if(log) outNum = log(outNum)
    return(list(time = outNum,
                count = as.numeric(outTab)))
  }else{
    out = matrix(0, nrow = n[1], ncol = n[2])
    for(i in 1:ceiling(n[1]/batchsize)){
      ind1 = ((i-1)*batchsize):min((i*batchsize), n[1])
      for(j in 1:ceiling(n[2]/batchsize)){
        ind2 = ((j-1)*batchsize):min((j*batchsize), n[2])
        out[ind1, ind2] = osrmTable(src = sam1[ind1, 1:2], dst = sam2[ind2, 1:2])$durations
      }
    }
    out = c(out)
    outTab = table(out)
    outNum = as.numeric(names(outTab)) + 0.1
    if(log) outNum = log(outNum)
    return(list(time = outNum,
                count = as.numeric(outTab)))
  }
}


gravity = function(origin,
                   destination,
                   pairwise = F,
                   log = F){
  n = c(nrow(origin), nrow(destination))
  nMin = min(n)
  if(!pairwise){
    out = distGeo(p1 = origin[1:nMin, 1:2], p2 = destination[1:nMin, 1:2])/1000
    out = cbind(out, sam1[, 3]/sam2[, 3])
    if(log) out = log(out)
    return(out)
  }else{
    out = matrix(0, nrow = n[1], ncol = n[2])
    for(i in 1:n[1]){
      out[i, ] = distGeo(origin[i, 1:2], destination[, 1:2])/1000
    }
    out = cbind(c(out), c(sam1[, 3]*t(sam2[, 3])))
    if(log) out = log(out)
    return(out)
  }
}


# Helper function for Circle and Radiation measure
replaceIndices = function(coords, rangeLat, rangeLong, dims){
  i = dims[1] - floor((coords[, 2] - rangeLong[1])/(rangeLong[2] - rangeLong[1])*dims[1])
  j = ceiling((coords[, 1] - rangeLat[1])/(rangeLat[2] - rangeLat[1])*dims[2])
  return(cbind(i, j))
}

# Additional input:
# popMat: Matrix from raster of population density
#  -Should have no NAs
# extent: Limits of population raster
#  - for example extent(pop) 

circle = function(origin,
                  destination,
                  popMat,
                  extent,
                  pairwise = F,
                  log = F){
  n = c(nrow(origin), nrow(destination))
  nMin = min(n)
  
  # Extend samples by indices of population matrix
  Iorigin = replaceIndices(origin[, 1:2],
                      rangeLong = c(ymin(extent), ymax(extent)),
                      rangeLat = c(xmin(extent), xmax(extent)),
                      dims = dim(popMat))
  Idestination = replaceIndices(destination[, 1:2],
                                rangeLong = c(ymin(extent), ymax(extent)),
                                rangeLat = c(xmin(extent), xmax(extent)),
                                dims = dim(popMat))
  
  if(!pairwise){
    out = rep(0, nMin)
    for(i in 1:nMin){
      ind1 = Iorigin[i, ]
      ind2 = Idestination[i, ]
      r = ceiling(sqrt(sum((ind1 - ind2)^2)))
      out[i] = computePop(ind1[1], ind1[2], r, popMat)
    }
    out = out/100000
    if(log) out = log(out)
    return(out)
  }else{
    out = matrix(0, nrow = n[1], ncol = n[2])
    for(i in 1:n[1]){
      for(j in 1:n[2]){
        ind1 = Iorigin[i, ]
        ind2 = Idestination[j, ]
        r = ceiling(sqrt(sum((ind1 - ind2)^2)))
        out[i, j] = computePop(ind1[1], ind1[2], r, popMat)
      }
    }
    out = c(out)/100000
    if(log) out = log(out)
    return(out)
  }
}


radiation = function(origin,
                     destination,
                     popMat,
                     extent,
                     pairwise = F,
                     log = F){
  n = c(nrow(origin), nrow(destination))
  nMin = min(n)
  
  # Extend samples by indices of population matrix
  Iorigin = replaceIndices(origin[, 1:2],
                           rangeLong = c(ymin(extent), ymax(extent)),
                           rangeLat = c(xmin(extent), xmax(extent)),
                           dims = dim(popMat))
  Idestination = replaceIndices(destination[, 1:2],
                                rangeLong = c(ymin(extent), ymax(extent)),
                                rangeLat = c(xmin(extent), xmax(extent)),
                                dims = dim(popMat))
  Ni = origin[1:nMin, 3]
  Nj = destination[1:nMin, 3]
  
  if(!pairwise){
    S = rep(0, nMin)
    for(i in 1:nMin){
      ind1 = Iorigin[i, ]
      ind2 = Idestination[i, ]
      r = ceiling(sqrt(sum((ind1 - ind2)^2)))
      S[i] = computePop(ind1[1], ind1[2], r, popMat)
    }
    out = ((Ni + S)*(Ni + Nj + S))/(Ni*Nj)
    out = out/1e09
    if(log) out = log(out)
    return(out)
  }else{
    S = matrix(0, nrow = n[1], ncol = n[2])
    for(i in 1:n[1]){
      for(j in 1:n[2]){
        ind1 = Iorigin[i, ]
        ind2 = Idestination[j, ]
        r = ceiling(sqrt(sum((ind1 - ind2)^2)))
        S[i, j] = computePop(ind1[1], ind1[2], r, popMat)
      }
    }
    Ni = matrix(rep(Ni, length(Nj)), nrow = length(Ni))
    Nj = matrix(rep(Nj, nrow(Ni)), ncol = length(Nj), byrow = T)
    out = ((Ni + S)*(Ni + Nj + S))/(Ni*Nj)
    out = c(out)/1e09
    if(log) out = log(out)
    return(out)
  }
}