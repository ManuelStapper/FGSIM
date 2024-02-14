source("getLNpars.R")

# Example for NUTS-3 districts in Germany 2011
countryCode = "DEU"
level = 3
year = 2011

# Download data
distr = getJSON(countryCode, level)
# Translate to spatial object
distrS = as_Spatial(distr$geometry)

# Download population data
pop = getPOP(countryCode, year, 1000)

# To actually work with the data, we need to inspect which polygons
# are water. The table IDlookup contains the RKI/NUTS, names of districts,
# IDs in the JSON object and IDs in JSON that are water.

# Disclaimer:
# A similar table needs to be created for other countries depending
# on other data sets. Here we match the districts to RKI infection data

IDlookup = read.csv("IDlookup.csv", header = T, fileEncoding = "UTF-8")
IDs = IDlookup$ID


### Unweighted

sam = list()
for(i in 1:401){
  set.seed(i)
  sam[[i]] = wspsample(distrS[IDlookup$JsonNr[i]], 1000, pop, weighted = F)
}

# Now fit LN and save parameters

gerPars = data.frame(iO   = rep(0, 80601),
                     iD   = rep(0, 80601),
                     dist = rep(0, 80601),
                     pO   = rep(0, 80601),
                     pD   = rep(0, 80601),
                     s11  = rep(0, 80601),
                     s12  = rep(0, 80601),
                     s13  = rep(0, 80601),
                     s22  = rep(0, 80601),
                     s23  = rep(0, 80601),
                     s33  = rep(0, 80601))

counter = 1
for(i in 1:401){
  for(j in i:401){
    temp = fitLN(sam[[i]], sam[[j]], IDlookup$JsonNr[i], IDlookup$JsonNr[j], distrS, pop, gravity = T)
    gerPars$iO[counter]   = i
    gerPars$iD[counter]   = j
    gerPars$dist[counter] = temp$mu[1]
    gerPars$pO[counter]   = temp$mu[2]
    gerPars$pD[counter]   = temp$mu[3]
    gerPars$s11[counter]  = temp$sigma2[1, 1]
    gerPars$s12[counter]  = temp$sigma2[1, 2]
    gerPars$s13[counter]  = temp$sigma2[1, 3]
    gerPars$s22[counter]  = temp$sigma2[2, 2]
    gerPars$s23[counter]  = temp$sigma2[2, 3]
    gerPars$s33[counter]  = temp$sigma2[3, 3]
    counter = counter + 1
    print(c(i, j))
  }
}

# write.csv(gerPars, "gerPars.csv")

### Weighted
samW = list()
for(i in 1:401){
  set.seed(i)
  samW[[i]] = wspsample(distrS[IDlookup$JsonNr[i]], 1000, pop, weighted = T)
}

gerParsW = gerPars

counter = 1
for(i in 1:401){
  for(j in i:401){
    temp = fitLN(samW[[i]], samW[[j]], IDlookup$JsonNr[i], IDlookup$JsonNr[j],
                 distrS, pop, gravity = T)
    gerParsW$iO[counter]   = i
    gerParsW$iD[counter]   = j
    gerParsW$dist[counter] = temp$mu[1]
    gerParsW$pO[counter]   = temp$mu[2]
    gerParsW$pD[counter]   = temp$mu[3]
    gerParsW$s11[counter]  = temp$sigma2[1, 1]
    gerParsW$s12[counter]  = temp$sigma2[1, 2]
    gerParsW$s13[counter]  = temp$sigma2[1, 3]
    gerParsW$s22[counter]  = temp$sigma2[2, 2]
    gerParsW$s23[counter]  = temp$sigma2[2, 3]
    gerParsW$s33[counter]  = temp$sigma2[3, 3]
    counter = counter + 1
    print(c(i, j))
  }
}

