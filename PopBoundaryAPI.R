library("raster")
library("sf")
library("geosphere")

############################################
### Goal: Prepare data for FGSIM model   ###
### We need: District boundary data      ###
###          Population density data     ###
############################################


# Population counts downloaded from worldpop

# Function objective: Download population data as raster object
# Checks if the data is already in the directory to avoid double downloads
# Disclaimer: Takes a bit to download, timeout should be set higher

# Input:
# countryCode:  ISO-Code for country, see
#               https://en.wikipedia.org/wiki/List_of_ISO_3166_country_codes
# year:         The year for which data is downloaded
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