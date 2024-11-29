# FGSIM - Fine Grid Spatial Interaction Matrix

This repo is a toolbox that serves as an addon to the [surveillance](https://cran.r-project.org/web/packages/surveillance/index.html) package and its hhh4 function for spatio-temporal modeling of disease spread.

In spatio-temporal models, we need to specify the transmission between districts in a weight matrix. Common approaches are based on the adjacency of districts / the neighbourhood structure, which depends on where district boundaries are. The methods provided here are based on different distance measures serve two purposes:

1) Transmission risk is modeled on the individual level instead of the district level to reduce the influence of boundary locations and to unify the definition of transmission risk for intra- and inter-district transmission

2) Detailed population density data is included in hope to improve forecast performance and to enable risk mapping on resolution finer than the observed district level

## Theory

When modeling the risk of transmitting disease in a spatio-temporal setting, we typically express the risk as a function of a distance measure. For example if we use the order of neighbourhood as distance $D$, it is common to assume a power law relationship between distance and transmission risk $D \mapsto D^{-d}$. A similar idea is used by the approach provided here. Instead of defining a distance measure on a district level, we define it on the individual level, i.e. if we select two individuals at random, one infected and one susceptible, we can for example compute the beeline distance $B$ between them. Is it reasonable to assume that the transmission risk decays with increasing distance between them, we may also assume a power law relationship, i.e. $B^{-d}$. Due to the random selection of an infected and a susceptible person, the distance measure between them is a random variable. Given that we selected individuals from certain districts, we can derive the distribution of the distance between individuals from population density data. Then, we can aggregate the random risk measure $B^{-d}$ to the district level by taking the expectation and use it to construct a weight matrix. As a distribution to fit, the truncated LogNormal distribution is a good choice, since it accounts for positivity and boundedness of distance measures and yields a closed form expectation of $B^{-d}$. Fitting a mixture of truncated LogNormals increases flexibility and preserves those properties.

This toolbox is designed for applications with spatio-temporal data, i.e. time series of infection counts for different districts. With such a data set the following steps are needed:

### 1) Downloading district boundaries
#### $\rightarrow$ PopBoundaryAPI.R

The function `getJSON` downloads shapefiles with district boudnaries.

### 2) Downloading population density data
#### $\rightarrow$ PopBoundaryAPI.R

The function `getPOP` downloads population density data from [Worldpop](https://worldpop.org) and saves it in a raster file.


### 3) Match boundary data and population density data to case data districts 

Needs to be done by hand.

### 4) Sampling individuals in each district according to population density #### $\rightarrow$ wspsample.R

The function `wspsample` takes a polygon from the shapefile (step 1), a sample size and population density data (step 2) as input and returns a sample of locations where individuals live. The output is a matrix where each row contains coordinates and the population estimate in the corresponding 100m x 100m grid cell. The boolean argument `weighted` specifies if the sampling should be weighted by population. 

### 5) Compute distance measure between individuals
#### $\rightarrow$ Distances.R

Different distance measures are implemented that each take two matrices of coordinates as input, one from an origin district and one from the destination district. The argument `pairwise` specifies if distances shall be computed pairwise, and the argument `log` if the distances shall be returned in logs.

Currently implemented are:
- `beeline`: Beeline distance
- `travelTime`: The time to travel from one point to another (by car, in minutes)
- `gravity`: Like the beeline distance but in addition the ratio of population densities at origin and destination is returned as well
- `circle`: Needs additional input. From the population (step 2), create a matrix (`as.matrix(pop)`) and an extent object (`extent(pop)`). It then computes the population in a circle around the origin that touches the destination. It makes use of the midpoint circle algorithm and C++ for speed-up.
- `radiation`: Similar to the `circle` measure $C$, but also includes population dennsity at origin $P_O$ and destination $P_D$. It is the inverse of the radiation formula for human migration such that larger $R$ indicates larger distance.
$$
R = \frac{(P_o + C)\cdot(P_O + P_D + C)}{P_O\cdot P_D}
$$

### 6) Fitting a distribution to each sample of distances
#### $\rightarrow$ Fitting.R

Using the sample of distances (step 5), a mixture of truncated LogNormal distributions is fitted. For the gravity model, a function fits a mixture of untruncated bivariate Normal distributions. It uses an EM algorithm for estimation and returns the estimates as well as the BIC to determine the number of mixture components.

### 7) Saving parameter estimates

To use the fitted distributions in a weight function for `hhh4`, they must be saved in a list with elements
- `mu`: Mean parameters
- `sigma`: Standard deviation parameters
- `w`: Weights
- `l`: Lower bound of distance measure (in logs)
- `u`: Upper bound of distance measure (in logs)
- `pop`: diagonal matrix of district population

The first three elements should be of dimension $N\times N \times k$, where $N$ is the number of districts and $k$ the number of mixture components. The bounds are assumed to be equal for all mixture compoments, thus of dimension $N\times N$, and the population matrix is also of dimension $N\times N$.

For the gravity model, the distance measure is two dimensional, thus we have list elements `mu1` and `mu2` instead of `mu` for mean parameters and 
`sigma1`, `sigma2` and `sigma12` instead of `sigma` for the elements of the variance covariance matrix.


### 8) Defining how weights are computed in hhh4
#### $\rightarrow$ W_fgsim

Finally, we define weights for `hhh4` with the function `W_fgsim`. It is similar to the `W_powerlaw` function in `surveillance` and returns a list of functions that are used to compute weights and their derivatives as well as initial values for the decay parameter. Its arguments are

- `pars`: List of parameters (step 7)
- `truncPL`: Shall the power law be truncated: Instead of translating distances $D$ to $D^{-d}$, it would translate to $\min\{D^{-d}, \delta^{-d}\}$, where $\delta$ is a threshold parameter that is estimated. 
- `gravity`: Set to TRUE if the distance measure is `gravity`
- `maxlag`: The maximum neighbourhood order with non-zero weights
- `normalize`: Shall rows of the weight matrix be normalised?
- `log`: Are parameters estimated in logs to ensure positivity?
- `initial`: Initial values for estimation
- `from0`: If `TRUE`, weights are also computed for intra-district transmission
- `popScale`: If `TRUE`, the weights are multiplied by the population size of the destination district to account for heterogeneity


### References

- [Barbosa H. et al. (2018) **Human mobility: Models and applications**, Physics Reports, 734:1–74.](https://doi.org/10.1016/j.physrep.2018.01.001)
- [Giraud, T. (2022) **osrm: Interface Between R and the OpenStreetMap-Based Routing Service OSRM**, Journal of Open Source Software, 7(78), 4574.](https://cran.r-project.org/web/packages/osrm/osrm.pdf)
- [Held L., Höhle M., Hofmann M. (2005) **A statistical framework for the analysis of multivariate infectious disease surveillance counts**, Statistical Modelling, 5(3):187–199.](https://doi.org/10.1191/1471082X05st098oa)
- [Meyer S., Held L. (2014) **Power-law models for infectious disease spread**, The Annals of Applied Statistics, 8(3), 1612–1639.](https:/doi.org/10.1214/14-AOAS743)
- [Meyer S., Held L., Höhle M. (2017) **Spatio-Temporal Analysis of Epidemic Phenomena Using the R Package surveillance**, Journal of Statistical Software, 77(11), 1–55.](https://cran.r-project.org/web/packages/surveillance/index.html)
- [Runfola D. et al. (2020) **geoboundaries: A global database of political administrative boundaries**, PLoS ONE, 15(4).](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0231866)
- [WorldPop (2018) **Global high resolution population denominators project**, opp1134076.](https://worldpop.org)