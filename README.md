# FGSIM - Fine Grid Spatial Interaction Matrix

The objective of this repo is to provide an extension to the hhh4 function of the surveillance package. It allows to use newly defined weight matrices that shall capture the closeness of administrative districts.

## Idea
To define the closeness of two districts $i$ and $j$, we define $D_{ij}$ as the distance between two randomly selected individuals, one from district $i$ and one from $j$. Due to the random selection of individuals, $D_{ij}$ is a random variable. The distribution of $D_{ij}$ depends on the selection procedure. We can either select individuals while taking into account how a district's population is distributed across the district or we can select a location inside the district assuming uniform distribution of individuals across the district. Both selection procedures yield a distribution of $D_{ij}$ that is well approximated by a LogNormal distribution. Let $\mu_{ij}$ and $\sigma^2_{ij}$ be the parameters of the LogNormal distribution. Following the idea of powerlaw models, we define raw weights as $$W_{ij} = \text{E}\left(D_{ij}^{-d}\right) = \exp\left(-d \mu_{ij} + \frac{d^2\sigma_{ij}^2}{2}\right)$$ where $d$ is a decay parameter. The parameters $\mu_{ij}$ and $\sigma_{ij}^2$ are estimated by simulation. By the definition of raw weights, we can express the closeness of a district to itself, i.e. $W_{ii}$ is well-defined.

Using a similar idea, gravity model weights can be derived. Let $P_i$ and $P_j$ be the population densities at locations of randomly selected individuals. It is assumed that not only the distance follows a LogNormal distribution, but that the triple $(D_{ij}, P_i, P_j)$ follows a multivariate LogNormal, i.e. the log of the three components follows a multivariate Normal. Let $\mu_{ij}$ and $\Sigma_{ij}$ be the parameters of that distribution. The gravity, or attraction with with an individual of district $i$ is pulled to district $j$, is defined as $$D_{ij}^{-d_1}P_i^{-d_2}P_j^{d_3} .$$
If all decay parameters are positive, the force of attraction decreases with increasing distance between districts and with increasing population density in the district of origin. Districts with larger population densities have a larger attraction on other individuals. Since $D_{ij}$, $P_i$ and $P_j$ are assumed to be random variables, we compute raw weights as
$$W_{ij} = \text{E}\left( D_{ij}^{-d_1}P_i^{-d_2}P_j^{d_3} \right) = \exp\left(d'\mu_{ij} + \frac{d'\Sigma_{ij}d}{2}\right)$$
where $d = (-d_1, -d_2, d_3)'$.

## Code Overview
The repo provides functions for data preparation, that is carried out before the application and functions to compute weight matrices.

`PopBoundaryAPI.R` contains two function to load data on district boundaries and population density on a fine grid. Population density data is downloaded from [Worldpop](https://www.worldpop.org/), district boundaries are downloaded from [geoBoundaries](https://github.com/wmgeolab/geoBoundaries).

`wspsample.R` contains a function to sample points inside a selected district. Users can select to sample either individuals according to population density data (`weighted = TRUE`) or locations uniformly across the district (`weighted = FALSE`).

`fitLN.R` takes samples of two districts (from `wspsample`) and computes estimates of the LogNormal parameters. If both samples are weighted, the resulting distribution gives the distance between individuals from districts. If the first sample is weighted and the second sample is uniform, the resulting distribution gives the distance between an individual from a district of origin $i$ and a location in the destination district $j$.

`RawWeights.R` contains three functions to compute raw weights and the first two derivatives w.r.t the decay parameter(s). Either the approach based only on distances, the gravity model or a powerlaw applied to neighbourhood order.

`DerivativeHelper.R` provides functions that carry out transformations of weights and their derivatives. Four transformations are considered:
1) Row normalisation: $W_{ij} \rightarrow \frac{W_{ij}}{\sum_{k = 1}^n W_{ik}}$
2) Log-transformation of parameters: $d \rightarrow exp(d)$
3) Scaling the columns by a diagonal matrix ($P$): $W \rightarrow W \cdot P$
4) Translating weights to contact probabilities: $W \rightarrow W \cdot M \cdot W^{T}$

`FGSIMdistance.R` provides the function that computes weights not including gravity effects. Similar to the `W_powerlaw` function of the `surveillance` package, it can be used as

`W_pdist(pars, maxlag, normalize, log, initial, from0, areaScale, popScale, contactScale)`

Setting `maxlag` to a finite integer sets weights to zero if two districts have a neighbouring order larger that the maximum lag selected. `pars` needs to be a list with at least two entries `dist` and `s11` that contain the mean and variance parameter of the LogNormal distribution. Each needs to be an $n\times n$ matrix where $n$ is the number of districts. `normalize` is set to `TRUE` if the rows of the matrix shall be normalized. If we want to estimate $d$ with the restriction that is positive, we set `log` to `TRUE`. `initial` shall provide a starting value of $d$ for estimation and `from0` can be set to `FALSE` if we want to set diagonal elements of the weight matrix to zero. The three transformations `areaScale`, `popScale` and `contactScale` are motivated as follows:

`areaScale`: If we treat the raw weights as attraction of a location in $j$ on an individual from $i$, we can account for different district areas, i.e. number of locations, by scaling the column of weights by the destination area. In that case, the argument `pars` needs to contain a daigonal matrix `area`. For numerical stability, is is reasonable to scale the matrix to have unit determinant. Normalising the rows afterwards, the resulting weights then approximate the probability that an individual from $i$ moves to district $j$.

`contactScale`: If we interpret initial weights $W_{ij}$ as probabilities that an individual from $i$ is in district $j$, we may want to account for two individuals having contact in a third district. The probability that two individuals are in the same district is computed as $W\dot W^{T}$, iterating over all districts the two individuals might meet in. Since districts may exhibit very different characteristics in area and population density, we can include the probability that two individuals have contact given that they are in the same district. If we define a diagonal matrix $M$ such that $M_{ii}$ gives the probability of two individuals having contact, given that they are in district $i$, we define contact probabilities as $W \dot M \dot W^{T}$. There are different ways of defining elements $M_{ii}$: Either we can make use of the distribution of $D_{ii}$ and compute $M_{ii} = P(D_{ii}  \le \epsilon)$ for some small value of $\epsilon$. On the other hand we can define "having contact" as being in the same 100m by 100m gridcell. If we set `contactScale` to `TRUE`, the list `pars` must contain a matrix `M`, which for stability shall be scaled to have unit determinant.

`popScale`: Starting with weights $W_{ij}$ that are interpreted as the probability that an individual from $i$ has contact with an individual from $j$, we can account for different population sizes in districts $j$ by applying the population scaling. That means, we multiply the columns of weights by the corresponding district's population size to obtain the expected number of contacts an individual from $i$ has had with individuals from $j$. If we then normalise rows of the weight matrix, $W_{ij}$ is interpreted as the probability that an individual from $i$ has had contact with an individual from district $j$ given that the individual from $i$ has had contact with one individual. 

`FGSIMgravity.R` provides a weight function that does include gravity effects. Similar to above case, we use is as

`W_gravity(pars, maxlag, normalize, log, initial, from0, areaScale, popScale, contactScale)`

In contrast to the `W_pdist` function, we provide the LogNormal parameters as a list containing entries: `dist, pO, pD, s11, s12, s13, s22, s23, s33`. The former three are the mean parameters for the distance, population in the district of origin and district of destination respectively. `s11`, ... are the elements of the covariance matrix. Each of the nine entries again needs to be an $n \times n$ matrix. The initial values `initial` now need to be a vector of length 3.

`PowerLaw.R`: This file contains a wrapper function that is similar to `W_powerlaw` from the `surveillance` package but allows to apply above three transformation steps `areaScale`, `popScale` and `contactScale`.

