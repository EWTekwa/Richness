# Richness
#### R package `Richness` for asymptotic species richness estimates (see https://github.com/EWTekwa/Richness/blob/main/R%20package%20Richness%20Documentation.pdf for full documentation)

This R Package uses spatial abundance data to produce estimated species richness and estimator precision while correcting for observation biases from low occupancies and abundances. The featured estimators include a set from Tekwa et al 2023: 𝛺 (exact estimator, recommended version), 𝛺ᴛ (second-order Taylor approximated estimator, useful for quantifying biases from means and variances in occupancy and observed abundance across species), and 𝛺o (zeroth-order approximated estimator). Additionally, standard estimates are produced using Chao1, Chao2, ACE, Jackknife-abundance, Jackknife-incidence, and observed (uncorrected) richness.

A wrapper function `estimateRichness.R` allows the user to provide community data as a matrix, data.frame, list, or array. Lists or arrays can be input to compare multiple communities, such as a time series.

The wrapper function calls two functions `RichnessEsts.R` and `bootstrapRichness.R` that respectively return point estimates and bootstrapped estimates for estimator precision quantification. These functions take in data matrix with species as columns, spatial units (e.g., transects or quadrats) as rows, and individual counts as values. Only bootstrapped estimates are provided, and it is up to the user to choose how to use these for estimating precision and possibly confidence bounds. For example, precision can be quantified by taking the standard deviation of the bootstrapped estimates.


See `/Matlab` subfolder for Matlab code.


## Installation
To download package in R, type:
```r
library(devtools)
install_github("EWTekwa/Richness")
library(Richness)
```

To open help files for functions in R, type:
```r
?estimateRichness
?RichnessEsts
?bootstrapRichness
```

## References

Tekwa, E.W., M.A. Whalen, P.T. Martone, & M.I. O’Connor. 2023. Theory and application of an improved species richness estimator. _Philosophical Transactions of the Royal Society B_. https://doi.org/10.1098/rstb.2022.0181.2022-0187
