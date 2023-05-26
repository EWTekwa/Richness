# Richness
R-package "Richness" for asymptotic species richness estimates
Eden Tekwa and Matt Whalen

Summary:
The R-Package takes spatial abundance data and return estimated species richness, including observed richness, Chao1, Chao2, ACE, Jackknife-abundance, Jackknife-incidence, 𝛺o (base-version of estimator developed in citation), and 𝛺T (recommended estimator developed in citation).

A wrapper function estimateRichness.R allows the user to provide community data as a matrix, data.frame, list, or array. Lists or arrays can be input to compare multiple communities, such as a time series.

The wrapper function calls two functions RichnessEsts.R and bootstrapRichness.R that respectively return point estimates and bootstrapped estimates for estimator precision/uncertainty quantification. These functions take in data matrix with species as columns, spatial units (e.g., transects or quadrats) as rows, and individual counts as values.


See "R package Richness Documentation.pdf" for installation and operation in R.
See /Matlab subfolder for Matlab code.

Citation:	Tekwa, Whalen, Martone, O’Connor, Theory and application of an improved species richness estimator 2023. Philosophical Transactions B. 10.1098/rstb.2022.0181.2022-0187.

# Installation
To download package in R, type:
```
library(devtools)
install_github("EWTekwa/Richness")
library(Richness)
```

To open help files
for functions
in R, type:
```
?RichnessEsts
?estimate_richness
?bootstrapRichness
```
