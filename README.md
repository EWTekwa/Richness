# Richness
R-package "Richness" for asymptotic species richness estimates
Eden Tekwa and Matt Whalen

Version: 1.1
Depends: R 4.2.2
Published: 2023-04-27
Last update: 2023-05-26
Authors: Matt Whalen and Eden Tekwa
Maintainers: Matt Whalen (mawhal@gmail.com)
Eden Tekwa (ewtekwa@gmail.com)
License:
CC BY-NC-SA
URL: https://github.com/EWTekwa/Richness
Citation: 
Tekwa, Whalen, Martone, O’Connor, Theory and application of an improved species richness estimator 2023. Philosophical Transactions of the Royal Society B. https://doi.org/10.1098/rstb.2022.0181.2022-0187

Summary:
The R-Package uses spatial abundance data to produce estimated species richness and estimator precision while correcting for observation biases from low occupancies and abundances. The featured estimators include a set from Tekwa et al 2023: 𝛺 (exact estimator, recommended version), 𝛺T (second-order Taylor approximated estimator, useful for quantifying biases from means and variances in occupancy and observed abundance across species), and 𝛺o (zeroth-order approximated estimator). Additionally, standard estimates are produced using Chao1, Chao2, ACE, Jackknife-abundance, Jackknife-incidence, and observed (uncorrected) richness.

A wrapper function estimateRichness.R allows the user to provide community data as a matrix, data.frame, list, or array. Lists or arrays can be input to compare multiple communities, such as a time series.

The wrapper function calls two functions RichnessEsts.R and bootstrapRichness.R that respectively return point estimates and bootstrapped estimates for estimator precision quantification. These functions take in data matrix with species as columns, spatial units (e.g., transects or quadrats) as rows, and individual counts as values. Only bootstrapped estimates are provided, and it is up to the user to choose how to use these for estimating precision and possibly confidence bounds. For example, precision can be quantified by taking the standard deviation of the bootstrapped estimates.

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
?estimateRichness
?RichnessEsts
?bootstrapRichness
```
