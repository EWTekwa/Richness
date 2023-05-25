\name{bootstrapRichness}
\alias{bootstrapRichness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrapped Richness Estimates
%%  ~~function to do ... ~~
}
\description{
This function takes spatial abundance data and returns several richness estimates, including two defined by Tekwa et al. 2023.  The function uses a bootstrapping routine that can be used to obtain confidence intervals around richness point estimates.
}
\usage{
bootstrapRichness(Community, numBoot = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Community}{ Data frame or matrix of community data with species as columns and spatial sampling units as rows. }
  \item{numBoot}{ The number of bootstrapped samples to return. }
}

\details{
Point estimates are first calculated using \code{\link{RichnessEsts}} and used in the bootstrapping routine.

The bootstrapping routine first resamples rows (sampling units) with replacement, with the number of sampled units the same as the number of rows in the dataset.  Next, columns (e.g., species, taxa) are resampled with replacement such that the total number of columns is equal to the number of types calcuated for each richness estimate, generating a data frame for each richness estimate.  Then, each richness estimate is calcualted on each resampled data frame, and the bootstap routine repeats until  \code{numBoot} is reached.
}
\value{
A data frame of dimensions c(\code{numBoot}, 8), where the columns represent the different richness estimates and rows are bootstrapped samples.

Note that the first row contains point estimates calculated using \code{\link{RichnessEsts}}.
}
\references{
Tekwa, E.W., Whalen, M.A., Martone, P.T. and O'Connor, M.I. (2023) Theory and application of an improved species richness estimator. \emph{Philosophical Transactions of the Royal Society B.} \doi{10.1098/rstb.2022.0181.2022-0187}
}
\author{
E.W. Tekwa and M.A. Whalen
}
\note{
%%  ~~further notes~~
}
\seealso{
 \code{\link{estimateRichness}}, \code{\link{RichnessEsts}}.
}
\examples{
## install the package
require(devtools)
install_github("EWTekwa/Richness")
library(Richness)

## Use Barro Colorado Island data from the vegan package
require(vegan)
data("BCI")

## run the function
bootstrapRichness(BCI, numBoot = 10)
}
