\name{RichnessEsts}
\alias{RichnessEsts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Richness Estimates with Correction Terms and Detection Probabilities
%%  ~~function to do ... ~~
}
\description{
This function takes spatial abundance data and returns several richness estimates, including three defined by Tekwa et al. 2023.  In addition to richness point estimates, the function returns mean states of correction terms and approximated detection probabilities used in Tekwa's estimation method.
}
\usage{
RichnessEsts(Community)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Community}{ Data frame or matrix of community data with species as columns and spatial sampling units as rows. }
}
\details{
% can include plenty of details about the method introduced and the other richness estimators.
}
\value{
A list of 3 elements.
\item{element1 }{point estimates for all richness estimators. Estimators include raw richness, Chao1, Gamma-Poisson, Chao2, ACE, Jackknife-abundance, Jackknife-incidence, \eqn{\Omega}, \eqn{\Omega{T}}, \eqn{\Omega{T0}}.}
\item{element2 }{mean states of the correction terms for \eqn{\Omega}, \eqn{\Omega{T}}, and \eqn{\Omega{T0}}.}
\item{element3 }{approximated detection probabilities for \eqn{\Omega}, \eqn{\Omega{T}}, and \eqn{\Omega{T0}}.}
}
\references{
Tekwa, E.W., Whalen, M.A., Martone, P.T. and O'Connor, M.I. (2023) Theory and application of an improved species richness estimator. \emph{Philosophical Transactions of the Royal Society B.} \doi{10.1098/rstb.2022.0181.2022-0187}

References for Chao1, Gamma-Poisson, Chao2, ACE, and Jackknife richness estimators:

Chao A. 1984. Nonparametric estimation of the number of classes in a population. \emph{Scandinavian J. Stat.} 11, 265–270. https://www.jstor.org/stable/4615964

Chiu CH. 2023. A more reliable species richness estimator based on the
Gamma–Poisson model. \emph{PeerJ.} 11, e14540. \doi{10.7717/peerj.14540}

Chao A. 1987 Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics.} 43, 783–791. \doi{10.2307/2531532}

Chao A, Lee SM. 1992 Estimating the number of classes via sample coverage.
\emph{J. Am. Stat. Assoc.} 87, 210–217. \doi{10.1080/01621459.1992.10475194}

Burnham KP, Overton WS. 1979. Robust estimation of population size when capture probabilities vary among animals. \emph{Ecology.} 60, 927–936. \doi{10.2307/1936861}
}
\author{
E.W. Tekwa and M.A. Whalen
}
\note{
% other notes about the package. Appears between Value and Authors
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{estimateRichness}}, \code{\link{bootRichnessEsts}}.
}
\examples{
## install the package
require(devtools)
install_github("EWTekwa/Richness")
library(Richness)

## Use Barro Colorado Island data from the vegan package
require(vegan)
data("BCI")

## run the function, which returns a list
RichnessEsts(BCI)

## isolate the richness estimates
RichnessEsts(BCI)[[1]]
}
