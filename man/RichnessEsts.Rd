\name{RichnessEsts}
\alias{RichnessEsts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Richness Estimates with Correction Terms and Detection Probabilities
%%  ~~function to do ... ~~
}
\description{
This function takes spatial abundance data and returns several richness estimates, including two defined by Tekwa et al. 2023.  In addition to richness point estimates, the function returns mean states of correction terms and approximated detection probabilities used in Tekwa's estimation method.
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
\item{element1 }{point estimates for all richness estimators.}
\item{element2 }{mean states of the correction terms for \eqn{\Omega{T}} and \eqn{\Omega{0}}.}
\item{element3 }{approximated detection probabilities for \eqn{\Omega{T}} and \eqn{\Omega{0}}.}
}
\references{
Tekwa, E.W., Whalen, M.W., Martone, P.T. and O'Connor, M.I. (2023) Theory and application of an improved species richness estimator. in revision.
%% \emph{Philosophical Transactions of the Royal Society B.} \doi{}
}
\author{
E.W. Tekwa and Matt Whalen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{estimateRichness}}, \code{\link{bootRichnessEsts}}.
}
\examples{
## Use Barro Colorado Island data from the vegan package
require(vegan)
data("BCI")

## run the function, which returns a list
RichnessEsts(BCI)

## isolate the richness estimates
RichnessEsts(BCI)[[1]]
}
