\name{estimateRichness}
\alias{estimateRichness}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate Richness in Community Data}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A wrapper function that runs \code{RichnessEsts} or  \code{bootstrapRichness} using a variety of data input types, including time series community data.
}
\usage{
estimateRichness(Community, boot = FALSE, numBoot = 100,
                             meanStates = FALSE, Apx_detectP_terms = FALSE )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Community}{ Community data with species as columns and spatial sampling units as rows. Data can take the form of a data frame, matrix, list, or array of 3 dimensions.  List and array inputs are useful when calculating richness for community time series or for different communities who richness is to be compared. }
  \item{boot}{ logical. Should bootstrapped samples be computed? }
  \item{numBoot}{ The number of bootstrapped samples to return. }
  \item{meanStates}{ logical. Should mean states of correction terms bootstrapped samples be returned? Only used when \code{boot = FALSE} }
  \item{Apx_detectP_terms}{ logical. Should approximated detection probabilities be returned? Only used when \code{boot = FALSE} }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{RichnessEsts}}, \code{\link{bootRichnessEsts}}.
}

\examples{
## install the package
require(devtools)
install_github("EWTekwa/Richness")
library(Richness)

## Use Barro Colorado Island data from the vegan package
require(vegan)
data("BCI")

## run the function for a single community data frame without bootstrapping
# point estimates only
estimateRichness(BCI)
# with mean states for correction terms and approximated detection probabilities
estimateRichness(BCI, meanStates = TRUE, Apx_detectP_terms = TRUE)

## run the function for a single community data frame with bootstrapping
estimateRichness(BCI, boot = TRUE, numBoot = 10)

## run the function on a list
# Note: this may not be the most sensible comparison
require(vegan)
data("pyrifos")
week <- gl(11, 12, labels = c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
communitylist <- split(pyrifos, f = week)
estimateRichness(communitylist)
# bootstrapping on a list returns a tidy data frame
require(tidyverse)
estimateRichness(communitylist, boot = TRUE, numBoot = 10)
}

