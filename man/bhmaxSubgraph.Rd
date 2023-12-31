\name{bhmaxSubgraph}
\alias{bhmaxSubgraph}


\title{Find maximal BH-complete subgraph}

\description{
Given an adjacency matrix of bait-hit AP-MS protein data, this function finds the maximal BH-complete subgraphs and reports them as an affiliation matrix.
}
\usage{
bhmaxSubgraph(adjMat,VBs=NULL,VPs=NULL,unrecip=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjMat}{\code{adjMat} is an N by (N+M) adjacency matrix with N equal
  to the number of bait proteins and M equal to the number of hit-only
  proteins.  \code{adjMat} should have row and column names corresponding to
  the proteins in the experiment.  An entry of "1" in the ith row and jth
  column of \code{adjMat} corresponds to bait protein i finding protein j as a
  hit.  All other entries should be 0. }

  \item{VBs}{\code{VBs} is an optional vector of viable baits.}

  \item{VPs}{\code{VPs} is an optional vector of viable prey.}

  \item{unrecip}{By default set to 1 so that unreciprocated bait-bait edges
  are treated as present.  If set to 0, unreciprocated bait-bait edges will be
  treated as absent.}  
}
\details{

A BH-complete subgraph with n bait nodes and m hit-only nodes for AP-MS data is defined as a subgraph for which all n*(n-1)+nm directed edges exist.  A maximal BH-complete subgraph is a BH-complete subgraph which is not contained in any other BH-complete subgraph.

If \code{VBs} and/or \code{VPs} are not specified, then by default \code{VBs}
will be assigned the set of baits that detect at least one prey and \code{VPs}
the set of prey that are detected by at least one bait.

By default, unreciprocated bait-bait observations will be
treated as present.  If \code{unrecip} is set to 0, they will be treated as
absent.  If the sensitivity of the AP-MS technology is believed to be less
than the specificity, then it is suggested that \code{unrecip}=1.

This function calls \code{maxCliques} from the RBGL package.

}
\value{
A list of length one named 'maxCliques' which is itself a list of 
character vectors containing the names of the elements in the cliques.
}

\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).

}
\author{Denise Scholtens}

\seealso{\code{\link{mergeComplexes}},\code{\link{findComplexes}}}
\examples{

data(apEX)
PCMG0 <- bhmaxSubgraph(apEX)
PCMG1 <- mergeComplexes(PCMG0,apEX,sensitivity=.7,specificity=.75)
}
\keyword{array}
