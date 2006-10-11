\name{sortComplexes}
\alias{sortComplexes}

\title{Sort complex estimates}
\description{
Sorts complexes recorded in PCMG list into three separate
lists containing MBME, SBMH, and UnRBB complexes.
}
\usage{
sortComplexes(PCMG, adjMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PCMG}{Current PCMG estimate}
  \item{adjMat}{Adjacency matrix of bait-hit data from an AP-MS experiment.  Rows correspond to baits and columns to hits.}
}
\details{

MBME complexes contain multiple bait proteins and multiple edges.  SBMH
complexes contain one bait and a collection of hit-only proteins.  UnRBB
complexes contain only two baits (no hit-only proteins) that are connected by
an unreciprocated edge.

}
\value{
A list of lists representing the MBME, SBMH, and UnRBB complex
estimates. 
}

\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).

}
\author{Denise Scholtens}

\seealso{\code{\link{findComplexes}}}

\examples{

data(apEX)
PCMG2 <- findComplexes(apEX,sensitivity=.7,specificity=.75)
sortComplexes(PCMG2,apEX)

}
\keyword{graphs}

