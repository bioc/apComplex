\name{findComplexes}
\alias{findComplexes}

\title{Estimate a Protein Complex Membership Graph (PCMG) using protein complex comembership data from AP-MS technology}
\description{
Performs all steps in the local modeling algorithm described by Scholtens and
Gentleman (2004) and Scholtens, Vidal, and Gentleman (submitted), beginning with an adjacency matrix recording bait-hit AP-MS data.
}
\usage{
findComplexes(adjMat, simMat = NULL, sensitivity = 0.75, specificity =
0.995, Beta = 0, wsVal=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjMat}{Adjacency matrix of bait-hit data from an AP-MS experiment.  Rows correspond to baits and columns to hits.}
  \item{simMat}{An optional square matrix with entries between 0 and 1.  Rows and columns correspond to the proteins in the experiment, and should be reported in the same order as the columns of \code{adjMat}.  Higher values in this matrix are interpreted to mean higher similarity for protein pairs.}
  \item{sensitivity}{Believed sensitivity of AP-MS technology.}
  \item{specificity}{Believed specificity of AP-MS technology.}
  \item{Beta}{Optional additional parameter for the weight to give data
    in \code{simMat} in the logistic regression model.}
  \item{wsVal}{A numeric. This is the value assigned as the work-space
    in the call to fisher.test}
}
\details{
\code{findComplexes} performs all steps in the complex estimation algorithm using the apComplex package functions \code{bhmaxSubgraph}, \code{LCdelta}, and \code{mergeComplexes}.  These steps can also be performed separately by the user.  
}
\value{
An affiliation matrix representing the estimated protein complex memberships. 
}
\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).
}
\author{Denise Scholtens}

\seealso{\code{\link{bhmaxSubgraph}},code{\link{LCdelta}},\code{\link{mergeComplexes}}}
\examples{

data(apEX)
PCMG2 <- findComplexes(apEX,sensitivity=.7,specificity=.75)


}
\keyword{graphs}

