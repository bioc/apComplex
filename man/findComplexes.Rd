\name{findComplexes}
\alias{findComplexes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate a Protein Complex Membership Graph (PCMG) using protein complex comembership data from AP-MS technology}
\description{
Performs all steps in the PCMG algorithm described by Scholtens and Gentleman (2004), beginning with an adjacency matrix recording bait-hit AP-MS data.
}
\usage{
findComplexes(adjMat, simMat = NULL, sensitivity = 0.75, specificity = 0.995, beta = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjMat}{Adjacency matrix of bait-hit data from an AP-MS experiment.  Rows correspond to baits and columns to hits.}
  \item{simMat}{An optional square matrix with entries between 0 and 1.  Rows and columns correspond to the proteins in the experiment, and should be reported in the same order as the columns of dataMat.  Higher values in this matrix are interpreted to mean higher similarity for protein pairs.}
  \item{sensitivity}{Believed sensitivity of AP-MS technology.}
  \item{specificity}{Believed specificity of AP-MS technology.}
  \item{beta}{Optional additional parameter for the weight to give data in \code{simMat} in the logistic regression model.}
}
\details{
\code{findComplexes} performs all steps in the PCMG algorithm using the apComplex package functions \code{bhmaxSubgraph}, \code{LCdelta}, and \code{mergeComplexes}.  These steps can also be performed separately by the user.  
}
\value{
An affiliation matrix representing the estimated PCMG. 
}
\references{Scholtens, D. and Gentleman, R.  Making sense of high throughput protein-protein interaction data. (2004).}
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{bhmaxSubgraph}},code{\link{LCdelta}},\code{\link{mergeComplexes}}}
\examples{

# to be added

}
\keyword{graphs}

