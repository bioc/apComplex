\name{mergeComplexes}
\alias{mergeComplexes}

\title{Iteratively combine columns in initial PCMG estimate}
\description{
Repeatedly applies the function \code{LCdelta} to make combinations of columns in the affiliation matrix representing the protein complex membership graph (PCMG) for AP-MS data.
}
\usage{
mergeComplexes(PCMG, adjMat, simMat = NULL, sensitivity = 0.75, specificity = 0.995, Beta = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PCMG}{Current PCMG estimate}
  \item{adjMat}{Adjacency matrix of bait-hit data from an AP-MS experiment.  Rows correspond to baits and columns to hits.}
  \item{simMat}{An optional square matrix with entries between 0 and 1.  Rows and columns correspond to the proteins in the experiment, and should be reported in the same order as the columns of \code{adjMat}.  Higher values in this matrix are interpreted to mean higher similarity for protein pairs.}
  \item{sensitivity}{Believed sensitivity of AP-MS technology.}
  \item{specificity}{Believed specificity of AP-MS technology.}
  \item{Beta}{Optional additional parameter for the weight to give data in \code{simMat} in the logistic regression model.}
}
\details{

The local modeling algorithm for AP-MS data described by Scholtens and
Gentleman (2004) and Scholtens, Vidal, and Gentleman (submitted) uses a
two-component measure of protein complex estimate quality, namely P=LxC.
Columns in \code{cMat} represent individual complex estimates.  The algorithm
works by starting with a maximal BH-complete subgraph estimate of \code{cMat},
and then improves the estimate by combining columns such that P=LxC increases.  

When proposing combinations of columns \code{comp1} and \code{comp2} in
\code{cMat}, the proposal is accepted if the output from \code{LCdelta} (the
log of LxC) is greater than zero.   \code{mergeComplexes} performs all column
combinations until no more combinations result in an output from
\code{LCdelta} greater than zero.

}
\value{
An affiliation matrix representing the estimated PCMG.  The number of rows and
the row labels of the matrix will be the same as \code{adjMat}.  The number of columns will be less than or equal to the number of columns in \code{adjMat}.
}

\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Submitted.

}
\author{Denise Scholtens}

\seealso{\code{\link{LCdelta}},\code{\link{bhmaxSubgraph}},\code{\link{findComplexes}}}

\examples{

data(apEX)
PCMG0 <- bhmaxSubgraph(apEX)
PCMG1 <- mergeComplexes(PCMG0,apEX,sensitivity=.7,specificity=.75)

}
\keyword{graphs}

