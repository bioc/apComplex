\name{mergeComplexes}
\alias{mergeComplexes}
\title{Iteratively combine columns in initial PCMG estimate}
\description{
Repeatedly applies the function \code{LCdelta} to make combinations of columns in the affiliation matrix representing the protein complex membership graph (PCMG) for AP-MS data.
}
\usage{
mergeComplexes(PCMG, adjMat, simMat = NULL, sensitivity = 0.75, specificity = 0.995, beta = 0)
}
\arguments{
  \item{PCMG}{Current PCMG estimate}
  \item{adjMat}{Adjacency matrix of bait-hit data from an AP-MS experiment.  Rows correspond to baits and columns to hits.}
  \item{simMat}{An optional square matrix with entries between 0 and 1.  Rows and columns correspond to the proteins in the experiment, and should be reported in the same order as the columns of dataMat.  Higher values in this matrix are interpreted to mean higher similarity for protein pairs.}
  \item{sensitivity}{Believed sensitivity of AP-MS technology.}
  \item{specificity}{Believed specificity of AP-MS technology.}
  \item{beta}{Optional additional parameter for the weight to give data in \code{simMat} in the logistic regression model.}
}
\details{

The Protein Complex Membership Graph (PCMG) algorithm for AP-MS data described by Scholtens and Gentleman (2004) uses a two-component measure of protein complex estimate quality, namely P=LxC.  Columns in \code{cMat} represent individual complex estimates.  The PCMG algorithm works by starting with a maximal BH-complete subgraph estimate of \code{cMat}, and then improves the estimate by combining columns.  

When proposing combinations of columns \code{comp1} and \code{comp2} in the PCMG estimate \code{cMat}, the proposal is accepted if the output from LCdelta (the log of LxC) is greater than zero.   \code{mergeComplexes} performs all column combinations until no more combinations result in an output from LCdelta greater than zero.

}
\value{
An affiliation matrix representing the estimated PCMG.  The number of rows and row labels of the matrix will be the same as \code{adjMat}.  The number of columns will be less than or equal to the number of columns in \code{adjMat}.
}

\references{Scholtens, D. and Gentleman, R.  Making sense of high throughput protein-protein interaction data. (2004).}
\author{Denise Scholtens}
\seealso{\code{\link{LCdelta}},\code{\link{findComplexes}}}
\examples{

# to be added

}
\keyword{graphs}

