\name{LCdelta}
\alias{LCdelta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes change in LxC measure}
\description{
Computes the change in the LxC measure for AP-MS protein data when two protein complex estimates are combined into one complex.    
}
\usage{
LCdelta(comp1, comp2, cMat, dataMat, baitList, simMat, mu, alpha, beta, wsVal = 2e+07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{comp1}{Column index in \code{cMat}.}
  \item{comp2}{Column index in \code{cMat}.}
  \item{cMat}{Current protein complex membership estimate affiliation matrix.}
  \item{dataMat}{Adjacency matrix of bait-hit data from an AP-MS experiment.  Rows correspond to baits and columns to hits.}
  \item{baitList}{A vector of the names of the proteins used as baits.}
  \item{simMat}{An optional square matrix with entries between 0 and 1.  Rows and columns correspond to the proteins in the experiment, and should be reported in the same order as the columns of dataMat.  Higher values in this matrix are interpreted to mean higher similarity for protein pairs.}
  \item{mu}{Parameter specification equal to log((1-specificitiy)/specificity).}
  \item{alpha}{Parameter specification equal to log(sensitivity/(1-sensitivity)).}
  \item{beta}{Optional additional parameter for the weight to give data in \code{simMat} in the logistic regression model.}
  \item{wsVal}{Workspace value to be used for computing Fisher's exact test.}
}
\details{

The Protein Complex Membership Graph (PCMG) algorithm for AP-MS data described by Scholtens and Gentleman (2004) uses a two-component measure of protein complex estimate quality, namely P=LxC.  Columns in \code{cMat} represent individual complex estimates.  The PCMG algorithm works by starting with a maximal BH-complete subgraph estimate of \code{cMat}, and then improves the estimate by combining columns.  

When proposing combinations of columns \code{comp1} and \code{comp2} in the PCMG estimate \code{cMat}, the proposal is accepted if the output from LCdelta is greater than zero.   
}
\value{
The numeric value of the change in P=LxC when columns \code{comp1} and \code{comp2} in \code{cMat} are combined into one column.  
}
\references{Scholtens, D. and Gentleman, R.  Making sense of high throughput protein-protein interaction data. (2004).}
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

 %~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{bhmaxSubgraph}},\code{\link{mergeComplexes}},\code{\link{findComplexes}}}
\examples{

#to be added

}
\keyword{graphs}

