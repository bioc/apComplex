\name{bhmaxSubgraph}
\alias{bhmaxSubgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find maximal BH-complete subgraph}
\description{
Given an adjacency matrix of bait-hit AP-MS protein data, this function finds the maximal BH-complete subgraphs and reports them as an affiliation matrix.
}
\usage{
bhmaxSubgraph(adjMat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{adjMat}{\code{adjMat} is an N by (N+M) adjacency matrix with N equal to the number of bait proteins and M equal to the number of hit-only proteins.  \code{adjMat} should have row and column names corresponding to the proteins in the experiment.  An entry of "1" in the ith row and jth column of \code{adjMat} corresponds to bait protein i finding protein j as a hit.  All other entries should be 0. }
}
\details{

A BH-complete subgraph with n bait nodes and m hit-only nodes for AP-MS data is defined as a subgraph for which all n*(n-1)+nm directed egdes exist.  A maximal BH-complete subgraph is a BH-complete subgraph which is not contained in any other BH-complete subgraph.

This function assumes that the first N columns of \code{adjMat} form a symmetric matrix.  For actual AP-MS data, this may not be true.  It is up to the user to decide how best to symmetrize this part of the matrix.

}
\value{
A matrix with (N+M) rows and a number of columns equal to the number of maximal BH-complete subgraphs.  An entry of "1" in the ith row and jth column indicates membership of the ith protein in the jth maximal BH-complete subgraph.
}

\references{Scholtens, D. and Gentleman R.  Making Sense of High Throughput Protein-Protein Interaction Data.}
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

 %~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }
\examples{

a <- matrix(c(1,1,1,1,1,0,1,0,1,1,1,0,1,1,1),nrow=3)
rownames(a) <- c("B1","B2","B3")
colnames(a) <- c("B1","B2","B3","H1","H2")

bhmaxSubgraph(a)

}
\keyword{array}
