\name{reduceMat}
\alias{reduceMat}

\title{Remove matrix columns}
\description{
This function eliminates columns that are identically equal to or strictly less than other columns in a matrix.
}
\usage{
reduceMat(mat, compare = "equal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{\code{mat} is a matrix.}
  \item{compare}{\code{compare} is a character equal to "equal" (default) or "less".}
}
\details{
If \code{compare} is set to "equal", sets of identically equal columns are found.  The first of these columns is kept in the matrix and the rest are eliminated.

If \code{compare} is set to "less", then a column is removed from the matrix if all of its entries are less than all entries of another column in the matrix. 

}
\value{
A matrix with columns eliminated according to the criteria specified in \code{compare}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

 %~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{reduceMatNew}},\code{\link{vecInMat}},\code{\link{whichVecInMat}}}
\examples{
a <- matrix(c(1,0,1,1,1,0,1,0,1,1,0,0),nrow=3)
reduceMat(a)
reduceMat(a,compare="less")

}
\keyword{array}

