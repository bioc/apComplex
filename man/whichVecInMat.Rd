\name{whichVecInMat}
\alias{whichVecInMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Which Vector in Matrix}
\description{
This functions checks if a vector is identical to or strictly less than at least one column in a matrix, and returns the indices of the respective columns.
}
\usage{
whichVecInMat(x, mat, compare = "equal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector.}
  \item{mat}{A matrix with number of rows equal to the length of \code{x}.}
  \item{compare}{A character equal to "equal" (default) or "less".}
}
\details{
 If \code{compare} is set to "equal", the function checks to see if \code{x} is identically equal to at least one column in \code{mat} and records the indices of the matching columns.

If \code{compare} is set to "less", the function checks to see if all elements of \code{x} are less than all elements of at least one column in \code{mat} and records the corresponding column indices.

}
\value{
A vector of column names of \code{mat}, or indices if the columns are not named, that are identical to or greater than \code{x}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{reduceMat}},\code{\link{reduceMatNew}},\code{\link{vecInMat}}}
\examples{

a <- matrix(c(1,0,1,1,1,0),nrow=3)
x1 <- c(1,0,1)
x2 <- c(1,1,1)
x3 <- c(1,0,0)

whichVecInMat(x1,a)
whichVecInMat(x2,a)
whichVecInMat(x3,a,compare="less")


}
\keyword{array}

