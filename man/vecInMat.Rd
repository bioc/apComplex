\name{vecInMat}
\alias{vecInMat}
\alias{whichVecInMat}

\title{Vector in Matrix}
\description{
These functions check if a vector is identical to or less than at
least one column in a matrix and return either a logical value or the indices
of the respective columns.
}
\usage{
vecInMat(x, mat, compare = "equal")
whichVecInMat(x, mat, compare="equal")
}

\arguments{
  \item{x}{A vector.}
  \item{mat}{A matrix with number of rows equal to the length of \code{x}.}
  \item{compare}{A character equal to "equal" (default) or "less".}
}
\details{
If \code{compare} is set to "equal", the function checks to see if \code{x} is identically equal to at least one column in \code{mat}.

If \code{compare} is set to "less", the function checks to see if all elements of \code{x} are less than all elements of at least one column in \code{mat}.


}
\value{
\code{vecInMat} returns a logical value, indicating if \code{x} is equal to
less than at least one column in \code{mat}.

\code{whichVecInMat} returns a vector of column names of \code{mat}, or
indices if the clumns are not named, that are identical to or greater than \code{x}.
}

\author{Denise Scholtens}


\seealso{\code{\link{reduceMat}}}
\examples{
a <- matrix(c(1,0,1,1,1,0),nrow=3)
x1 <- c(1,0,1)
x2 <- c(1,1,1)
x3 <- c(1,0,0)

vecInMat(x1,a)
vecInMat(x2,a)
vecInMat(x3,a,compare="less")

whichVecInMat(x1,a)
whichVecInMat(x2,a)
whichVecInMat(x3,a,compare="less")

}
\keyword{array}

