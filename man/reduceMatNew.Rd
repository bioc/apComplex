\name{reduceMatNew}
\alias{reduceMatNew}

\title{Remove matrix columns}
\description{
This function eliminates columns that are identically equal to or strictly less than other columns in a matrix.  This is an alternate version of reduceMat.  It allows a distinction between previously reduced columns of a matrix and new columns that have not yet been reduced.
}
\usage{
reduceMatNew(mat, newvec, compare = "equal")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{\code{mat} is a matrix.}
  \item{newvec}{\code{newvec} is a vector with length equal to the number of columns of \code{mat}.  \code{newvec} contains a 1 for  new column and 0 for a previously reduced column.  The order of the entries in \code{newvec} should correspond to the columns of \code{mat}.}
  \item{compare}{\code{compare} is a character equal to "equal" (default) or "less".}
}
\details{
If \code{compare} is set to "equal", then all new columns, indicated in \code{newvec}, which are identically equal to another column in the previously reduced columns of \code{mat} are removed from the matrix.  

If \code{compare} is set to "less", then a new column is removed from the matrix if all of its entries are less than all entries of another column in the matrix. 
}
\value{
A matrix with columns eliminated according to the criteria specified in \code{compare}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{reduceMat}},\code{\link{vecInMat}},\code{\link{whichVecInMat}}}
\examples{
a <- matrix(c(1,0,1,1,1,0,1,0,1,1,0,0),nrow=3)
reduceMatNew(a,newvec=c(0,0,1,1))
reduceMatNew(a,newvec=c(0,0,1,1),compare="less")

}
\keyword{array}

