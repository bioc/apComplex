\name{apEX}
\alias{apEX}
\alias{apEXG}
\docType{data}
\title{Example data set for apComplex package}
\description{
A matrix and corresponding graph of AP-MS purifications from a small hypothetical experiment.
}
\usage{
data(apEX)
data(apEXG)
}
\details{

\code{apEX} is a matrix of hypothetical AP-MS data.  The 4 rows correspond to bait proteins and the 8 columns correspond to
proteins found as hits in the hypothetical experiment and are named accordingly.
The first 4 column names are the same as the 4 row names; bait proteins
can also be found as hits by other baits, hence their inclusion as columns in
\code{apEX}.  An entry of "1" in the ith row and jth column of \code{apEX}
indicates that bait protein i found protein j as a hit.  All other entries are
"0".  There are a total of 9 "1" entries in the matrix, corresponding to 9 comemberships detected in the experiment.

\code{apEXG} is a \code{graphNEL} object in which 8 nodes
represent proteins and 9 directed edges represent comemberships detected in
the purification.  Each directed edge originates at the bait and ends at
the hit.  Each edge in \code{apEXG} corresponds to an entry of "1" in
the \code{apEX} matrix. 

These data are used in the \code{apComplex} vignette.

}
%\source{
%}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(apEX)
apEX

data(apEXG)
adj(apEXG,"P1")
}
\keyword{datasets}
