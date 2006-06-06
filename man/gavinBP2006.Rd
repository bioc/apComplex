\name{gavinBP2006}
\alias{gavinBP2006}
\alias{gavinBP2006graph}

\docType{data}
\title{Tandem Affinity Purification (TAP) Data from Gavin et al. (2006)}
\description{
TAP data published by Gavin, et al. (2006).
}
\usage{
data(gavinBP2006)
}
\details{

\code{gavinBP2006} is a matrix of the TAP data published by Gavin, et al. (2006).  The 1752 rows correspond to bait proteins and the 2551 columns correspond to
proteins found as hits in the TAP experiment and are named accordingly.
The first 1752 column names are the same as the 1752 row names; bait proteins
can also be found as hits by other baits, hence their inclusion as columns in
\code{gavinBP2006}.  An entry of "1" in the ith row and jth column of
\code{gavinBP2006}
indicates that bait protein i found protein j as a hit.  All other entries are
"0".  There are a total of 19105 "1" entries in the matrix, corresponding to
the 19105 comemberships detected in the experiment.

These data are available in the IntAct repository under Gavin et al - 2006.

}
\source{
Gavin, et al.  Functional Organization of the Yeast Proteome by
Systematic Analysis of the Yeast Complex.
Nature 440, 631-636 (2002).
}

\seealso{\code{\link{yTAP}},\code{\link{MBMEcTAP}},\code{\link{SBMHcTAP}},\code{\link{UnRBBcTAP}},\code{\link{HMSPCI}}
}

\examples{
data(gavinBP2006)
which(gavinBP2006[1,]==1)

}
\keyword{datasets}
