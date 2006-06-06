\name{kroganBPMat2006}
\alias{kroganBPMat2006}


\docType{data}
\title{Tandem Affinity Purification (TAP) Data from Krogan et al. (2006)}
\description{
TAP data published by Krogan, et al. (2006).
}
\usage{
data(kroganBPMat2006)
}
\details{

\code{kroganBPMat2006} is a matrix of the TAP data published by Krogan, et al. (2006).  The 2264 rows correspond to bait proteins and the 5361 columns correspond to
proteins found as hits in the TAP experiment and are named accordingly.
The first 2264 column names are the same as the 2264 row names; bait proteins
can also be found as hits by other baits, hence their inclusion as columns in
\code{kroganBPMat2006}.  An entry of "1" in the ith row and jth column of
\code{kroganBPMat2006}
indicates that bait protein i found protein j as a hit.  All other entries are
"0".  There are a total of  "1" entries in the matrix, corresponding to
the 63360 comemberships detected in the experiment.

These data were obtained from the Primary Source - 

}
\source{
  N Krogan et al. Global Landscape of Protein Complexes in the Yeast,
  S. Cerevisiae. Nature, 440:667-643, 2006.
}

\seealso{\code{\link{yTAP}},\code{\link{MBMEcTAP}},\code{\link{SBMHcTAP}},\code{\link{UnRBBcTAP}},\code{\link{HMSPCI}}
}

\examples{
data(kroganBPMat2006)
which(kroganBPMat2006[1,]==1)

}
\keyword{datasets}
