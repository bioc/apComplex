\name{yTAP}
\alias{yTAP}
\docType{data}
\title{yTAP Complexes published Gavin, et al. (2002).}
\description{
A list of the 232 'yTAP' complexes published by Gavin, et al. (2002) using their TAP data.
}
\usage{data(yTAP)}
\details{
Each element of the yTAP list contains a vector of the proteins reported as part of the complex, as well as the functional annotation category published by Gavin, et al. (2002).

These data are available in Supplementary Table S3 of Gavin, et al. (2002) at
http://www.nature.com and at http://yeast.cellzome.com.

}
\source{
Gavin, et al.  Functional organization of the yeast proteome by systematic analysis of protein complexes.  Nature 415, 141-147 (2002).

}

\seealso{\code{\link{TAP}},\code{\link{TAPgraph}},\code{\link{MBMEcTAP}},\code{\link{SBMHcTAP}},\code{\link{UnRBBcTAP}},\code{\link{HMSPCI}}
}


\examples{
data(yTAP)
yTAP[1]
}
\keyword{datasets}
