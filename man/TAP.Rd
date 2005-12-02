\name{TAP}
\alias{TAP}
\alias{TAPgraph}

\docType{data}
\title{Tandem Affinity Purification (TAP) Data from Gavin et al. (2002)}
\description{
TAP data published by Gavin, et al. (2002).
}
\usage{
data(TAP)
data(TAPgraph)
}
\details{

\code{TAP} is a matrix of the TAP data published by Gavin, et al. (2002).  The 455 rows correspond to bait proteins and the 1364 columns correspond to
proteins found as hits in the TAP experiment and are named accordingly.
The first 455 column names are the same as the 455 row names; bait proteins
can also be found as hits by other baits, hence their inclusion as columns in
\code{TAP}.  An entry of "1" in the ith row and jth column of \code{TAP}
indicates that bait protein i found protein j as a hit.  All other entries are
"0".  There are a total of 3420 "1" entries in the matrix, corresponding to
the 3420 comemberships detected in the experiment.

\code{TAPgraph} is a \code{graphNEL} object in which 1364 nodes
represent proteins and 3420 directed edges represent comemberships detected in
the purification.  Each directed edge originates at the bait and ends at
the hit.  Each edge in \code{TAPgraph} corresponds to an entry of "1" in
the \code{TAP} matrix. 

These data are available in Supplementary Table S1 of Gavin et al. at
http://www.nature.com.

}
\source{
Gavin, et al.  Functional organization of the yeast proteome by systematic analysis of protein complexes.  Nature 415, 141-147 (2002).
}

\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).
}

\seealso{\code{\link{yTAP}},\code{\link{MBMEcTAP}},\code{\link{SBMHcTAP}},\code{\link{UnRBBcTAP}},\code{\link{HMSPCI}}
}

\examples{
data(TAP)
which(TAP["Abd1",]==1)

data(TAPgraph)
adj(TAPgraph,"Abd1")
}
\keyword{datasets}
