\name{HMSPCI}
\alias{HMSPCI}
\alias{HMSPCIgraph}

\docType{data}
\title{High-Throughput Mass Spectromic Protein Complex Identification (HMS-PCI) Data from Ho, et al. (2002)}
\description{
HMS-PCI data published by Ho, et al. (2002).
}
\usage{
data(HMSPCI)
data(HMSPCIgraph)
}
\details{

\code{HMSPCI} is a matrix of the HMS-PCI data published by Ho, et al. (2002).
The 493 rows correspond to bait proteins and the 1578 columns correspond to
proteins found as hits in the HMSPCI experiment and are named accordingly.
The first 493 column names are the same as the 493 row names; bait proteins
can also be found as hits by other baits, hence their inclusion as columns in
\code{HMSPCI}.  An entry of "1" in the ith row and jth column of \code{HMSPCI}
indicates that bait protein i found protein j as a hit.  All other entries are
"0".  There are a total of 3687 "1" entries in the matrix, corresponding to
the 3687 comemberships detected in the experiment.
 
\code{HMSPCIgraph} is a \code{graphNEL} object in which 1578 nodes
represent proteins and 3687 directed edges represent comemberships detected in
the purification.  Each directed edge originates at the bait and ends at
the hit.  Each edge in \code{HMSPCIgraph} corresponds to an entry of "1" in
the \code{HMSPCI} matrix. 

These data are available at http://www.mdsp.com/yeast.

}
\source{
Ho, et al.  Systematic identification of protein complexes in Saccharomyces cerevisiae by mass spectrometry.  Nature 415, 180-183 (2002).
}
\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).
}

\seealso{\code{\link{MBMEcHMSPCI}},\code{\link{SBMHcHMSPCI}},\code{\link{UnRBBcHMSPCI}},\code{\link{TAP}}
}

\examples{
data(HMSPCI)
which(HMSPCI["YAL015C",]==1)

data(HMSPCIgraph)
adj(HMSPCIgraph,"YAL015C")
}
\keyword{datasets}
