\name{HMSPCIgraph}
\alias{HMSPCIgraph}
\docType{data}
\title{High-Throughput Mass Spectromic Protein Complex Identification (HMSPCI) Data from Ho, et al. (2002)}
\description{
A graph of the HMSPCI data published by Ho, et al. (2002).
}
\usage{data(HMSPCIgraph)}
\details{
A graphNEL object with 1578 nodes and 3687 directed edges.  The 1578 nodes represent the 493 bait proteins in the HMSPCI data and the additional 1085 proteins that were found as hits in the purifications but never used as baits.  Directed edges connect bait proteins to the proteins found as hits in their purifications. 

This data is also available at http://www.mdsp.com/yeast.

}
\source{
Ho, et al.  Systematic identification of ptoein complexes in Saccharomyces cerevisiae by mass spectrometry.  Nature 415, 180-183 (2002).
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(HMSPCIgraph)

}
\keyword{datasets}
