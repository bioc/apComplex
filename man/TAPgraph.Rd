\name{TAPgraph}
\alias{TAPgraph}
\docType{data}
\title{Tandem Affinity Purification (TAP) Data from Gavin et al. (2002)}
\description{
A graph of the TAP data published by Gavin, et al. (2002).
}
\usage{data(TAPgraph)}
\details{
A graphNEL object with 1364 nodes and 3420 directed edges.  The 1364 nodes represent the 455 bait proteins in the TAP data and the additional 909 proteins that were found as hits in the purifications but never used as baits.  Directed edges connect bait proteins to the proteins found as hits in their purifications. 

This data is also available in Supplementary Table S1 of Gavin et al. at http://www.nature.com. 
}
\source{
Gavin, et al.  Functional organization of the yeast proteome by systematic analysis of protein complexes.  Nature 415, 141-147 (2002).
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(TAPgraph)
%## maybe str(TAPgraph) ; plot(TAPgraph) ...
}
\keyword{datasets}
