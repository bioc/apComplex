\name{HMSPCI}
\alias{HMSPCI}
\docType{data}
\title{High-Throughput Mass Spectromic Protein Complex Identification (HMSPCI) Data from Ho, et al. (2002)}
\description{
A matrix of the HMSPCI data published by Ho, et al. (2002).
}
\usage{data(HMSPCI)}
\details{

The 493 rows of \code{HMSPCI} correspond to bait proteins in the HMSPCI experiment and are named accordingly.  The 1578 columns of \code{HMSPCI} correspond to proteins found as hits in the HMSPCI purifications and are named accordingly.  The first 493 column names are the same as the 493 row names.  Bait proteins can also be found as hits by other baits, hence their inclusion as columns in \code{HMSPCI}.

An entry of "1" in the ith row and jth column of \code{HMSPCI} corresponds to bait protein i finding protein j as a hit.  All other entries are "0".  

This data is also available at http://www.mdsp.com/yeast.

}
\source{
Ho, et al.  Systematic identification of ptoein complexes in Saccharomyces cerevisiae by mass spectrometry.  Nature 415, 180-183 (2002).
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(HMSPCI)

}
\keyword{datasets}
