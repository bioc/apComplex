\name{TAP}
\alias{TAP}
\docType{data}
\title{Tandem Affinity Purification (TAP) Data from Gavin et al. (2002)}
\description{
A matrix of the TAP data published by Gavin, et al. (2002).
}
\usage{data(TAP)}
\details{

The 455 rows of \code{TAP} correspond to bait proteins in the TAP experiment and are named accordingly.  The 1364 columns of \code{TAP} correspond to proteins found as hits in the TAP purifications and are named accordingly.  The first 455 column names are the same as the 455 row names.  Bait proteins can also be found as hits by other baits, hence their inclusion as columns in \code{TAP}.

An entry of "1" in the ith row and jth column of \code{TAP} corresponds to bait protein i finding protein j as a hit.  All other entries are "0".  

This data is also available in Supplementary Table S1 of Gavin et al. at http://www.nature.com.

}
\source{
Gavin, et al.  Functional organization of the yeast proteome by systematic analysis of protein complexes.  Nature 415, 141-147 (2002).
}
%\references{
%  ~~ possibly secondary sources and usages ~~
%}
\examples{
data(TAP)

}
\keyword{datasets}
