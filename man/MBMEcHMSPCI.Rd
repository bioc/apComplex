\name{HMSPCI complexes}
\alias{MBMEcHMSPCI}
\alias{SBMHcHMSPCI}
\alias{UnRBBcHMSPCI}


\docType{data}
\title{HMSPCI data complex estimates}
\description{
Affiliation matrices with rows corresponding to proteins and columns
corresponding to complexes. 
}
\usage{
data(MBMEcHMSPCI)
data(SBMHcHMSPCI)
data(UnRBBcHMSPCI)

}
\details{
These are the results from an analysis of the HMS-PCI data (Ho et al., 2002)
described by Scholtens and Gentleman (2004) and Scholtens, Vidal, and
Gentleman (submitted).  These 
estimates were constructed using \code{findComplexes} with a sensitivity
parameter of .75, specificity of .99, and Beta=-0.2 for externally derived
similarity measure based on Gene Ontology cellular component annotation (see
Scholtens and Gentleman (2004)).

\code{MBMEcHMSPCI} contains 242 multi-bait-multi-edge complex estimates.
\code{SBMHcHMSPCI} contains 437 single-bait-multi-hit complex estimates.
\code{UnRBBcHMSPCI} contains 329 unreciprocated bait-bait complex
estimates.  These data are also available at
http://www.bioconductor.org/Docs/Papers/2003/apComplex. 

}
\source{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).
}
\references{
Ho, et al. Systematic identification of protein complexes in
Saccharomyces cerevisiae by mass spectrometry.  Nature 415, 180-183
(2002). 
}

\seealso{\code{\link{HMSPCI}},\code{\link{HMSPCIgraph}},\code{findComplexes}}


\examples{
data(MBMEcHMSPCI)
MBMEcHMSPCI[1:4,1:4]
which(MBMEcHMSPCI[,"MBME1"]==1)

}
\keyword{datasets}
