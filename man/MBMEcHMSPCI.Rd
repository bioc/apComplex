\name{HMSPCI complexes}
\alias{MBMEcHMSPCI}
\alias{SBMHcHMSPCI}
\alias{UnRBBcHMSPCI}


\docType{data}
\title{HMSPCI data complex estimates}
\description{
Affiliation matrices with rows corresponding to proteins and columns corresponding to complexes.
}
\usage{
data(MBMEcHMSPCI)
data(SBMHcHMSPCI)
data(UnRBBcHMSPCI)

}
\details{
These are the results from an analysis of the HMS-PCI data (Ho et al., 2002) using the Protein Complex Membership Graph (PCMG) algorithm proposed by Scholtens and Gentleman (2004).  These estimates were constructed using a sensitivity parameter of .75, specificity of .99, and external data based on Gene Ontology cellular component annotation.

\code{MBMEcHMSPCI} contains multi-bait-multi-edge complex estimates.  \code{SBMHcHMSPCI} contains single-bait-multi-hit complex estimates.  \code{UnRBBcHMSPCI} contains unreciprocated bait-bait complex estimates.


}
\source{
Scholtens, D. and Gentleman, R.  Making sense of high throughput protein-protein interaction data. (2004).
}
\references{
Ho, et al. Systematic identification of protein complexes in Saccharomyces cerevisiae by mass spectrometry.  Nature 415, 180-183 (2002).
}
\examples{
data(MBMEcHMSPCI)
}
\keyword{datasets}
