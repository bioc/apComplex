\name{TAP complexes}
\alias{MBMEcTAP}
\alias{SBMHcTAP}
\alias{UnRBBcTAP}


\docType{data}
\title{TAP data complex estimates}
\description{
Affiliation matrices with rows corresponding to proteins and columns corresponding to complexes.
}
\usage{
data(MBMEcTAP)
data(SBMHcTAP)
data(UnRBBcTAP)
}
\details{
These are the results from an analysis of the TAP data (Gavin et al., 2002) using the Protein Complex Membership Graph (PCMG) algorithm proposed by Scholtens and Gentleman (2004).  These estimates were constructed using a sensitivity parameter of .75, specificity of .995, and external data based on Gene Ontology cellular component annotation.

\code{MBMEcTAP} contains multi-bait-multi-edge complex estimates.  \code{SBMHcTAP} contains single-bait-multi-hit complex estimates.  \code{UnRBBcTAP} contains unreciprocated bait-bait complex estimates.


}
\source{
Scholtens, D. and Gentleman, R.  Making sense of high throughput protein-protein interaction data. (2004).
}
\references{
Gavin, et al. Functional organization of the yeast proteome by systematic analysis of protein complexes.  Nature 415, 141-147 (2002).

}
\examples{
data(MBMEcTAP)

}
\keyword{datasets}
