\name{Krogan complexes}
\alias{MBMEcKrogan}

\docType{data}
\title{Krogan data complex estimates}
\description{
Affiliation matrices with rows corresponding to proteins and columns corresponding to complexes.
}
\usage{
data(MBMEcKrogan)

}
\details{
These are the results from an analysis of the AP-MS data (Krogan et al., 2004).
These 
estimates were constructed using \code{findComplexes} with a sensitivity
parameter of .75, specificity of .99, and Beta=-0.2 for externally derived
similarity measure based on Gene Ontology cellular component annotation (see
Scholtens and Gentleman (2004)).

\code{MBMEcHMSPCI} contains 82 multi-bait-multi-edge complex estimates.

}
\source{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).
}
\references{
High-Definition Macromolecular Composition of Yeast RNA-Processing
Complexes; Molecular Cell, Vol 13, 225-239, 30 January 2004
}


\examples{
data(MBMEcKrogan)
MBMEcKrogan[1:4,1:4]

}
\keyword{datasets}
