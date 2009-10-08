\name{yNameTAP complexes}
\alias{yNameTAP}


\docType{data}
\title{TAP data complex estimates using standard gene names}
\description{
Affiliation matrices with rows corresponding to proteins and columns
corresponding to complexes. 
}
\usage{
data(yNameTAP)
}
\details{
These are the results from an analysis of the TAP data (Gavin et al., 2002) by
Scholtens and Gentleman (2004) and Scholtens, Vidal, and Gentleman
(submitted). These estimates were constructed using \code{findComplexes}
with a sensitivity 
parameter of .75, specificity of .995, and Beta=-0.2 for externally derived
similarity measure based on Gene Ontology cellular component annotation (see
Scholtens and Gentleman (2004)).


\code{yNameTAP} contains 260 multi-bait-multi-edge complex estimates.

}
\source{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).

}
\references{
Gavin, et al. Functional organization of the yeast proteome by
systematic analysis of protein complexes.  Nature 415, 141-147 (2002). 

}

\seealso{\code{\link{TAP}}, \code{\link{TAPgraph}},
  \code{\link{yTAP}}, \code{findComplexes}}


\examples{
data(yNameTAP)
yNameTAP[1:10,1:3]


}
\keyword{datasets}
