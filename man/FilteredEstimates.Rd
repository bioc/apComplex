\name{AP-MS complexes}
\alias{gavin06FilteredEstimates}
\alias{krogan06FilteredEstimates}

\docType{data}
\title{AP-MS data complex estimates}
\description{
List of complex estimates after filtering baits prone to systematic bias

}
\usage{
data(gavin06FilteredEstimates)
data(krogan06FilteredEstimates)
}
\details{
\code{gavin06FilteredEstimates} contains results from an analysis of the AP-MS
data published by Gavin et al. (2006). 
These estimates were constructed using \code{findComplexes} with a sensitivity
parameter of .70, specificity of .999, and commonFrac=0.5  

\code{krogan06FilteredEstimates} contains results from an analysis of the AP-MS
data published by Krogan et al. (2006). 
These estimates were constructed using \code{findComplexes} with a sensitivity
parameter of .70, specificity of .999, and commonFrac=1/3. 

Both sets of estimates are reported as lists of length three, corresponding to
MBME, SBMH, and UnRBB complex estimate types (see Scholtens et al., 2005).
Each of the three elements contains a list of character vectors of
estimated complex members. 
}

\source{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).

}
\references{
Gavin, et al. (2006)

Krogan, et al. (2006)

}

\seealso{\code{\link{gavinBP2006}},
  \code{\link{kroganBPMat2006}},
  \code{findComplexes}
}

\examples{
data(gavin06FilteredEstimates)
lapply(gavin06FilteredEstimates,FUN=length)

data(krogan06FilteredEstimates)
lapply(krogan06FilteredEstimates,FUN=length)


}
\keyword{datasets}
