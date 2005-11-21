\name{Krogan}
\alias{Krogan}


\docType{data}
\title{High-Definition Macromolecular Composition of Yeast
  RNA-Processing Complexes (2004)}
\description{
AP-MS data published by Krogan et al. (2004)
}
\usage{
data(Krogan)

}
\details{

\code{Krogan} is a matrix of the AP-MS data published by Krogan et al (2004).
The 153 rows correspond to bait proteins and the 485 columns correspond to
proteins found as hits in the AP-MS experiment and are named accordingly.
The first 153 column names are the same as the 153 row names; bait proteins
can also be found as hits by other baits, hence their inclusion as columns in
\code{Krogan}.  An entry of "1" in the ith row and jth column of \code{Krogan}
indicates that bait protein i found protein j as a hit.  All other entries are
"0".  There are a total of 1132 "1" entries in the matrix, corresponding to
the 1132 comemberships detected in the experiment.
 
These data are available at http://www.molecule.org/cgi/content/full/13/2/225/DC1/.

}
\source{
Krogan, et al.  High-Definition Macromolecular Composition of Yeast
RNA-Processing Complexes; Molecular Cell, Vol 13, 225-239, 30 January 2004.
}
\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Submitted.
}

\seealso{\code{\link{MBMEcKrogan}}}


\examples{
data(Krogan)

}
\keyword{datasets}
