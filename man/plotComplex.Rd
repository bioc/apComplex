\name{plotComplex}
\alias{plotComplex}

\title{Render complex estimates}
\description{
plotComplex renders complex estimates from the apComplex algorithm using Rgraphviz.
}
\usage{
plotComplex(complexMembers,g,VBs,VPs,geneName=FALSE,baitColor="yellow",preyColor="white",recipLineColor="red",unrecipBBLineColor="blue",unrecipBPLineColor="gray",y="neato")

}

\arguments{
  \item{complexMembers}{A character vector of proteins composing a complex estimate.}
  \item{g}{An object of class graph, the full bait-prey graph of AP-MS data
  used in analysis. complexMembers must be a subset of the node names of g.}
  \item{VBs}{A vector of viable baits used in the AP-MS experiment.}
  \item{VPs}{A vector of viable prey used in the AP-MS experiment.}
  \item{geneName}{A logical indicating whether or not nodes should be plotted
  with common gene names as labels rather than systematic names.}
  \item{baitColor}{Color of bait nodes.}
  \item{preyColor}{Color of prey nodes.}
  \item{recipLineColor}{Color of edges connecting baits which both detected
  each other as prey}.
  \item{unrecipBBLineColor}{Color of edges connecting baits in which one bait
  finds the other as prey but not vice versa.}
  \item{unrecipBPLineColor}{Color of edges extending from baits to proteins
  that were only used as prey, hence reciprocity is not possible.}
  \item{y}{Layout of plot}

}
\details{

This is a simple function for plotting complex estimates resulting from the
apComplex algorithm.  Giving the upcoming changes in Rgraphviz, it will likely
be changed substantially.

}
\value{
A plotted graph of the complex estimate subgraph.
}

\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).

}
\author{Denise Scholtens}

\seealso{\code{\link{findComplexes}}}

\examples{

data(apEX)
data(apEXG)
PCMG2 <- findComplexes(apEX,sensitivity=.7,specificity=.75)
PCMG2sorted <- sortComplexes(PCMG2,apEX)

VBs <- rownames(apEX)
VPs <- setdiff(colnames(apEX),VBs)

plotComplex(PCMG2sorted$MBME[[1]],g=apEXG,VBs=VBs, VPs=VPs)

}
\keyword{graphs}

