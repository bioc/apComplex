\name{LCjoin}
\alias{LCjoin}
\alias{adjBinFUN}
\alias{fisherFUN}
\alias{LCjoinadjBin}
\alias{LCjoinfisher}
\alias{LCjoinLchange}
\alias{adjBinFUN2}
\alias{fisherFUN2}
\alias{LCjoinadjBin2}
\alias{LCjoinfisher2}
\alias{LCjoinLchange2}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes change in LxC measure}
\description{
Computes the change in the P=LxC measure for AP-MS protein data when two protein complex estimates are combined into one complex.    
}
%\usage{
%LCdelta(comp1, comp2, cMat, dataMat, baitList, simMat, mu, alpha, Beta, 
%wsVal = 2e+07)
%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%  \item{comp1}{Column index in \code{cMat}.}
%  \item{comp2}{Column index in \code{cMat}.}
%  \item{cMat}{Current protein complex membership estimate affiliation 
%matrix.}
%  \item{dataMat}{Adjacency matrix of bait-hit data from an AP-MS 
%experiment.  Rows correspond to baits and columns to hits.}
%  \item{baitList}{A vector of the names of the proteins used as baits.}
%  \item{simMat}{An optional square matrix with entries between 0 and 1.  
%Rows and columns correspond to the proteins in the experiment, and should 
%be reported in the same order as the columns of \code{dataMat}.  Higher 
%values in this matrix are interpreted to mean higher similarity for 
%protein pairs.}
%  \item{mu}{Parameter specification equal to 
%log((1-specificitiy)/specificity).}
%  \item{alpha}{Parameter specification equal to 
%log(sensitivity/(1-sensitivity)).}
%  \item{Beta}{Optional additional parameter for the weight to give data 
%in \code{simMat} in the logistic regression model.}
%  \item{wsVal}{Workspace value to be used for computing Fisher's exact 
%test.}
%}
\details{
These functions are used to evaluate changes in the penalized likelihood 
when two complexes are combined.  They are not meant to be directly used.

}
\value{
The numeric value of the change in P=LxC when two complexes are combined.
}
\references{
Scholtens D and Gentleman R.  Making sense of high-throughput protein-protein
interaction data.  Statistical Applications in Genetics and Molecular Biology
3, Article 39 (2004).

Scholtens D, Vidal M, and Gentleman R.  Local modeling of global interactome
networks.  Bioinformatics 21, 3548-3557 (2005).

}
\author{Denise Scholtens}
%\note{ ~~further notes~~ }

 %~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{bhmaxSubgraph}},\code{\link{mergeComplexes}},\code{\link{findComplexes}}}
\examples{

data(apEX)
PCMG0 <- bhmaxSubgraph(apEX)
PCMG1 <- mergeComplexes(PCMG0,apEX,sensitivity=.7,specificity=.75)


}
\keyword{graphs}

