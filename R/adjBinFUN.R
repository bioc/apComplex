#x is a character vector of node names in a complex

adjBinFUN <- function(x,adjMat,VBPs,VBOs,VPOs,mu,alpha,Beta,simMat){

	  xVBP <- intersect(x,VBPs)
	  xVBO <- intersect(x,VBOs)
	  xVPO <- intersect(x,VPOs)	
	  
	  nVBP <- length(xVBP)
	  nVBO <- length(xVBO)
	  nVPO <- length(xVPO)

	  temp <- matrix(adjMat[c(xVBP,xVBO),c(xVBP,xVPO)],ncol=(nVBP+nVPO))
	  sim <- matrix(simMat[c(xVBP,xVBO),c(xVBP,xVPO)],ncol=(nVBP+nVPO))

	  adjBin <- sum(temp*(mu+alpha+Beta*sim)-
		 log(1+exp(mu+alpha+Beta*sim)))-
		 nVBP*(mu+alpha+Beta-log(1+exp(mu+alpha+Beta)))


	  X <- sum(temp)-nVBP

	  cX <- lgamma(nVBP*(nVBP+nVPO+nVBO-1)+nVBO*nVPO+1) - 
	     lgamma(X+1) - lgamma(nVBP*(nVBP+nVPO+nVBO-1)+nVBO*nVPO-X+1)

	  #compute penalty contribution
	  lK <- cX + adjBin	  
	  lK
}
