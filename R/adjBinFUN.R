#x is a character vector of node names in a complex

adjBinFUN <- function(x,adjMat,bNames,mu,alpha,Beta,simMat){

	  xB <- intersect(x,bNames)	
	  
	  nB <- length(xB)
	  nH <- length(x)-nB

	  temp <- matrix(adjMat[xB,x],ncol=(nB+nH))
	  sim <- matrix(simMat[xB,x],ncol=(nB+nH))

	  adjBin <- sum(temp*(mu+alpha+Beta*sim)-
		 log(1+exp(mu+alpha+Beta*sim)))-
		 nB*(mu+alpha+Beta-log(1+exp(mu+alpha+Beta)))


	  X <- sum(temp)-nB

	  cX <- lgamma(nB*(nB+nH-1)+1) - lgamma(X+1) - lgamma(nB*(nB+nH-1)-X+1)

	  #compute penalty contribution
	  lK <- cX + adjBin	  
	  lK
}
