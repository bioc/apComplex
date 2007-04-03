
#adjMat = 1 for bait-bait entries

fisherFUN <- function(x,adjMat,bNames,nMax,wsVal=2e7){

	  
	  xB <- intersect(x,bNames)	
	  
	  nB <- length(xB)
	  nH <- length(x)-nB

	  if(nB+nH<nMax){

	  temp <- matrix(adjMat[xB,x],ncol=(nB+nH))
	  rownames(temp) <- xB
	  colnames(temp) <- x

	  bh1 <- colSums(temp)
	  bh1[xB] <- bh1[xB] - 1
	  bh0 <- colSums(1-temp)

	  if(length(bh0)>15) wsVal <- 2e9

	  ans <- log(fisher.test(rbind(bh0,bh1),workspace=wsVal)$p.value)

	  } else ans <- NA

	  ans
}


fisherFUN2 <- function(x,adjMat,VBPs,VBOs,VPOs,nMax,wsVal=2e7){

	  xVBP <- intersect(x,VBPs)
	  xVBO <- intersect(x,VBOs)
	  xVPO <- intersect(x,VPOs)	
	  
	  nVBP <- length(xVBP)
	  nVBO <- length(xVBO)
	  nVPO <- length(xVPO)


	  if(length(x)<nMax){

	  temp <- matrix(adjMat[c(xVBP,xVBO),c(xVBP,xVPO)],ncol=(nVBP+nVPO))
	  rownames(temp) <- c(xVBP,xVBO)
	  colnames(temp) <- c(xVBP,xVPO)

	  bh1 <- colSums(temp)
	  bh1[xVBP] <- bh1[xVBP] - 1
	  bh0 <- colSums(1-temp)

	  if(length(bh0)>15) wsVal <- 2e9

	  ans <- log(fisher.test(rbind(bh0,bh1),workspace=wsVal)$p.value)

	  } else ans <- NA

	  ans
}
