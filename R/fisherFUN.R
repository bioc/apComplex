
#adjMat = 1 for bait-bait entries

fisherFUN <- function(x,adjMat,bNames,nMax,wsVal=20000000){

	  
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

	  ans <- log(fisher.test(rbind(bh0,bh1),workspace=wsVal)$p.value)

	  } else ans <- NA

	  ans
}
