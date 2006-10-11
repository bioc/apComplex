#revised LCdelta function using lists not matrices
#will only calculate for combination


LCjoinadjBin <- function(x,comp,adjMat,bNames,simMat,mu,alpha,Beta){

       y <- unique(c(x,comp))
       
       adjBinC <- adjBinFUN(y,adjMat=adjMat,bNames=bNames,mu=mu,alpha=alpha,Beta=Beta,simMat=simMat)

       adjBinC
}


LCjoinLchange <- function(x,comp,complexes,adjMat,ccMat,bNames,simMat,mu,alpha,Beta){

	      y <- unique(c(x,comp))
	      yB <- intersect(y,bNames)
	      
	      nB <- length(yB)
	      nH <- length(y) - nB 
	

	      exCadjMat <- ccMat[yB,y]
	      diag(exCadjMat) <- 1
	
	      LchangeC <-  sum((1-exCadjMat)*
		      (adjMat[yB,y]*alpha-
		      log(1+exp(mu+alpha+Beta*simMat[yB,y])) + 
		      log(1+exp(mu+Beta*simMat[yB,y]))))
	  	      
 
	      LchangeC

}

LCjoinfisher <- function(x,comp,adjMat,bNames,nMax,wsVal=20000000){


       y <- unique(c(x,comp))
       fisherC <-
       fisherFUN(y,adjMat=adjMat,bNames=bNames,nMax=nMax,wsVal=wsVal)
       fisherC	     

}
