#revised LCdelta function using lists not matrices
#will only calculate for combination


LCjoinadjBin <- function(x,comp,adjMat,bNames,simMat,mu,alpha,Beta){

       y <- unique(c(x,comp))
       
       adjBinC <- adjBinFUN(y,adjMat=adjMat,bNames=bNames,mu=mu,alpha=alpha,Beta=Beta,simMat=simMat)

       adjBinC
}

LCjoinadjBin2 <- function(x,comp,adjMat,VBPs,VBOs,VPOs,simMat,mu,alpha,Beta){

       y <- unique(c(x,comp))
       
       adjBinC <- adjBinFUN2(y,adjMat=adjMat,VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,mu=mu,alpha=alpha,Beta=Beta,simMat=simMat)

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

LCjoinLchange2 <- function(x,comp,complexes,adjMat,ccMat,VBPs,VBOs,VPOs,simMat,mu,alpha,Beta){

	      y <- unique(c(x,comp))
	      yVBPs <- intersect(y,VBPs)
	      yVBOs <- intersect(y,VBOs)
	      yVPOs <- intersect(y,VPOs)

	      exCadjMat <- ccMat[c(yVBPs,yVBOs),c(yVBPs,yVPOs)]
	      exCadjMat[yVBPs,yVBPs] <- 1
	
	      LchangeC <-  sum((1-exCadjMat)*
		      (adjMat[c(yVBPs,yVBOs),c(yVBPs,yVPOs)]*alpha-
		      log(1+exp(mu+alpha+Beta*
				simMat[c(yVBPs,yVBOs),c(yVBPs,yVPOs)])) + 
		      log(1+exp(mu+Beta*
				simMat[c(yVBPs,yVBOs),c(yVBPs,yVPOs)]))))
	  	      
 
	      LchangeC

}

LCjoinfisher <- function(x,comp,adjMat,bNames,nMax,wsVal=2e7){


       y <- unique(c(x,comp))
       fisherC <-
       fisherFUN(y,adjMat=adjMat,bNames=bNames,nMax=nMax,wsVal=wsVal)
       fisherC	     

}

LCjoinfisher2 <- function(x,comp,adjMat,VBPs,VBOs,VPOs,nMax,wsVal=2e7){


       y <- unique(c(x,comp))
       fisherC <-
       fisherFUN2(y,adjMat=adjMat,VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,nMax=nMax,wsVal=wsVal)
       fisherC	     

}
