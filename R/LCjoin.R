#revised LCdelta function using lists not matrices
#will only calculate for combination


LCjoinadjBin <- function(x,comp,adjMat,VBPs,VBOs,VPOs,simMat,mu,alpha,Beta){

       y <- unique(c(x,comp))
       
       adjBinC <- adjBinFUN(y,adjMat=adjMat,VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,mu=mu,alpha=alpha,Beta=Beta,simMat=simMat)

       adjBinC
}

LCjoinLchange <- function(x,comp,complexes,adjMat,ccMat,VBPs,VBOs,VPOs,simMat,mu,alpha,Beta){

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

LCjoinfisher <- function(x,comp,adjMat,VBPs,VBOs,VPOs,nMax,wsVal=2e7){


       y <- unique(c(x,comp))
       fisherC <-
       fisherFUN(y,adjMat=adjMat,VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,nMax=nMax,wsVal=wsVal)
       fisherC	     

}
