##propose complex combinations and compare to LxC measure

##PCMG is an initial estimate for the PCMG, presumably obtained 
##by finding bhMaximal subgraphs using bhmaxSubgraph

mergeComplexes <- function(PCMG,adjMat,simMat=NULL,sensitivity=.75,specificity=.995,Beta=0){

	   
	bNames <- rownames(adjMat)
	diag(adjMat) <- 1

	mu <- log((1-specificity)/specificity)
	alpha <- log(sensitivity/(1-sensitivity))-mu

	N <- dim(adjMat)[1]
	M <- dim(adjMat)[2]-N

	#make simMat with entries 0 and diagonal 1 if simMat not specified
	if(is.null(simMat)){
	simMat <- matrix(0,N,N+M)
	diag(simMat) <- 1
	rownames(simMat) <- rownames(adjMat)
	colnames(simMat) <- colnames(adjMat)}

	i <- 1 
	K <- dim(PCMG)[2]

	keepgoing <- i < K
   
 
	while(keepgoing){

	keepgoing2 <- TRUE

	while(keepgoing2){

	testset <- which(colSums(PCMG[,i]*PCMG)>0)
	testset <- testset[-which(testset==i)]
	Ktemp <- length(testset)

	
	if(Ktemp>0 & Ktemp!=0){
	LCIncs <- rep(0,Ktemp)

	for (m in 1:Ktemp){
	
	LCIncs[m] <- LCdelta(i,testset[m],PCMG,dataMat=adjMat,
				baitList=bNames,simMat=simMat,
				mu=mu,alpha=alpha,Beta=Beta)
	
	
	
	}
	
	same <- sum(LCIncs>0)>0
	
	if(same){
		thisone <- testset[which.max(LCIncs)]
		PCMG[,i] <- pmax(PCMG[,i],PCMG[,thisone])
		PCMG <- PCMG[,-thisone]
		K <- dim(PCMG)[2]
		if(thisone<i) i <- i-1
		

	}
	}else same <- FALSE
	keepgoing2 <- same
	}

	
	i <- i+1
	keepgoing <- i < K
}
nC <- dim(PCMG)[2]
colnames(PCMG) <- paste("Complex",1:nC,sep="")
return(PCMG)
}
