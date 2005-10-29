
# a function to run the entire algorithm at once


findComplexes <- function(adjMat,simMat=NULL,sensitivity=.75,specificity=.995, Beta=0){

	##find number of baits and number of hits

	N <- dim(adjMat)[1]
	M <- dim(adjMat)[2]-N

	##set parameters for logistic regression model

	mu <- log((1-specificity)/specificity)
	alpha <- log(sensitivity/(1-sensitivity))-mu

	##create simMat of zeroes with diagonal of ones if one is not specified

	if(is.null(simMat)) {
		simMat <- matrix(0,N,N+M)
		diag(simMat) <- 1
		colnames(simMat) <- colnames(adjMat)
		rownames(simMat) <- rownames(adjMat)
	}

	##find maximal BH-complete subgraphs for initial 
	##protein complex membership graph estimate

	print("Finding Initial Maximal BH-complete Subgraphs")
	PCMG <- bhmaxSubgraph(adjMat,unrecip=1*(sensitivity<specificity))

	##combine complexes using LC measure

	#put PCMG in order by number of baits in complex
	
	baitOrder <- order(colSums(PCMG[1:N,,drop=FALSE]),decreasing=TRUE)
	PCMGo <- PCMG[,baitOrder,drop=FALSE]
	
	#merge complex estimates using LCdelta criteria
	print("Combining Complex Estimates")
	PCMG2 <-
	mergeComplexes(PCMGo,adjMat=adjMat,simMat=simMat,Beta=Beta,
			sensitivity=sensitivity,specificity=specificity)

	return(PCMG2)

}



 
