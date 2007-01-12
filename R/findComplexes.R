

# a function to run the entire algorithm at once


findComplexes <- 
function(adjMat,simMat=NULL,sensitivity=.75,specificity=.995, 
Beta=0,commonFrac=2/3,wsVal = 2e7){

	##find number of baits and number of hits

	N <- dim(adjMat)[1]
	M <- dim(adjMat)[2]-N

	baits <- rownames(adjMat)

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
	numBaitsFUN <- function(x) sum(x %in% baits)
	numBaits <- unlist(lapply(PCMG$maxCliques,FUN=numBaitsFUN))
	
	baitOrder <- order(numBaits,decreasing=TRUE)
	PCMGo <- PCMG
	PCMGo$maxCliques <- PCMGo$maxCliques[baitOrder]
	
	#merge complex estimates using LCdelta criteria
	print("Combining Complex Estimates")
	PCMG2 <-
	mergeComplexes(PCMGo,adjMat=adjMat,simMat=simMat,Beta=Beta,
			sensitivity=sensitivity,specificity=specificity,commonFrac=commonFrac,wsVal = wsVal)

	return(PCMG2)

}



 
