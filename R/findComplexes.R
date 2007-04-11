

# a function to run the entire algorithm at once

findComplexes <- 
function(adjMat,VBs=NULL,VPs=NULL,simMat=NULL,sensitivity=.75,specificity=.995, 
Beta=0,commonFrac=2/3,wsVal = 2e7){

	##find number of baits and number of hits

        !is.null(colnames(adjMat)) || stop("Columns of adjMat must be named")
        !is.null(rownames(adjMat))|| stop("Rows of adjMat must be named")
        
        if(!is.null(VBs)) stopifnot(all(VBs %in% rownames(adjMat)))
        if(!is.null(VPs)) stopifnot(all(VPs %in% colnames(adjMat)))
        
        #create viable bait and prey sets if not specified
        if(is.null(VBs)) VBs <- rownames(adjMat)[rowSums(adjMat)>0]
        if(is.null(VPs)) VPs <- colnames(adjMat)[colSums(adjMat)>0]
        
        VBPs <- intersect(VBs,VPs)
        VBOs <- setdiff(VBs,VBPs) 
        VPOs <- setdiff(VPs,VBPs) 
    
	##set parameters for logistic regression model

	mu <- log((1-specificity)/specificity)
	alpha <- log(sensitivity/(1-sensitivity))-mu

	##create simMat of zeroes with diagonal of ones if one is not specified

	if(is.null(simMat)) {
		simMat <- matrix(0,dim(adjMat)[1],dim(adjMat)[2])
		diag(simMat) <- 1
		colnames(simMat) <- colnames(adjMat)
		rownames(simMat) <- rownames(adjMat)
	}

	##find maximal BH-complete subgraphs for initial 
	##protein complex membership graph estimate

	print("Finding Initial Maximal BH-complete Subgraphs")
	PCMG <- bhmaxSubgraph(adjMat,VBs=VBs,VPs=VPs,unrecip=1*(sensitivity<specificity))

	##combine complexes using LC measure

	#put PCMG in order by number of baits in complex
	numBaitsFUN <- function(x) sum(x %in% VBPs)
	numBaits <- unlist(lapply(PCMG$maxCliques,FUN=numBaitsFUN))
	
	baitOrder <- order(numBaits,decreasing=TRUE)
	PCMGo <- PCMG
	PCMGo$maxCliques <- PCMGo$maxCliques[baitOrder]
	
	#merge complex estimates using LCdelta criteria
	print("Combining Complex Estimates")
	PCMG2 <-
	mergeComplexes(PCMGo,adjMat=adjMat,VBs=VBs,VPs=VPs,simMat=simMat,
			Beta=Beta,sensitivity=sensitivity,
			specificity=specificity,commonFrac=commonFrac,
			wsVal = wsVal)

	return(PCMG2)

}



 
