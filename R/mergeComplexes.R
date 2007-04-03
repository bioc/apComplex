#revised function for mergeComplexes using lists rather than matrices

#bhmax is a list of length one names 'maxCliques'
#bhmax is output from the bhmaxSubgraph function
#bhmax$maxCliques is a list of character vectors containing the clique members

#this uses the following functions: adjBinFUN, fisherFUN,
#LCjoinadjBin, LCjoinfisher, and LCjoinLchange

mergeComplexes <- 
function(bhmax,adjMat,simMat=NULL,sensitivity=.75,specificity=.995,Beta=0,commonFrac=2/3,wsVal 
= 2e7){

	stopifnot("maxCliques" %in% names(bhmax))
	#initial complex estimates
	complexes <- bhmax$maxCliques

	#set parameters
	bNames <- rownames(adjMat)
	diag(adjMat) <- 1
	mu <- log((1-specificity)/specificity)
	alpha <- log(sensitivity/(1-sensitivity))-mu
	N <- dim(adjMat)[1]
	M <- dim(adjMat)[2]-N

	#make complex co-membership matrix
	ccMat <- adjMat
	ccMat[bNames,bNames] <-
		pmax(adjMat[bNames,bNames],t(adjMat[bNames,bNames]))


	#make simMat with entries 0 and diagonal 1 if simMat not specified
	if(is.null(simMat)){
	simMat <- matrix(0,N,N+M)
	diag(simMat) <- 1
	rownames(simMat) <- rownames(adjMat)
	colnames(simMat) <- colnames(adjMat)}

	#calculate adj binomial and fisher test for initial complexes
	print("calculating initial penalty terms")

	adjBin <-unlist(lapply(complexes,FUN=adjBinFUN,adjMat=adjMat,
	bNames=bNames,mu=mu,alpha=alpha,Beta=Beta,simMat=simMat))

	fisher <- unlist(lapply(complexes,FUN=fisherFUN,adjMat=adjMat,
	bNames=bNames,nMax=20,wsVal=wsVal))

	#start looking at combinations
	print("looking at complex combinations")
	i <- 1 
	K <- length(complexes)

	keepgoing <- i < K
   
 
	while(keepgoing){

	keepgoing2 <- TRUE

	while(keepgoing2){

	#for complex under consideration, 
	#narrow combination candidates to those with common members

	commonFracFUN <- function(x) length(intersect(x,complexes[[i]]))>floor(length(complexes[[i]])*commonFrac)

	testL <- lapply(complexes,FUN=commonFracFUN)
	testset <- which(unlist(testL))
	testset <- testset[-which(testset==i)]
	Ktemp <- length(testset)

	#if there are candidates, then test to see 
	#if penalized likelihood increases when combined	
	if(Ktemp>0){

	#find adjusted binomial for two complexes
	lK2 <- adjBin[i] + adjBin[testset]

	#find adjusted binomial for joined complex
	lK1adjBin <- unlist(lapply(complexes[testset],FUN=LCjoinadjBin,comp=complexes[[i]],adjMat=adjMat,bNames=bNames,simMat=simMat,mu=mu,alpha=alpha,Beta=Beta))

	#find change in likelihood when adding new edges to graph
	Lchange <-
	unlist(lapply(complexes[testset],FUN=LCjoinLchange,comp=complexes[[i]],complexes=complexes,adjMat=adjMat,ccMat=ccMat,bNames=bNames,simMat=simMat,mu=mu,alpha=alpha,Beta=Beta))

	#find total change in likelihood
	LCInc1 <- lK1adjBin+Lchange-lK2

	#for combinations where it could reasonably make a difference,
	#look at change attributable to fisher's exact test
	lK1fisher <- rep(0,length(testset))
	dofisher <- which(LCInc1>-20)
	
	if(length(dofisher)>0){

	lKfisher <-
	unlist(lapply(complexes[testset[dofisher]],FUN=LCjoinfisher,comp=complexes[[i]],adjMat=adjMat,bNames=bNames,nMax=20,wsVal=wsVal))
	

	lK1fisher[dofisher] <- lKfisher - fisher[i] - fisher[testset[dofisher]]
	}
	
	
	#add in fisher's exact component
	LCIncs <- LCInc1+lK1fisher

	#some fisher tests will result in NA if too many proteins
	#replace these with likelihood change without fisher component
	LCIncs[which(is.na(LCIncs))] <- LCInc1[which(is.na(LCIncs))]
	
	same <- sum(LCIncs>0)>0

	#if any combinations increase the likihood, then find the maximum
	#and make the combination
	
	if(same){
		wm <- which.max(LCIncs)
		thisone <- testset[wm]

		combo <- unique(c(complexes[[i]],complexes[[thisone]]))
		complexes[[i]] <- combo		

		#remove combined complex
		complexes <- complexes[-thisone]

		#make corresponding changes to adjBin and fisher vectors
	        adjBin[i] <- lK1adjBin[wm]
		adjBin <- adjBin[-thisone]

		fisher[i] <- lK1fisher[wm]	
		if(is.na(lKfisher[wm])) fisher[i] <- NA
		fisher <- fisher[-thisone]

		#make changes in complex comembership matrix
		comboB <- intersect(combo,bNames)
		ccMat[comboB,combo] <- 1


		K <- length(complexes)
		if(thisone<i) i <- i-1
		

	}
	}else same <- FALSE
	keepgoing2 <- same
	}

	
	i <- i+1
	keepgoing <- i < K

	print(paste("i",i,"K",K))

}
nC <- length(complexes)
names(complexes) <- paste("Complex",1:nC,sep="")
return(complexes)

}


#new function adjusting for viability status


mergeComplexes2 <- 
function(bhmax,adjMat,VBs=NULL,VPs=NULL,simMat=NULL,sensitivity=.75,specificity=.995,Beta=0,commonFrac=2/3,wsVal 
= 2e7){

	stopifnot("maxCliques" %in% names(bhmax))
	if(!is.null(VBs)) stopifnot(all(VBs %in% rownames(adjMat)))
	if(!is.null(VPs)) stopifnot(all(VPs %in% colnames(adjMat)))

	#create viable bait and prey sets if not specified
	if(is.null(VBs)) VBs <- rownames(adjMat)[rowSums(adjMat)>0]
	if(is.null(VPs)) VPs <- rownames(adjMat)[rowSums(adjMat)>0]

	VBPs <- intersect(VBs,VPs)
	VBOs <- setdiff(VBs,VBPs)
	VPOs <- setdiff(VPs,VBPs)

	#initial complex estimates
	complexes <- bhmax$maxCliques

	adjMat <- adjMat[c(VBPs,VBs),c(VBPs,VPs)]

	#set parameters
	adjMat[VBPs,VBPs] <- 1
	mu <- log((1-specificity)/specificity)
	alpha <- log(sensitivity/(1-sensitivity))-mu

	#make complex co-membership matrix
	ccMat <- adjMat
	ccMat[VBPs,VBPs] <- pmax(adjMat[VBPs,VBPs],t(adjMat[VBPs,VBPs]))


	#make simMat with entries 0 and diagonal 1 if simMat not specified
	if(is.null(simMat)){
	simMat <- matrix(0,dim(adjMat)[1],dim(adjMat)[2])
	rownames(simMat) <- rownames(adjMat)
	colnames(simMat) <- colnames(adjMat)
	simMat[VBPs,VBPs] <- 1
	}

	#calculate adj binomial and fisher test for initial complexes
	print("calculating initial penalty terms")

	adjBin <-unlist(lapply(complexes,FUN=adjBinFUN2,adjMat=adjMat,
	VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,
	mu=mu,alpha=alpha,Beta=Beta,simMat=simMat))

	fisher <- unlist(lapply(complexes,FUN=fisherFUN2,adjMat=adjMat,
	VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,
	nMax=20,wsVal=wsVal))

	#start looking at combinations
	print("looking at complex combinations")
	i <- 1 
	K <- length(complexes)

	keepgoing <- i < K
   
 
	while(keepgoing){

	keepgoing2 <- TRUE

	while(keepgoing2){

	#for complex under consideration, 
	#narrow combination candidates to those with common members

	commonFracFUN <- function(x) length(intersect(x,complexes[[i]]))>floor(length(complexes[[i]])*commonFrac)

	testL <- lapply(complexes,FUN=commonFracFUN)
	testset <- which(unlist(testL))
	testset <- testset[-which(testset==i)]
	Ktemp <- length(testset)

	#if there are candidates, then test to see 
	#if penalized likelihood increases when combined	
	if(Ktemp>0){

	#find adjusted binomial for two complexes
	lK2 <- adjBin[i] + adjBin[testset]

	#find adjusted binomial for joined complex
	lK1adjBin <- unlist(lapply(complexes[testset],FUN=LCjoinadjBin2,comp=complexes[[i]],adjMat=adjMat,VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,simMat=simMat,mu=mu,alpha=alpha,Beta=Beta))

	#find change in likelihood when adding new edges to graph
	Lchange <-
	unlist(lapply(complexes[testset],FUN=LCjoinLchange2,comp=complexes[[i]],complexes=complexes,adjMat=adjMat,ccMat=ccMat,VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,simMat=simMat,mu=mu,alpha=alpha,Beta=Beta))

	#find total change in likelihood
	LCInc1 <- lK1adjBin+Lchange-lK2

	#for combinations where it could reasonably make a difference,
	#look at change attributable to fisher's exact test
	lK1fisher <- rep(0,length(testset))
	dofisher <- which(LCInc1>-20)
	
	if(length(dofisher)>0){

	lKfisher <-
	unlist(lapply(complexes[testset[dofisher]],FUN=LCjoinfisher2,comp=complexes[[i]],adjMat=adjMat,VBPs=VBPs,VBOs=VBOs,VPOs=VPOs,nMax=20,wsVal=wsVal))
	

	lK1fisher[dofisher] <- lKfisher - fisher[i] - fisher[testset[dofisher]]
	}
	
	
	#add in fisher's exact component
	LCIncs <- LCInc1+lK1fisher

	#some fisher tests will result in NA if too many proteins
	#replace these with likelihood change without fisher component
	LCIncs[which(is.na(LCIncs))] <- LCInc1[which(is.na(LCIncs))]
	
	same <- sum(LCIncs>0)>0

	#if any combinations increase the likihood, then find the maximum
	#and make the combination
	
	if(same){
		wm <- which.max(LCIncs)
		thisone <- testset[wm]

		combo <- unique(c(complexes[[i]],complexes[[thisone]]))
		complexes[[i]] <- combo		

		#remove combined complex
		complexes <- complexes[-thisone]

		#make corresponding changes to adjBin and fisher vectors
	        adjBin[i] <- lK1adjBin[wm]
		adjBin <- adjBin[-thisone]

		fisher[i] <- lK1fisher[wm]	
		if(is.na(lKfisher[wm])) fisher[i] <- NA
		fisher <- fisher[-thisone]

		#make changes in complex comembership matrix
		comboB <- intersect(combo,bNames)
		ccMat[comboB,combo] <- 1


		K <- length(complexes)
		if(thisone<i) i <- i-1
		

	}
	}else same <- FALSE
	keepgoing2 <- same
	}

	
	i <- i+1
	keepgoing <- i < K

	print(paste("i",i,"K",K))

}
nC <- length(complexes)
names(complexes) <- paste("Complex",1:nC,sep="")
return(complexes)

}
