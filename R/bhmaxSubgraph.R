
#function to find bhmaxSubgraphs from a bait-hit adjacency matrix

#by default, unreciprocated bait-bait edges will be treated as observed

#adjMat has dimensions N by (N+M) corresponding to N baits and M hits
#adjMat is named with row and column names corresponding to proteins

#this function uses another function called reduceMat 
#any columns that are subsets of other columns are eliminated


bhmaxSubgraph <- function(adjMat,unrecip=1){

	!is.null(colnames(adjMat)) || stop("Columns of adjMat must be named")
	!is.null(rownames(adjMat))|| stop("Rows of adjMat must be named")

	Nb <- dim(adjMat)[1]
	Nh <- dim(adjMat)[2] - Nb
	
	identical(rownames(adjMat),colnames(adjMat)[1:Nb]) || stop("rownames
	and first Nb colnames of adjMat must be identical")
	
	#make adjMat[,1:Nb] symmetric
	if(unrecip==0){
	 adjMat[,1:Nb] <- pmin(adjMat[,1:Nb],t(adjMat[,1:Nb]))
	} else adjMat[,1:Nb] <- pmax(adjMat[,1:Nb],t(adjMat[,1:Nb]))

	#record the order of the columns of adjMat so the order of the
	#rows of the affiliation matrix will match the original column order
	rowOrder <- colnames(adjMat)
	hNames <- rowOrder[!rowOrder %in% rownames(adjMat)]

	#make diagonal entries equal to 1
	diag(adjMat) <- 1
	
	#first find baits that only have hit pairs  
	#will put these in at the end

	hitComps <- which(rowSums(adjMat[,1:Nb])==1)
	baitComps <- rownames(adjMat)
	if(length(hitComps)>0) 	baitComps <- baitComps[-hitComps]
	
	#now reorder by complex size - smaller first
	baitComps <- baitComps[order(rowSums(adjMat[baitComps,]))]


	Nbait <- length(baitComps)

	#reorder so that hit only complexes are last
	adjMat <- adjMat[c(baitComps,names(hitComps)),
				c(baitComps,names(hitComps),hNames)]


	
	M <- as.matrix(adjMat[1,])

	for (i in 2:Nbait){

	g <- as.matrix(c(rep(0,i-1),adjMat[i,i:(Nb+Nh)]))
	G <- NULL

	V <- which(M[i,]==1)
	n2c <- length(V)

	if(n2c>0){

	lose <- rep(FALSE,n2c)

	for (k in 1:n2c){
		
		v <- M[,V[k]]

		if(sum(v[(i+1):(Nb+Nh)]>g[(i+1):(Nb+Nh)])>0){

		lose[k] <- TRUE		

		v1 <- v
		v1[i] <- 0
		G <- cbind(G,v1)

		v2 <- c(v[1:i],pmin(v[(i+1):(Nb+Nh)],g[(i+1):(Nb+Nh)]))
		G <- cbind(G,v2)
		}

	}

	if(sum(lose)>0) M <- as.matrix(M[,-V[lose]])

	}

	G <- cbind(G,g)
	if(dim(G)[2]>1) G <- reduceMat(G,compare="less")
		
	for (g in 1:dim(G)[2]){
		
		if(dim(M)[2]>0){
		Msub <- as.matrix(M[,which(colSums(M)>=sum(G[,g]))])
					
		if(dim(Msub)[2]>0){
		if(!vecInMat(G[,g],Msub,compare="less")) M <- cbind(M,G[,g])
		} else M <- cbind(M,G[,g])
		} else M <- cbind(M,G[,g])
		
		
	}
	
	}


	if (Nbait<Nb){
	if (Nbait+1==Nb) { 
	   M <- cbind(M,as.matrix(adjMat[Nb,]))
	} else  M <- cbind(M,t(adjMat[(Nbait+1):Nb,]))
	}

	newcolnames <- paste("bhmax",1:dim(M)[2],sep="")
	colnames(M) <- newcolnames
	M <- M[rowOrder,]

	return(M)

}

