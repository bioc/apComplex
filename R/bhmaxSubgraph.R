
#function decompose a bait-hit adjacency matrix (made symmetric for baits)
#into complex memberships

#adjMat has dimensions N by (N+M) corresponding to N baits and M hits
#adjMat is named with row and column names corresponding to proteins

#this function uses another function called reduceMat 
#any columns that are subsets of other columns are eliminated


bhmaxSubgraph <- function(adjMat){

	#record the order of the columns of adjMat so the order of the
	#rows of the affiliation matrix will match the original column order
	rowOrder <- colnames(adjMat)

	N <- dim(adjMat)[1]
	M <- dim(adjMat)[2] - N

	#make diagonal entries equal to 1
	diag(adjMat) <- 1
	
	#first find baits that only have hit pairs  
	#will put these in at the end

	hitComps <- which(rowSums(adjMat[,1:N])==1)
	baitComps <- c(1:N)
	if(length(hitComps)>0) 	baitComps <- c(1:N)[-hitComps]
	

	Nbait <- length(baitComps)

	#reorder so that hit only complexes are last
	adjMat <- adjMat[c(baitComps,hitComps),
				c(baitComps,hitComps,(N+1):(N+M))]

	compMat <- NULL

	for (n in 1:Nbait){

	newComp <- as.matrix(c(rep(0,n-1),adjMat[n,n:(N+M)]))

	changeThese <- which(compMat[n,]==1)
	n2c <- length(changeThese)

	if(n2c>0){

	newVecs <- NULL
	lose <- rep(FALSE,n2c)

	for (i in 1:n2c){
		
		maybeChange <- compMat[,changeThese[i]]

		if(sum(maybeChange[(n+1):(N+M)]>newComp[(n+1):(N+M)])>0){

		lose[i] <- TRUE		

		newComp1 <- maybeChange
		newComp1[n] <- 0

		newComp2 <- c(maybeChange[1:n],
		pmin(maybeChange[(n+1):(N+M)],newComp[(n+1):(N+M)]))

		newVecs <- cbind(newVecs,cbind(newComp1,newComp2))
		
		
		}

	}
	newComp <- as.matrix(cbind(newVecs,newComp))
	if(sum(lose)>0) compMat <- as.matrix(compMat[,-changeThese[lose]])

	}

	compMat <- cbind(compMat,newComp)
	
	}

	print("reducing matrix")
	compMat <- reduceMat(compMat,compare="less")

	if (Nbait<N) compMat <- cbind(compMat,t(adjMat[(Nbait+1):N,]))

	newcolnames <- paste("bhmax",1:dim(compMat)[2],sep="")
	colnames(compMat) <- newcolnames
	compMat <- compMat[rowOrder,]

	return(compMat)
}

