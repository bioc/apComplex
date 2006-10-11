
#function to find bhmaxSubgraphs from a bait-hit adjacency matrix

#by default, unreciprocated bait-bait edges will be treated as observed

#adjMat has dimensions N by (N+M) corresponding to N baits and M hits
#adjMat is named with row and column names corresponding to proteins

#this function uses 'maxClique' from RBGL

bhmaxSubgraph <- function(adjMat,unrecip=1){

	!is.null(colnames(adjMat)) || stop("Columns of adjMat must be named")
	!is.null(rownames(adjMat))|| stop("Rows of adjMat must be named")

	Nb <- dim(adjMat)[1]
	Nh <- dim(adjMat)[2] - Nb
	
	baits <- rownames(adjMat)
	prey <- setdiff(colnames(adjMat),baits)

	#make sure first Nb colnames match rownames
	adjMat <- adjMat[baits,c(baits,prey)]
	diag(adjMat) <- 0

	#for hit-only proteins, insert an edge if found by a common bait
	adjMatAppend <- rbind(adjMat,matrix(0,nrow=Nh,ncol=(Nb+Nh)))
	rownames(adjMatAppend)[(Nb+1):(Nb+Nh)] <- prey
	adjMatAppend[prey,prey] <- (1*(t(adjMat) %*% adjMat > 0))[prey,prey]
	diag(adjMatAppend) <- 0

	#make corresponding undirected graph
	g <- as(adjMatAppend,"graphNEL")
	ug <- ugraph(g)

	#find maximal cliques
	mcs <- maxClique(ug)

	#remove prey-only cliques
	poFUN <- function(x) sum(x %in% baits)==0
	poc <- which(unlist(lapply(mcs$maxCliques,FUN=poFUN)))
	if(length(poc)>0) mcs$maxCliques <- mcs$maxCliques[-poc]

	#remove cliques with 1 members -- since the diagonal=0
	#this shouldn't happen, but it seems to - bug in maxClique?
	mem1 <- which(unlist(lapply(mcs$maxCliques,FUN=length))==1)
	if(length(mem1)>0) mcs$maxCliques <- mcs$maxCliques[-mem1]
	
	mcs
}

