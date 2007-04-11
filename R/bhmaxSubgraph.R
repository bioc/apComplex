
#function to find bhmaxSubgraphs from a bait-hit adjacency matrix

#by default, unreciprocated bait-bait edges will be treated as observed

#adjMat has dimensions N by (N+M) corresponding to N baits and M hits
#adjMat is named with row and column names corresponding to proteins

#this function uses 'maxClique' from RBGL

bhmaxSubgraph <- function(adjMat,VBs=NULL,VPs=NULL,unrecip=1){

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
	
	allProts <- c(VBPs,VBOs,VPOs)
	nProts <- length(allProts)

	#reorder adjMat rows and columns
	adjMat <- adjMat[c(VBPs,VBOs),c(VBPs,VPOs)]
	diag(adjMat) <- 0

	#for VPOs, insert an edge if found by a common bait
	#for VBOs, insert ad edge if they find a common prey
	
	adjMatAppend <- matrix(0,nProts,nProts)
	rownames(adjMatAppend) <- allProts
	colnames(adjMatAppend) <- allProts
	adjMatAppend[VBPs,VBPs] <- adjMat[VBPs,VBPs]
	adjMatAppend[VBPs,VPOs] <- adjMat[VBPs,VPOs]
	adjMatAppend[VBOs,VBPs] <- adjMat[VBOs,VBPs]
	adjMatAppend[VBOs,VPOs] <- adjMat[VBOs,VPOs]
	adjMatAppend[VPOs,VPOs] <- (1*(t(adjMat) %*% adjMat > 0))[VPOs,VPOs]
	adjMatAppend[VBOs,VBOs] <- (1*(adjMat %*% t(adjMat) > 0))[VBOs,VBOs]
	diag(adjMatAppend) <- 0

	#make corresponding undirected graph
	g <- as(adjMatAppend,"graphNEL")
	ug <- ugraph(g)

	#find maximal cliques
	mcs <- maxClique(ug)

	#remove cliques containing only VPOs or VBOs
	vbovpoFUN <- function(x) all(x %in% VBOs) | all(x %in% VPOs) 
	vbovpoc <- which(unlist(lapply(mcs$maxCliques,FUN=vbovpoFUN)))
	if(length(vbovpoc)>0) mcs$maxCliques <- mcs$maxCliques[-vbovpoc]

	#remove cliques with 1 members -- since the diagonal=0
	#this shouldn't happen, but it seems to - bug in maxClique?
	mem1 <- which(unlist(lapply(mcs$maxCliques,FUN=length))==1)
	if(length(mem1)>0) mcs$maxCliques <- mcs$maxCliques[-mem1]
	
	mcs
}

