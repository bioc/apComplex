#function to report which of the columns in a matrix mat x is equal to (or strictly less than)


#will return the column name, or the index if the matrix is unnamed

whichVecInMat <- function(x,mat,compare="equal"){

	A <- colnames(mat)

	if(compare=="equal") compFun <- function(y) sum(x != y)==0
	if(compare=="less") compFun <- function(y) sum(x > y)==0
	if(compare=="greater") compFun <- function(y) sum(x < y)==0

	v <- apply(mat,2,FUN = compFun)

	w <- which(v==1)
	if(!is.null(A)) w <- A[w]

	return(w)
}

