#function to see if a vector x is identical to (or strictly less than) at least one of the columns in a matrix mat


vecInMat <- function(x,mat,compare="equal"){

	(length(x)==dim(mat)[1]) || stop("vector length must equal row dimension of matrix")

	if(compare=="equal") condFunc <- function(temp) sum(x != temp)==0
	if(compare=="less") condFunc <- function(temp) sum(x > temp)==0
	if(compare=="greater") condFunc <- function(temp) sum(x < temp)==0



	n <- dim(mat)[2]

	test <- FALSE
	i <- 1
	while (!test & i<=n){

	temp <- mat[,i]
	
	test <- condFunc(temp)

	i <- i+1
	}

	return(test)
}

