#function to eliminate columns that are subsets of other columns
#this is a modified version of reduceMat
#it requires a vector that says which columns have not been compared 
#to the previously reduced matrix

#new vec contains a 1 for a new column and 0 for an previously reduced one
#the order of the entries in newvec should correspond to the columns of mat

reduceMatNew <- function(mat,newvec,compare="equal"){

	oldCols <- which(newvec==0)
	newCols <- which(newvec==1)

	matKeep <- as.matrix(mat[,oldCols])
	
	if(compare=="equal"){
		
		lose <- NULL
		for (i in 1:length(newCols)){

		Vec <- mat[,newCols[i]]
		test <- vecInMat(Vec,matKeep,compare="equal")
		
		if(test) lose <- c(lose,newCols[i])
		}
		if(length(lose)>0) mat <- mat[,-lose]

	}

	if(compare=="less"){

		lose <- NULL

		for (i in 1:length(newCols)){

		Vec <- mat[,newCols[i]]
		test <- vecInMat(Vec,matKeep,compare="less")
		
		if(test){
			lose <- c(lose,newCols[i])
		}

		}

		if(!is.null(lose)) {mat <- mat[,-lose]; newvec <- newvec[-lose]}

		#check to see if columns in the previously reduced matrix
		#are strictly less than a new column

		if(sum(newvec)>0){
	
			lose <- NULL
			oldCols <- which(newvec==0)
			newCols <- which(newvec==1)
			matNew <- as.matrix(mat[,newCols])

			for (j in 1:length(oldCols)){
			temp <- mat[,oldCols[j]]
			test <- vecInMat(temp,matNew,compare="less")
			if(test) lose <- c(lose,oldCols[j])
			}

			if(!is.null(lose)) mat <- mat[,-lose]
		}
	}
	return(mat)
}
