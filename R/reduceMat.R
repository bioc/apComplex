#function to eliminate columns that are subsets or identically equal 
#to other columns

#if working in identical mode, will keep the column that appears first 

reduceMat <- function(mat,compare="equal"){

	#first order columns from largest column sum to smallest column sum
	mat <- mat[,order(colSums(mat),decreasing=FALSE)]

	#take off colnames of mat
	colnames(mat) <- NULL

	if(compare=="equal"){

		i <- 1
		nCols <- dim(mat)[2]
		
		continue1 <- i<nCols

		while(continue1){
	
		Vec <- mat[,i]
		lose <- whichVecInMat(Vec,mat)
		if(length(lose)>1){
			lose <- lose[-1]
			mat <- mat[,-lose]
			nCols <- dim(mat)[2]
		}
		i <- i+1
		continue1 <- i <= nCols
		}	
	matfull <- mat
	}


	if(compare=="less"){

		i <- 1
		nCols <- dim(mat)[2]

		#find columns with singly affiliated members
		temp <- which(rowSums(mat)==1)

		matfull <- mat
		if(length(temp)>0){			

			if(length(temp)==1){
			notestset <- which(mat[temp,]==1)
			}
			if(length(temp)>1){
			notestset <- which(colSums(mat[temp,])>0)
			}

			mat <- mat[-temp,]
		}
				

		continue1 <- i<nCols
		lose <- NULL

		while(continue1){

		Vec <- mat[,i]
		test <- vecInMat(Vec,as.matrix(mat[,(i+1):nCols]),
							compare="less")
		if(test) lose <- c(lose,i) 
		i <- i+1		
		continue1 <- i<nCols
		#print(paste("i =",i))
		}

		if(length(temp)>0) lose <- lose[!(lose %in% notestset)]
		mat <- mat[,-lose]

		
	
	}	

	matfull <- matfull[,-lose]

	return(matfull)
}
