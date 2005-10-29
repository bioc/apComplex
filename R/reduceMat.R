#function to eliminate columns that are less than or identically equal 
#to other columns

#if working in identical mode, will keep the column that appears first 

reduceMat <- function(mat,compare="equal"){


	#take off colnames of mat
	colnames(mat) <- NULL
	
	if(compare=="equal"){
	
		cS <- colSums(mat)
		cSt <- table(cS)
		cStt <- as.numeric(names(cSt[cSt>1]))

		#find columns with singly affiliated members
		#these columns will not be tested or removed

		temp <- which(rowSums(mat)==1)

		if(length(temp)>0){			

			if(length(temp)==1){
			notest <- which(mat[temp,]>0)
			}
			if(length(temp)>1){
			notest <- which(colSums(mat[temp,,drop=FALSE])>0)
			}
		} else {notest <-  NULL}
				
	
		keep <- rep(TRUE,length(cS))
		if(length(notest)>0) keep[notest] <- TRUE

		for (j in cStt){

			testset <- which(cS==j)

			for (i in testset){
			k <- which(testset==i)
			if(keep[i]){
			
			l <- whichVecInMat(mat[,i],
				as.matrix(mat[,testset[-k]]))
 			if(length(l)>0){
				keep[testset[-k][l]] <- FALSE
			}
			}
			}
		}	

	matRed <- mat[,keep,drop=FALSE]

	}

	if(compare=="less"){

		#order columns from largest to smallest column sum

		ord <- order(colSums(mat),decreasing=TRUE)
		matOrd <- mat[,ord,drop=FALSE]

		nCols <- dim(mat)[2]

		#find columns with singly affiliated members
		#these columns will not be tested or removed
		temp <- which(rowSums(matOrd)==1)

		if(length(temp)>0){			

			if(length(temp)==1){
			testset <- which(matOrd[temp,,drop=FALSE]==0)
			}
			if(length(temp)>1){
			testset <- which(colSums(matOrd[temp,,drop=FALSE])==0)
			}
		} else {testset <- 1:nCols}
				
		#don't directly test first column
		if(1 %in% testset) testset <- testset[-which(testset==1)]
		lose <- NULL

		for (i in testset){

		Vec <- matOrd[,i,drop=FALSE]
		test <- vecInMat(Vec,as.matrix(matOrd[,1:(i-1),drop=FALSE]),
							compare="less")
		if(test) lose <- c(lose,i) 

		}		
	
	if(length(lose)>0){
		matRed <- mat[,-ord[lose],drop=FALSE]
	}else matRed <- mat
	}	


	return(matRed)
}
