
#A function to decide whether or not to combine complexes based on changes
# in the product of L and C

#currently this function only includes the Fisher's exact portion for complexes with less than 20 proteins. For larger complexes, the binomial criteria is adequate.

LCdelta <- function(comp1,comp2,cMat,dataMat,baitList,simMat,mu,alpha,beta,wsVal=20000000){

N <- dim(dataMat)[1]
M <- dim(dataMat)[2]-N


#compute likelihood contribution for separate complexes

#first for comp1
tP1 <- names(which(cMat[,comp1]==1))   
tB1 <- tP1[which(tP1 %in% baitList)]

nB1 <- length(tB1)
nH1 <- length(tP1) - nB1

temp1 <- matrix(dataMat[tB1,tP1],ncol=(nB1+nH1))

sim1 <- matrix(simMat[tB1,tP1],ncol=(nB1+nH1))
adjBin1 <- sum(temp1*(mu+alpha+beta*sim1)-log(1+exp(mu+alpha+beta*sim1)))-
		nB1*(mu+alpha+beta-log(1+exp(mu+alpha+beta)))

X1 <- sum(temp1)-nB1

b11 <- colSums(matrix(temp1[,1:nB1],ncol=nB1))
b11[1:nB1] <- b11-1
b01 <- (nB1-1)-b11
h11 <- NULL
h01 <- NULL

if(nH1>0){
	h11 <- colSums(matrix(temp1[,(nB1+1):(nB1+nH1)],ncol=nH1))
	h01 <- nB1 - h11
}

bh01 <- c(b01,h01)
bh11 <- c(b11,h11)



#repeat for comp2
tP2 <- names(which(cMat[,comp2]==1))
tB2 <- tP2[which(tP2 %in% baitList)]

nB2 <- length(tB2)
nH2 <- length(tP2) - nB2

temp2 <- matrix(dataMat[tB2,tP2],ncol=(nB2+nH2))

sim2 <- matrix(simMat[tB2,tP2],ncol=(nB2+nH2))
adjBin2 <- sum(temp2*(mu+alpha+beta*sim2)-log(1+exp(mu+alpha+beta*sim2)))-
		nB2*(mu+alpha+beta-log(1+exp(mu+alpha+beta)))



X2 <- sum(temp2)-nB2

b12 <- colSums(matrix(temp2[,1:nB2],ncol=nB2))
b12[1:nB2] <- b12-1
b02 <- (nB2-1)-b12
h12 <- NULL
h02 <- NULL

if(nH2>0){
	h12 <- colSums(matrix(temp2[,(nB2+1):(nB2+nH2)],ncol=nH2))
	h02 <- nB2 - h12
}

bh02 <- c(b02,h02)
bh12 <- c(b12,h12)



cX1 <- lgamma(nB1*(nB1+nH1-1)+1) - lgamma(X1+1) - lgamma(nB1*(nB1+nH1-1)-X1+1)
cX2 <- lgamma(nB2*(nB2+nH2-1)+1) - lgamma(X2+1) - lgamma(nB2*(nB2+nH2-1)-X2+1)


#compute likelihood contribution if left as 2 complexes

lK <- cX1+cX2+adjBin1+adjBin2


#combine complexes and compute likelihood contribution

tP <- names(which(rowSums(cMat[,c(comp1,comp2)])>0))
tB <- tP[which(tP %in% baitList)]

nB <- length(tB)
nH <- length(tP) - nB

temp <- matrix(dataMat[tB,tP],ncol=(nB+nH))

sim <- matrix(simMat[tB,tP],ncol=(nB+nH))
adjBin <- sum(temp*(mu+alpha+beta*sim)-log(1+exp(mu+alpha+beta*sim)))-
		nB*(mu+alpha+beta-log(1+exp(mu+alpha+beta)))



X <- sum(temp)-nB

b1 <- colSums(matrix(temp[,1:nB],ncol=nB))
b1[1:nB] <- b1-1
b0 <- (nB-1)-b1
h1 <- NULL
h0 <- NULL

if(nH>0){
	h1 <- colSums(matrix(temp[,(nB+1):(nB+nH)],ncol=nH))
	h0 <- nB - h1
}


bh0 <- c(b0,h0)
bh1 <- c(b1,h1)


cX <- lgamma(nB*(nB+nH-1)+1) - lgamma(X+1) - lgamma(nB*(nB+nH-1)-X+1)

#compute likelihood contribution if combined into 1 complex
lKm1 <- cX+adjBin

#find adjustments for edges that used to be in "no complex" bin
#make sure dataMat and cMat have their rows and columns in the same order

theseComps <- which(colSums(cMat[tP,])>0)

exCadjMat <- (1*(cMat[tP,theseComps] %*% t(cMat[tP,theseComps])) >0)[tB,]
diag(exCadjMat) <- 1


part4 <-  sum((matrix(1,nB,(nB+nH))-exCadjMat)*(dataMat[tB,tP]*alpha-
	log(1+exp(mu+alpha+beta*simMat[tB,tP])) + 
	log(1+exp(mu+beta*simMat[tB,tP]))))


ans <- lKm1-lK+part4

#if total size of complex is less than 20 proteins, 
#include Fisher's exact component
if((nB+nH)<20){

	ans <- ans + log(fisher.test(rbind(bh0,bh1),workspace=wsVal)$p.value) - 
		log(fisher.test(rbind(bh01,bh11),workspace=wsVal)$p.value) -
		log(fisher.test(rbind(bh02,bh12),workspace=wsVal)$p.value)

}
return(ans)
}


 
