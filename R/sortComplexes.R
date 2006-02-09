#function to sort complex estimates into MBME, SBMH, UnRBB

sortComplexes <- function(PCMG,adjMat){

   diag(adjMat) <- 0
   bNames <- rownames(adjMat)
   nComps <- dim(PCMG)[2]
	    

   SBMHs <- which(colSums(PCMG[bNames,])==1)

   UnRBBs <- which(colSums(PCMG)==2 & colSums(PCMG[bNames,])==2)
   keep <- rep(TRUE,length(UnRBBs))
   for (i in UnRBBs){
       tBs <- which(PCMG[,i]==1)
       keep[which(UnRBBs==i)] <- sum(adjMat[tBs,tBs])==1
   }   
   UnRBBs <- UnRBBs[keep]

   MBMEs <- c(1:nComps)[!(1:nComps) %in% c(SBMHs,UnRBBs)]

   if(length(MBMEs)>0){
   MBME <- as.matrix(PCMG[,MBMEs])
   colnames(MBME) <- paste("MBME",1:length(MBMEs),sep="")
   } else MBME <- NA

   if(length(SBMHs)>0){
   SBMH <- as.matrix(PCMG[,SBMHs])
   colnames(SBMH) <- paste("SBMH",1:length(SBMHs),sep="")
   } else SBMH <- NA

   if(length(UnRBBs)>0){
   UnRBB <- as.matrix(PCMG[,UnRBBs])
   colnames(UnRBB) <- paste("UnRBB",1:length(UnRBBs),sep="")
   } else UnRBB <- NA


   sComps <- list(MBME=MBME,SBMH=SBMH,UnRBB=UnRBB)

   return(sComps)

}