#function to sort complex estimates into MBME, SBMH, UnRBB

sortComplexes <- function(PCMG,adjMat){

   diag(adjMat) <- 0
   bNames <- rownames(adjMat)
   nComps <- length(PCMG)

   nBFUN <- function(x) sum(x %in% bNames)
   nBs <- unlist(lapply(PCMG,FUN=nBFUN))
   nT <- unlist(lapply(PCMG,FUN=length))
	    
   SBMHi <- which(nBs==1)
   SBMH <- PCMG[SBMHi]

   UnRBBi <- which(nBs==2 & nT==2)
   keep <- rep(TRUE,length(UnRBBi))
   for (i in 1:length(UnRBBi)){
       tBs <- PCMG[[UnRBBi[i]]]
       keep[i] <- sum(adjMat[tBs,tBs])==1
   }   
   UnRBBi <- UnRBBi[keep]
   UnRBB <- PCMG[UnRBBi]

   MBMEi <- c(1:nComps)[!(1:nComps) %in% c(SBMHi,UnRBBi)]
   MBME <- PCMG[MBMEi]

   if(length(MBMEi)>0){
   names(MBME) <- paste("MBME",1:length(MBMEi),sep="")
   } else MBME <- NA

   if(length(SBMHi)>0){
   names(SBMH) <- paste("SBMH",1:length(SBMHi),sep="")
   } else SBMH <- NA

   if(length(UnRBBi)>0){
   names(UnRBB) <- paste("UnRBB",1:length(UnRBBi),sep="")
   } else UnRBB <- NA


   sComps <- list(MBME=MBME,SBMH=SBMH,UnRBB=UnRBB)

   return(sComps)

}