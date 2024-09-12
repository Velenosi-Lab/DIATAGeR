#' Generating a list of TGs passed through a pre-determined FDR 
#'
#' @param identified.lipids List of combined TG annotations 
#' @param cutoff Empirical FDR. Default 0.1
#' @param standard.file list of standard TGs (Optional)
#' @param iteration #Default 10
#'
#' @return list of TGs that pass through the FDR cutoff
#' @export
#'
#' @examples
#'  \dontrun{
#' ResultsTG_FDR <- FDRCutoff(identified.lipids=ResultsTG,
#' standard.file = "Standards.csv",
#' iteration=10,
#' cutoff = 0.1)
#' }

FDRCutoff<- function(identified.lipids, 
                     standard.file=NULL, 
                     iteration=3, 
                     cutoff = 0.1) {
  
  FDRCutoff_file<- function(identified.lipids, cutoff = 0.1){
    
    #List of scores to try empirically 
    scores <- round(seq(0.7, 1,  length.out=6500), 5)
    
    #Score for empirical FDR to be below
    FDRcut <- cutoff+0.001
    
    #Function tries scores as cutoff, to find empirical FDRs
    getFDR <- function(){
      for (i in 1:length(scores)){
        subset_temp <- identified.lipids[identified.lipids$Score>scores[i],]
        FDR <- sum(subset_temp$Class=="Fp")/sum(subset_temp$Class=="Tp")
        
        if(is.na(FDR)){
          subset_temp <- identified.lipids[identified.lipids$Score>scores[i-1],]
          FDR <- sum(subset_temp$Class=="Fp")/sum(subset_temp$Class=="Tp")
          Cutoff <- scores[i-1]
          
          print(paste("Empirical FDR:", round(FDR, 4)))
          assign("CutoffFDR", FDR, envir = parent.frame())
          
          return(Cutoff)
          break
        }
        
        if(FDR < FDRcut){
          
          Cutoff <- scores[i]
          
          print(paste("Empirical FDR:", round(FDR, 4)))
          assign("CutoffFDR", FDR, envir = parent.frame())
          
          return(Cutoff)
          break
        }
      }
    }
    
    #It may be possible for the score to be 0.7, so we try with a lower range 
    test <- try(Cutoffx <- getFDR(), silent = F)
    
    if(inherits(test, what = "try-error")){
      scores <- round(seq(0.4, 0.8,  length.out=8000), 5)
      Cutoffx <- getFDR()
    }
    
    # tryCatch({
    #   Cutoffx <- getFDR()
    # },
    # error = function(e) {
    #   scores <- round(seq(0.4, 0.8,  length.out=8000), 5)
    #   Cutoffx <- getFDR()
    # })
    
    #If the empirical cutoff is below the desired cutoff, we keep increasing the cutoff 
    #until it is above. 
    while(CutoffFDR < cutoff){
      FDRcut <- FDRcut+0.001
      test <- try(Cutoffx <- getFDR(), silent = F)
      if(inherits(test, what = "try-error")){
        scores <- round(seq(0.4, 1,  length.out=8000), 5)
        Cutoffx <- getFDR()
      }
      if(FDRcut>cutoff+0.02){
        break
      }
    }
    
    ret <- identified.lipids[identified.lipids$Score>Cutoffx,]
    return(ret)
  }
  
  
  std <- c()
  srt <- c()
  lng <- c()
  auc <- c()
  fdr <- c()
  FDRout<-c()
  p=1
  for (p in 1:iteration) {
    identified.lipids <- getScoreNew(identified.lipids, standard.file= standard.file,
                                     rt.type=c("truncated"), rttol=2)
    
    identified.lipids.FDR <- FDRCutoff_file(identified.lipids = identified.lipids,
                                            cutoff = cutoff)
    
    auc <- c(auc, AUC)
    srt <- c(srt, unique(identified.lipids.FDR$Short.name) %>% length() )
    lng <- c(lng, unique(identified.lipids.FDR$Name) %>% length() )
    FDRout[[p]]<-identified.lipids.FDR
    
    if (!is.null(standard.file)) {
      fdr <- c(fdr,standardsFDR(identified.lipids, standard.file= standard.file, rt.type=c("truncated"), rttol=2))
      sds <- standardsCheck(identified.lipids.FDR, standard.file= standard.file, rt.type=c("truncated"), rttol=2)
      std <- c(std, nrow(sds))
    }
  }
  
  #Average value of AUC, species and molecular species found
  cat("Results after ",iteration,"iteration(s):","\n")
  cat("The average number of species found: ",mean(srt),"\n")
  cat("The average number of molecular species found: ",mean(lng),"\n")
  cat("The average AUC found: ",mean(auc),"\n")
  
  if (!is.null(standard.file)){
    standard.file<-read.csv(standard.file)
    cat("\n")
    cat("For standard lipids:","\n")
    cat("Number of standards found after FDR cutoff: ",mean(std),"/",nrow(standard.file), "\n")
    cat("The FDR required to find all standards: ",round(mean(fdr), 4),"\n")
  }
  else{
    cat("Standard file was not provided. FDR required to find all standards and number of standards 
      found were not calculated.")
  }
  
  cat("The output is the file with highest number of molecular species.","\n")
  
  FDRout_1<-FDRout[[which.max(lng)]]
  
  return(FDRout_1)
  
}

