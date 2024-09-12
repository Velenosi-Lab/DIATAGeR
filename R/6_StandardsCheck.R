####Check if standards exist 
#rt.type is the format of the rt in the standards file, whether it's an exact rt estimate, or a truncated rt time. 
#In the case exact rt estimate, a rttol range (+/-) will be used. In the case of a truncated rt time, the number will
#be searched ignoring the next digit (eg. 15.8 will find lipids from 15.800 to 15.899).


standardsCheck <- function(TAGsAll, standard.file= "NewData/MSe_Pos_Data/Standards.csv", rt.type=c("truncated", "exact"), rttol=2){
  ret <- c()
  standards <- read.csv(standard.file)
  Lipids <- TAGsAll[TAGsAll$Name %in% standards$name,]
  Lipids <- Lipids[Lipids$Class=="Tp",]
  
  
  if (rt.type=="truncated") {
    standards$rt <- format(round(  standards$rt , 1), nsmall = 1)
    
    for (q in 1:nrow(standards)){
      temp <-  Lipids[grep(standards$rt[q], Lipids$rt)[(grep(standards$rt[q], Lipids$rt) %in% which(Lipids$Name %in% standards$name[q]))],]
      if (nrow(temp) > 1) {
        temp <- temp[order(temp$FragIntensity, decreasing = T),][1,]
      }
      ret <- rbind(ret,temp) 
    }
  }
  
  if (rt.type=="exact") {
    for (q in 1:nrow(standards)){
      temp <-  Lipids[grep(standards$name[q], Lipids$Name),]
      temp <- temp[(temp$rt > standards$rt[q]-rttol/60 & temp$rt < standards$rt[q]+rttol/60),]
      if (nrow(temp) > 1) {
        temp <- temp[order(temp$FragIntensity, decreasing = T),][1,]
      }
      ret <- rbind(ret,temp)
    }
  }
  ret
}


####Check if standards exist 
#rt.type is the format of the rt in the standards file, whether it's an exact rt estimate, or a truncated rt time. 
#In the case exact rt estimate, a rttol range (+/-) will be used. In the case of a truncated rt time, the number will
#be searched ignoring the next digit (eg. 15.8 will find lipids from 15.800 to 15.899).

standardsFDR <- function(TAGsAll, standard.file= "NewData/MSe_Pos_Data/Standards.csv", rt.type=c("truncated", "exact"), rttol=2){
  stds <- standardsCheck(TAGsAll, standard.file= standard.file, rt.type=rt.type, rttol=rttol)
  min_score <- stds[order(stds$Score),][1,]
  subset_temp <- TAGsAll[TAGsAll$Score>=min_score$Score,]
  FDR <- sum(subset_temp$Class=="Fp")/sum(subset_temp$Class=="Tp")
  FDR
  print(round(FDR, 4))
}








# Standards <- data.frame(name=c( "TAG 16:0-16:0-18:0", "TAG 16:0-16:0-18:1",
#                                  "TAG 18:1-18:1-18:1", "TAG 16:0-16:0-16:0",
#                                  "TAG 16:0-18:1-18:1", "TAG 18:2-18:2-18:2", 
#                                  "TAG 18:0-18:0-18:1", "TAG 16:0-18:1-18:2", 
#                                  "TAG 16:1-16:1-16:1"),
#                         rt=c(18.3, 17.9, 18.0, 17.9, 17.9, 17.1, 18.6, 17.6, 16.8)
# )
# write.csv(Standards, "NewData/MSe_Pos_Data/Standards.csv",row.names = F)


