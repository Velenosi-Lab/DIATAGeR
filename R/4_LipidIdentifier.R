plot_spectra <- function(){
  
  Lipid_Identifed$MSMSref[[1]]
  Merged_Find_all_spectra[j,]$MSMS
  
  Complete <-  Merged_Find_all_spectra[j,]$MSMS[[1]]
  Matched <-  Lipid_Identifed$MSMSref[[1]]
  Matched <- Matched[grep("Fragment|First",Matched$name),]
  textMatched <- -Matched$intensity
  plotMatched <- Matched[!duplicated(Matched$mzref),]
  
  xrange <- c(min(Matched$mz)-25, max(Matched$mz)+25)
  yrange <- max(Complete$intensity[(Complete$mz>xrange[1] & Complete$mz<xrange[2])])*1.05
  
  ggtitleName <- paste0("Feature ID: ", Lipid_Identifed$ID, "\n", Lipid_Identifed$Short.name," ",
                        Lipid_Identifed$Name, "\nReference ID: ", Lipid_Identifed$ID.Number)
  
  plotMatched$intensity <- -plotMatched$intensity
  text_height <-  -yrange/2
  text_shift <- (xrange[2]-xrange[1])/10
  
  label_text <- paste(Matched$name, round(Matched$mzref, 2)) %>% str_replace(., "_mz", "") %>% str_replace(., "_", " ")
  label_text <- paste0(label_text, "\nIntensity ", Matched$intensity)
  
  g <- ggplot() +
    geom_bar(data=Complete,aes(x=mz, y=intensity, fill="Measured"), width=(resolution(Complete$mz)*110+0.3), stat = "identity") +
    geom_bar(data=plotMatched,aes(x=mzref, y=intensity, fill="Reference " ), width=(resolution(Complete$mz)*110+0.3), stat = "identity") +
    geom_text_repel(aes(x = Matched$mzref, y = textMatched, label=label_text),
                    direction="y", nudge_y = text_height, 
                    nudge_x = text_shift,
                    segment.alpha=0.4, size =2.5, segment.size=0.3, hjust = 0) +
    scale_fill_manual("Fill", values=c("#ff3c2e","#5881e0"))+
    # xlim(xrange)+
    xlab("m/z")+
    ylab("Intensity")+
    # ylim(c(-yrange, yrange)) +
    ggtitle(paste0(ggtitleName))+
    theme(legend.direction = "horizontal", 
          legend.position = c(0.5, 0.98),
          legend.background  = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8),
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(family = 'Helvetica', color = '#666666', face = 'bold', size = 11), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 0.45),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line( size=.1, color="black" ),
          panel.grid.major =  element_blank(),
          plot.margin = unit(c(0.2,0.5,0.5,0.5), "cm"))
  
  lipid_name <- str_replace(Lipid_Identifed$Short.name, ":", ".") %>% str_replace(., " ", ".")
  pdf(file=paste0(graphdir,"/", lipid_name, "_", Lipid_Identifed$ID,"_", Lipid_Identifed$ID.Number,"_", ion.mode,"_",version,".pdf"),
      width = 7, height = 6)
  plot(g)
  dev.off()
}


plot_spectra_limitrange <- function(){
  
  Lipid_Identifed$MSMSref[[1]]
  Merged_Find_all_spectra[j,]$MSMS
  
  Complete <-  Merged_Find_all_spectra[j,]$MSMS[[1]]
  Matched <-  Lipid_Identifed$MSMSref[[1]][-grep("Precursor",Lipid_Identifed$MSMSref[[1]]$name),]
  textMatched <- -Matched$intensity
  plotMatched <- Matched[!duplicated(Matched$mzref),]
  
  xrange <- c(min(Matched$mz)-25, max(Matched$mz)+25)
  yrange <- max(Complete$intensity[(Complete$mz>xrange[1] & Complete$mz<xrange[2])])*1.05
  
  ggtitleName <- paste0("Feature ID: ", Lipid_Identifed$ID, "\n", Lipid_Identifed$Short.name," ",
                        Lipid_Identifed$Name, "\nReference ID: ", Lipid_Identifed$ID.Number)
  
  plotMatched$intensity <- -plotMatched$intensity
  text_height <-  -yrange/2
  text_shift <- (xrange[2]-xrange[1])/10
  
  label_text <- paste(Matched$name, round(Matched$mzref, 2)) %>% str_replace(., "_mz", "") %>% str_replace(., "_", " ")
  label_text <- paste0(label_text, "\nIntensity ", Matched$intensity)
  
  g <- ggplot() +
    geom_bar(data=Complete,aes(x=mz, y=intensity, fill="Measured"), width=(resolution(Complete$mz)*110+0.3), stat = "identity") +
    geom_bar(data=plotMatched,aes(x=mzref, y=intensity, fill="Reference " ), width=(resolution(Complete$mz)*110+0.3), stat = "identity") +
    geom_text_repel(aes(x = Matched$mzref, y = textMatched, label=label_text),
                    direction="y", nudge_y = text_height, 
                    nudge_x = text_shift,
                    segment.alpha=0.4, size =2.5, segment.size=0.3, hjust = 0) +
    scale_fill_manual("Fill", values=c("#ff3c2e","#5881e0"))+
    xlim(xrange)+
    xlab("m/z")+
    ylab("Intensity")+
    ylim(c(-yrange, yrange)) +
    ggtitle(paste0(ggtitleName))+
    theme(legend.direction = "horizontal", 
          legend.position = c(0.5, 0.98),
          legend.background  = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8),
          legend.key.size = unit(0.4, "cm"),
          plot.title = element_text(family = 'Helvetica', color = '#666666', face = 'bold', size = 11), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black", size = 0.45),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line( size=.1, color="black" ),
          panel.grid.major =  element_blank(),
          plot.margin = unit(c(0.2,0.5,0.5,0.5), "cm"))
  
  lipid_name <- str_replace(Lipid_Identifed$Short.name, ":", ".") %>% str_replace(., " ", ".")
  pdf(file=paste0(graphdir,"/", lipid_name, "_", Lipid_Identifed$ID,"_", Lipid_Identifed$ID.Number,"_", ion.mode,"_",version,".pdf"),
      width = 6, height = 6)
  plot(g)
  dev.off()
}


#' TG identification 
#'
#' @param DIADataObj The object stores alignment file and spectra
#' @param lipid = "TAG" #DIATAGeR is for triacylglycerols only. 
#' @param ion.mode Pos or Neg to indicate positive or negative ion mode
#' @param format tell the function if the alignment file is MSDIAL or Generic Format
#' @param ppmtol ppm tolerance. Default to 15
#' @param ppmtolMS1 Defaults to ppmtol if not specified. 
#' @param rttol retention time tolerance in seconds. Default to 5.
#' @param intensity.window Two numerical values, specifying the percentage multiplication
#' factor for the precursor intensity that the fragments must fall within. If set 
#' to c(0, 100), the fragment intensity cannot be greater than the precursor intensity. 
#' Defaults to 0 and infinity. 
#' @param which.frag Either 'any' or 'firsttwo', specifying whether any fragment 
#' should be used or if the first two fragments should be used for identification. 
#' Useful when there are three possible fragments in the reference and only two are 
#' required for a positive match. Default 'any'.
#' @param version Appends this character to the output file name. 
#' @param print.MS2spectra If TRUE, prints MS/MS spectra and mirrored reference peaks. Default FALSE. 
#' @param write.annotations If TRUE, writes RDS file of identified lipids. Default TRUE. 
#' @param ms1.precursors Number of precursors found in MS1 to identify compounds. Default to 1 
#' @param ms2.precursors Number of precursor found in MS2 to identify compounds. 
#' "0" is suggested for high collision energy and "1" for low collision energy.
#' @param max.tails list of desired number of carbons and the maximum number of double 
#' bonds for each fatty acyl chain length. 
#' Default = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1, 14.3, 15.3, 16.5, 17.3, 
#' 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0"
#' @param exact.tails library can be customizable for specific tails. (Eg: exact_tails = c("18.1","18.2) creates library of TGs made up 18:1 and 18:2 fatty acids))
#'
#' @return  A folder named "Pos" or "Neg" (depends on the ion mode) will be created and contain a list of identified triacylglycerols 
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' LipidIdentifier(DIADataObj = DIA_Pos,
#' lipid = "TAG", 
#' ion.mode = "Pos", 
#' format = "MSDIAL",
#' ppmtol=15, rttol=5, 
#' intensity.window = c(0,Inf),
#' version = "24", 
#' ms1.precursors= 1, ms2.precursors = 1,
#' max.tails = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1, 14.3, 15.3, 16.5, 17.3, 
#' 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0",
#' exact.tails=NULL,
#' which.frag = "any", 
#' print.MS2spectra = F,
#' write.annotations = T)
#' }

LipidIdentifier <- function(DIADataObj, lipid, ion.mode,
                            format = c("MSDIAL","Generic","Progenesis"), 
                            ppmtol=15, ppmtolMS1, rttol=5,
                            intensity.window = c(0, Inf), 
                            which.frag = "any", ms1.precursors= c(1),
                            ms2.precursors = c(0,1), version,
                            max.tails= "8.2, 9.0, 10.2, 11.0, 12.3, 13.1,14.3, 15.3, 16.5, 17.3, 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0",
                            exact.tails = NULL,
                            print.MS2spectra = FALSE, write.annotations=TRUE){

  ## Calculates mz difference from PPM for MS1
  ppmMS1 <- function(mz, mode, ppmtolMS1){
    if(mode == "-"){
      mz*-ppmtolMS1/1e6+mz
    }else{
      mz*ppmtolMS1/1e6+mz
    }
  }
  
  ## Calculates mz difference from PPM for MS2
  ppm <- function(mz, mode, ppmtol){
    if(mode == "-"){
      mz*-ppmtol/1e6+mz
    }else{
      mz*ppmtol/1e6+mz
    }
  }
  
  ## Unique numbered column names for tibbles
  make_unique = function(x, sep='_'){
    ave(x, x, FUN = function(a){paste(a, 1:length(a), sep=sep)})
  }

  ## MS1 and MS2 parameters the same if not specified
  if(missing(ppmtolMS1)){
    ppmtolMS1 <- ppmtol 
  }
  
  if (format=="Generic") {
    MSDialMinusMSMS <- DIADataObj@Results
    for (i in 1:nrow(MSDialMinusMSMS)){
    sample<-MSDialMinusMSMS$`Spectrum reference file name`[i]
    MSDialMinusMSMS$"Feature_Height"[i]<-MSDialMinusMSMS[i,grep(sample,names(MSDialMinusMSMS))]
    }
  }
  
  if(format=="MSDIAL"){
    MSDialMinusMSMS <- DIADataObj@Results
  } else if(format=="Progenesis"){
    MSDialMinusMSMS <- DIADataObj@AlignResultForProgen
  }
  
  if(format=="MSDIAL"|format=="Progenesis") {
  MSDialMinusMSMS <- MSDialMinusMSMS%>%
    mutate(Feature=gsub(" .*","",`MS1 isotopic spectrum`))%>%
    separate(., col = Feature, sep = ":",into = c("Feature_mz","Feature_Height"), convert = T)
  }
  
  #Load spectra into MSMS_list_Final_Deisotope
  #Spectra are called by list names. 
  if(length(DIADataObj@ResultsNestedSpectra)==0){
    stop("Please Deisotope First")
  }
  
    MSMS_list_Final_Deisotope <- DIADataObj@ResultsNestedSpectra$NestedDeisotopedSpectra
    names(MSMS_list_Final_Deisotope) <- DIADataObj@ResultsNestedSpectra$ID[!is.na(DIADataObj@ResultsNestedSpectra$NestedDeisotopedSpectra)]
  
  #Identify missing spectra
  all_spectra <- as.numeric(names(MSMS_list_Final_Deisotope))
  max_spectra <- max(all_spectra)
  no_spectra <- which(!min(MSDialMinusMSMS$ID):max_spectra %in% all_spectra) 
  no_spectra <- (min(MSDialMinusMSMS$ID):max_spectra)[no_spectra]

  #Import lipid precursor spreadsheet
  # LipidPrecurs <- read_csv(paste0("inst/extdata/",lipid,"_",ion.mode,"_Core.csv")) 
  LipidPrecurs <- getTAGs(tails = max.tails,exact_tails=exact.tails)
  #LipidIDFormat <- read_csv('inst/extdata/LipidIDFormat.csv')
  LipidIDFormat <- list.files(system.file("extdata", package="DIATAGeR"), full.names = T)
  LipidIDFormat <- read_csv(LipidIDFormat[grep("LipidIDFormat", LipidIDFormat)])
  
  NumberOfPrecursors <- as.numeric(LipidIDFormat[grep(paste0("^",lipid,"_",ion.mode),LipidIDFormat$LipidID),2])
  NumberOfFragments <- as.numeric(LipidIDFormat[grep(paste0("^",lipid,"_",ion.mode),LipidIDFormat$LipidID),3])
  

  MSMSPrecursorNumber <- NumberOfPrecursors
  NumberOfPrecursors <- ms1.precursors


  Identified_lipids_Final<-c()
  pb <- txtProgressBar(min = 1, max = dim(LipidPrecurs)[1], style = 3)
  
  #TAG
  # i=61900
  # i=488
  # i=7
  #Plot for poster: i=9292, j=7
  #Plot for BCPN TAG 50:3: LipidPrecurs<-LipidPrecurs[c(31066, 31114, 31537, 42786, 44010),]
  #Plot for BCPN TAG 54:2: LipidPrecurs<-LipidPrecurs[c(31938, 43194, 43232, 53385),]
  #Plot for BCPN TAG 56:3: LipidPrecurs<-LipidPrecurs[c(43246, 43623,44432,44847,53398,54177),]
  #Plot for BCPN TAG 44:1: LipidPrecurs<-LipidPrecurs[c(11167,19453,19871,21886,30614,30670),]
  i=9292
  for (i in 1:nrow(LipidPrecurs)){
    
    setTxtProgressBar(pb, i)
    All_Identified<-c()
    
    sublipid <- LipidPrecurs[i,] %>% .[,which(!is.na(.))]  
    First_Precursor_Known <- sublipid[,"First_Precursor_mz"]%>% as.numeric()
   
    #First precursor matched to feature list
    First_Precursor_Find <- MSDialMinusMSMS[between(MSDialMinusMSMS$mz,
                                                    ppmMS1(First_Precursor_Known,"-",ppmtolMS1),
                                                    ppmMS1(First_Precursor_Known,"+",ppmtolMS1)),]

    if (nrow(First_Precursor_Find)>0) {
      Lipid_Known_all <- sublipid[ c(grep("Precursor", names(sublipid)), 1)]  #Select all precursor m/z and ID 
      Lipid_Find_all <- list()
      
      Lipid_Find_all[[1]] <- First_Precursor_Find
      
      if (NumberOfPrecursors == 1) {
        
        Merged_Find_all <- tibble(Lipid_Find_all[[1]][,1:2], .name_repair = make_unique)
        
      } else { 
        
        for (k in 2:NumberOfPrecursors) {
          #Additional precursors matched from feature list
          Lipid_Find_all[[k]] <- MSDialMinusMSMS[between(MSDialMinusMSMS$mz, 
                                                         ppmMS1(Lipid_Known_all[[k]],"-",ppmtolMS1), 
                                                         ppmMS1(Lipid_Known_all[[k]],"+",ppmtolMS1)),]
        } 
        
        IDs <- list()
        
        if (nrow(Lipid_Find_all[[k]])==0) { 
          next
        } else {
          
          #We want to merge the precursors we have found by retention time. We may 
          #have found more than one feature for each precursor match. So we calculate 
          #whether every combination of precursor 1 and 2 and potentially 3 has a similar retention 
          #time within the rttol range. 
          for (e in 2:length(Lipid_Find_all)){
            temp <- NULL
            
            for (w in 1:nrow(Lipid_Find_all[[1]])){
              temp <- rbind(temp, abs(Lipid_Find_all[[e]]$rt-Lipid_Find_all[[1]]$rt[w]) < rttol/60)
            }
            
            IDs[[e-1]] <- temp
          }
        }
        
        if(e==2){
          #The case of 2 precursors.
          #Returns ID_1 and ID_2 for precursor 1 and 2, respectively 
          temp_1 <- which(IDs[[1]], arr.ind = T)
          Merged_Find_all <- tibble(Lipid_Find_all[[1]][temp_1[,1],][,1:2],
                                    Lipid_Find_all[[2]][temp_1[,2],][,1:2], .name_repair = make_unique)
          
        } else if (e==3) {
          #The case of 3 precursors.
          temp_1 <- which(IDs[[1]], arr.ind = T)
          temp_2 <- which(IDs[[2]], arr.ind = T)
          match_test <- temp_1[,1]==temp_2[,1]
          
          if(any(match_test)){
            Merged_Find_all <- tibble(Lipid_Find_all[[1]][temp_1[,1],][,1:2][match_test,],
                                      Lipid_Find_all[[2]][temp_1[,2],][,1:2][match_test,], 
                                      Lipid_Find_all[[3]][temp_2[,2],][,1:2][match_test,], .name_repair = make_unique)
          }
        }
      }
      
      if (!is.na(Merged_Find_all$ID_1[1])) { 
        
        #Looping through every selected feature 
        for (j in 1:nrow(Merged_Find_all)) {
          
          #Check that spectra is not missing
          if (all(Merged_Find_all$ID_1[j]!=no_spectra) & !Merged_Find_all$ID_1[j]>max_spectra) { 
            
            MSMS_Precursor1 <- MSMS_list_Final_Deisotope[[which(all_spectra==Merged_Find_all$ID_1[j])]]
            
            if (j==1) {
              Merged_Find_all_spectra <- tibble(Merged_Find_all[j,], MSMS=list(tibble(MSMS_Precursor1)))
            } else {
              Merged_Find_all_spectra[j,] <- tibble(Merged_Find_all[j,], MSMS=list(tibble(MSMS_Precursor1)))
            }
            
            Merged_Find_Precursor_spectra <- vector("list", nrow(Merged_Find_all))
            names(Merged_Find_Precursor_spectra) <- Merged_Find_all$ID_1
            
            Precursor_spectra <- NULL
            
            ##Find precursors in the MS/MS
            for (k in 1:MSMSPrecursorNumber) {
              
              Test_frag <- Merged_Find_all_spectra$MSMS[[j]][between(Merged_Find_all_spectra$MSMS[[j]]$mz,
                                                                     ppm(Lipid_Known_all[[k]],"-",ppmtol), 
                                                                     ppm(Lipid_Known_all[[k]],"+",ppmtol)),]
              Test_frag$name <-  names(Lipid_Known_all)[k]
              Test_frag$mzref <- Lipid_Known_all[[k]]
              
              #If more than 1 peak found in the window, the closest is selected
              if (nrow(Test_frag)>1){
                Test_frag <- Test_frag[which.min(abs(Lipid_Known_all[[k]]-Test_frag$mz)),]
              }
              Precursor_spectra <- bind_rows(Precursor_spectra, Test_frag)
            }
            ## 
            
            Merged_Find_Precursor_spectra[[j]] <- Precursor_spectra
            
            #If 1 precursor required in ms2, check if it is first precursor
            if(ms2.precursors==1){
              First_frag <- Merged_Find_Precursor_spectra[[j]][between(Merged_Find_Precursor_spectra[[j]]$mz,
                                                                       ppm(Lipid_Known_all[[1]],"-",ppmtol), 
                                                                       ppm(Lipid_Known_all[[1]],"+",ppmtol)),]
              if (nrow(First_frag)>=1) {
                PrecursorFound <- T
              } else{
                PrecursorFound <- F
              }
            }
            
            #If no precursors required in ms2, move to next.
            if(ms2.precursors==0){
              PrecursorFound <- T
            }
            
            if(PrecursorFound){
              
              ID_list <- Merged_Find_all_spectra[,grep("ID", names(Merged_Find_all_spectra))][j,]  #select by ID
              Lipid_Find_all_edited <- lapply(Lipid_Find_all, function(x) subset(x, x$ID %in% ID_list))
              Fragments_Known <- sublipid[,grep("ragment", names(sublipid))] 
              
              if(ms2.precursors!=0){
                #Determine the number of duplicate lipid tails
                tails <- duplicated(as.character(Fragments_Known[1,]))
                tails <- sum(tails) + 1 #multiplication factor 
                
                #Convert all MS/MS intensities to percentage of adduct (higher precursor)
                MSMS_all_Percent <- Merged_Find_all_spectra$MSMS[[j]] 
                MSMS_all_Percent$intensity <- MSMS_all_Percent$intensity/max(Merged_Find_Precursor_spectra[[j]]$intensity)*100
                
                #filter intensities based on precursor in proportion to tail amount
                Merged_Find_all_spectra$MSMS[[j]] <- Merged_Find_all_spectra$MSMS[[j]][between(MSMS_all_Percent$intensity, 
                                                                                               intensity.window[1],
                                                                                               intensity.window[2]*tails),]
              }
            }else{
              next
            }
            
            Find_MSMS_Fragments <- c()
            
            if (which.frag=="any"){
              for(t in 1:ncol(Fragments_Known)){
                temp <- Merged_Find_all_spectra$MSMS[[j]][abs(Merged_Find_all_spectra$MSMS[[j]]$mz-Fragments_Known[[t]])<Fragments_Known[[t]]*ppmtol/1e6,]  
                temp$name <- names(Fragments_Known)[t]
                temp$mzref <- Fragments_Known[[t]]
                
                
                if(nrow(temp>1)){
                  temp <- temp[which.min(abs(temp$mz-Fragments_Known[[t]])),]
                }
                Find_MSMS_Fragments <- bind_rows(Find_MSMS_Fragments, temp)
              }
              
            } else if(which.frag=="firsttwo"){
              for(t in 1:2){
                temp <- Merged_Find_all_spectra$MSMS[[j]][abs(Merged_Find_all_spectra$MSMS[[j]]$mz-Fragments_Known[[t]])<Fragments_Known[[t]]*ppmtol/1e6,]  
                temp$name <- names(Fragments_Known)[t]
                temp$mzref <- Fragments_Known[[t]]
                
                if(nrow(temp>1)){
                  temp <- temp[which.min(abs(temp$mz-Fragments_Known[[t]])),]
                }
                Find_MSMS_Fragments <- bind_rows(Find_MSMS_Fragments, temp)
              }
            }
            
            #If fragments found, we have a successful match
            if(nrow(Find_MSMS_Fragments) >= NumberOfFragments){
              
              ##Generates nested fragment and precursor names, reference mz, mz, and intensity
              
              matchedPrecursors <- data.frame(name = Precursor_spectra$name,
                                              mzref = Precursor_spectra$mzref,
                                              mz = Precursor_spectra$mz,
                                              intensity = Precursor_spectra$intensity)
              
              matchedFragments <- data.frame(name = Find_MSMS_Fragments$name,
                                             mzref = Find_MSMS_Fragments$mzref,
                                             mz = Find_MSMS_Fragments$mz,
                                             intensity = Find_MSMS_Fragments$intensity)
              
              matchedSpectraKnown <- bind_rows(matchedFragments, matchedPrecursors)
              ##
              
              Lipid_Identifed <- tibble(sublipid[,c(1:6)],
                                        ID.simple=Merged_Find_all_spectra$ID_1[[j]],
                                        ID=paste0(Merged_Find_all_spectra$ID_1[[j]],".",j),
                                        mz= Lipid_Find_all_edited[[1]]$mz,
                                        rt=Lipid_Find_all_edited[[1]]$rt,
                                        FeatureType="Lipid",
                                        FragIntensity=sum(Find_MSMS_Fragments$intensity),
                                        MSMSref=list(matchedSpectraKnown))
              
              ##Identifies potential adducts and fragments in the feature list
              InsourceFragTestDF <- expand.grid(mz=as.numeric(Fragments_Known[1,]),rt=Lipid_Identifed$rt)
              InsourceFrag <- c()
              
              for(t in 1:nrow(InsourceFragTestDF)){
                InsourceFrag <- bind_rows(InsourceFrag, 
                                          MSDialMinusMSMS[(between(MSDialMinusMSMS$mz, 
                                                                   ppm(InsourceFragTestDF$mz[t],"-",ppmtol), 
                                                                   ppm(InsourceFragTestDF$mz[t],"+",ppmtol)) &
                                                           between(MSDialMinusMSMS$rt,
                                                                    InsourceFragTestDF$rt[t]-rttol/60,
                                                                    InsourceFragTestDF$rt[t]+rttol/60)),])
              }
              
              Frags_Identified<-c()
              
              if(nrow(InsourceFrag)>0){
                Frags_Identified <- sublipid[,c(1:6)] %>%
                  bind_cols(ID = c(paste0(InsourceFrag$ID,".",rownames(InsourceFrag)))) %>%
                  bind_cols(InsourceFrag %>% select(c(rt, mz))) %>%
                  mutate(FeatureType = "Fragment", FragIntensity=InsourceFrag$Feature_Height, ID.simple=InsourceFrag$ID ) 
              }
              
              Adduct_Identified <- NULL
              
              if(length(Lipid_Find_all_edited) > 1){
                empty_df_test <- Lipid_Find_all_edited[-1] %>% bind_rows()
                
                if(dim(empty_df_test)[1]>0){
                  Adduct_Identified <- empty_df_test %>%
                    select(c(ID, rt, mz)) %>%
                    mutate( ID.simple=ID , ID = paste0(ID,".", j) ) %>%
                    mutate(FeatureType = "Adduct", FragIntensity = sapply(Lipid_Find_all_edited[-1], "[[", "Feature_Height")) %>% ##SUMMED FRAGINTENSITY
                    cbind(sublipid[,c(1:6)])
                }
              }
              
              ##
              All_Identified <-  bind_rows(All_Identified, Lipid_Identifed, Adduct_Identified, Frags_Identified) 
          
              
              
              if ( print.MS2spectra == TRUE){
                graphdir <- paste0(ion.mode,"/",lipid,"_",version)
                dir.create(ion.mode, showWarnings = F)
                dir.create(graphdir, showWarnings = F) 
                environment(plot_spectra) <- environment()
                suppressWarnings(plot_spectra())
                #environment(plot_spectra_limitrange) <- environment()
                #suppressWarnings(plot_spectra_limitrange())
              } 
            }
          }
        }
      }
      
      if(!is.null(All_Identified)){
        
        #Add all identified to final data frame
        Identified_lipids_Final <- bind_rows(Identified_lipids_Final, All_Identified %>% distinct(ID.simple, .keep_all = T) )
      }
    }
  }
  
  DIADataObj <- addLipidClass(DIADataObj, Lipid=lipid, value=Identified_lipids_Final)
  
  if(write.annotations==T){
    dir.create(ion.mode, showWarnings = F)
    #Gemerate unique names in case of duplicated
    number <- 1
    while (paste0(ion.mode,"/Identified_",lipid,"_",ion.mode,"_",version,"_",number) %in% list.files(ion.mode,full.names = T)){
      number <- number+1
    }
    
    saveRDS(Identified_lipids_Final, file = paste0(ion.mode,"/Identified_",lipid,"_",ion.mode,"_",version,"_",number))
  
    }
  
  invisible(DIADataObj)
}



# LipidIdentifierWrapper<-function(DIADataObj, cores=1, ion.mode, format = c("MSDIAL","Progenesis"), IDParam=NULL,ppmtol=15, rttol=5,
#                                  rt.range = c(0,Inf), intensity.window = c(5,300), version,
#                                  print.MS2spectra = TRUE, write.annotations=TRUE,...){
# 
#   IDFiles <- list.files("inst/extdata/")
#   IDFilesNames<-IDFiles[grep(ion.mode,IDFiles)]
#   IDFilesNames<-sapply(strsplit(IDFilesNames,"_"),"[[",1)
#   pb2 <- txtProgressBar(min = 1, max = length(IDFilesNames), style = 3)
#   
#   if(!is.null(IDParam)){
#     
#     IDParam<-read_csv(IDParam)
#     
#   }
#   
#   for (i in 1:length(IDFilesNames)){
#     print(paste0("Identify ",IDFilesNames[i]," Start"))
#     setTxtProgressBar(pb2,i)
#     suppressMessages(
# 
#       DIADataObj<-LipidIdentifier(DIADataObj, Lipid=IDFilesNames[i], format = format,ion.mode=ion.mode, ppmtol=ppmtol, rttol=rttol,rt.range = rt.range, intensity.window = intensity.window, version=version, print.MS2spectra = print.MS2spectra,
#                                   write.annotations=write.annotations)
#     )
#     print(paste0("Identify ",IDFilesNames[i]," Complete"))
#     
#   }
#   
#   return(DIADataObj)
# }
# 
