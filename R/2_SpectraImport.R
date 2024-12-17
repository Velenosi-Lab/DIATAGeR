#' MS2 spectra import
#'
#' @param DIADataObj object stores the feature list
#' @param fileORfoldername path to the folder that stores centroid MS2 spectra
#' @param ion.mode Pos or Neg to indicate positive or negative ion mode
#' @param results.file.type tell the function if the feature list is MSDIAL alignment file or Generic Format
#' @param spectra.type Only for centroid data. The user will have determined the spectra type
#' when exporting from MSDIAL so this is more of a reminder and will ensure proper functions are called.
#' @param spectra.file.type Specify txt, msp or mgf files. The user will have determined the spectra name
#' when exporting from MSDIAL so this is more of a reminder and will ensure proper functions are called.
#' @param rttol retention time tolerance in seconds.
#' @param ppmtol ppm tolerance. Default to 5. This parameter is for for generic format since it retrieve the spectra based on rt and m/z similarity
#'
#' @return MSMS spectra of feature list
#' @export
#'
#' @examples 
#' \dontrun{
#' DIA_Pos <- SpectraImport(DIADataObj = DIA_Pos, 
#' fileORfoldername =  "Data/centroid",
#' ion.mode = "Pos", 
#' spectra.file.type ="txt",
#' results.file.type = "MSDIAL",
#' rttol = 10)
#' }

SpectraImport <-
  function(fileORfoldername,
           DIADataObj,
           ion.mode,
           results.file.type = c("MSDIAL", "Generic","Progenesis"),
           spectra.type = "centroid",
           spectra.file.type = c("txt","msp","mgf","mgf_mzmine"),
           ppmtol = 5,
           rttol=2) {
    
    dir.create("Troubleshooting", showWarnings = F)
    
    ppm <- function(mz, mode, ppmtol){
      if(mode == "-"){
        mz*-ppmtol/1e6+mz
      }else{
        mz*ppmtol/1e6+mz
      }
    }
    
    parse_msp <- function(file_path) {
     gc()
      msp <- scan(file = file_path, what = "",
                  sep = "\n", quote = "",
                  allowEscapes = FALSE,
                  quiet = TRUE)
      # general indices
      indices <- grep("^SCANNUMBER:", msp)
      
      # rt
      rt_indices <- indices +1
      rt_temp <- as.numeric(sub("RETENTIONTIME: ", "", msp[rt_indices]))
      
      # mz
      mz_indices <- indices -2
      mz_temp <-as.numeric(sub("PRECURSORMZ: ", "", msp[mz_indices]))
      
      #msms
      start_indices <- indices + 5
      end_indices <- c(tail(indices, -1) - 5, length(msp))
      fragments_temp <- sapply(seq_along(start_indices), function(i) {
        paste(msp[start_indices[i]:end_indices[i]], collapse = " ")
      })
      
      
      spectra_file <- data.frame(rt = rt_temp, mz = mz_temp,msms=fragments_temp)
      
      return(spectra_file)
    }
    parse_mgf <- function(file_path) {
      gc()
      
      mgf <- readLines(file_path,warn=FALSE)
      # general indices
      indices <- grep("^SCANS=", mgf)
      
      # rt
      rt_indices <- indices +1
      rt_temp <- as.numeric(sub("RTINMINUTES=", "", mgf[rt_indices]))
      
      # mz
      mz_indices <- indices + 2
      mz_temp <-as.numeric(sub("PEPMASS=", "", mgf[mz_indices]))
      
      #msms
      start_indices <- indices + 5
      last_ion <- length(mgf)-2
      end_indices <- c(tail(indices, -1) - 5, last_ion)
      fragments_temp <- sapply(seq_along(start_indices), function(i) {
        paste(mgf[start_indices[i]:end_indices[i]], collapse = " ")
      })
      
      
      spectra_file <- data.frame(rt = rt_temp, mz = mz_temp,msms=fragments_temp)
      
      return(spectra_file)
    }
    
    parse_mgf_mzmine <- function(file_path) {
      gc()
      
      mgf <- readLines(file_path,warn=FALSE)
      # general indices
      indices <- grep("Title: ", mgf)
      
      # rt
      rt_temp <- as.numeric(sub(".*RT: ([0-9.]+) min.*", "\\1", mgf[indices]))
      
      #msms
      start_indices <- indices + 1
      last_ion <- length(mgf)-2
      end_indices <- c(tail(indices, -1) - 7, last_ion)
      fragments_temp <- sapply(seq_along(start_indices), function(i) {
        paste(mgf[start_indices[i]:end_indices[i]], collapse = " ")
      })
      
      
      spectra_file <- data.frame(rt = rt_temp,msms=fragments_temp)
      return(spectra_file)
    }
    
    if (length(spectra.type) > 1)
      stop("Spectratype must be centroid")

    ##START IF CENTROID
    if (spectra.type == "centroid") {
      
      #if (spectra.file.type == "msp")
        #stop("spectra.file.type must be txt")
      if (dir.exists(fileORfoldername) == F)
        stop("PeakResults Folder does not exist")
      
      #Extract the results from DIADataObj
      if (results.file.type == "Progenesis") {
        featureDF <- DIADataObj@AlignResultForProgen %>%
          select("ID", "rt", "mz", "Spectrum reference file name", "MS/MS spectrum")
        featureDFFinal <- DIADataObj@AlignResultForProgen
      }  
      if (results.file.type == "MSDIAL") {
        featureDF <- DIADataObj@Results %>%
          select("ID", "rt", "mz", "Spectrum reference file name", "MS/MS spectrum") %>% 
          rename(msms= `MS/MS spectrum`)
        featureDFFinal <- DIADataObj@Results 
      }
      if(results.file.type == "Generic") {
        featureDF <- DIADataObj@Results %>%
          select("ID", "rt", "mz", "Spectrum reference file name")
        featureDFFinal <- DIADataObj@Results
      }

      AllMatch <- featureDF
      AllMatch$fullmsms <- NA
      
      if (spectra.file.type=="msp"|spectra.file.type=="mgf"|spectra.file.type=="txt"){
        peakListFiles <-list.files(fileORfoldername,pattern = paste0(".", spectra.file.type),full.names = T)  
      }
      
      if (spectra.file.type=="mgf_mzmine"){
        peakListFiles <-list.files(fileORfoldername,pattern = paste0(".", "mgf"),full.names = T)  
      }
      
      lst <- AllMatch %>% group_by(`Spectrum reference file name`) %>% group_map(~.x,.keep = T)
      names(lst) <- sapply(lst, function(x) x$`Spectrum reference file name`[[1]])

      pb <- txtProgressBar(min = 1, max = length(lst) , style = 3)
p=4
      for (p in 1:length(lst)){
        sampl_num <- lst[[p]]
        if (spectra.file.type=="msp"){
          spectra_file <- parse_msp(peakListFiles[grep(names(lst[p]), peakListFiles)])
        }
        if (spectra.file.type=="mgf"){
          spectra_file <- parse_mgf(peakListFiles[grep(names(lst[p]), peakListFiles)])
        }
        if (spectra.file.type=="txt"){
          spectra_file <- read_delim(peakListFiles[grep(names(lst[p]), peakListFiles)], delim = "\t", col_types = cols(),progress = F, show_col_types = F) %>%  
            select("RT (min)", "Precursor m/z", "MSMS spectrum") %>%
            setNames(c("rt", "mz", "msms"))
        }
        if (spectra.file.type=="mgf_mzmine"){
          spectra_file <- parse_mgf_mzmine(peakListFiles[grep(names(lst[p]), peakListFiles)])
        }
        
        ## Get MSMS spectra based on 'mz' in alignment file
          ## Generic feature list with spectra exported from msdial 
        if (results.file.type == "Generic" && (spectra.file.type == "txt" || spectra.file.type == "msp" || spectra.file.type== "mgf")) {
          k=1
          for (k in 1:nrow(sampl_num)){
            
            setTxtProgressBar(pb, p+k/nrow(sampl_num))
            
            ppmtolMS1<-ppmtol
            rttolMS1<-rttol
            
            match <- spectra_file[
              (between(spectra_file$mz,
                       ppm(sampl_num$mz[k], "-", ppmtolMS1),
                       ppm(sampl_num$mz[k], "+", ppmtolMS1))) &
                (between(spectra_file$rt,
                         sampl_num$rt[k] - rttolMS1 / 60,
                         sampl_num$rt[k] + rttolMS1 / 60))
              ,]
            
            if (nrow(match)==1) {
              lst[[p]]$fullmsms[k] <- match$msms
            }
            
            if (nrow(match)>1) {
              lst[[p]]$fullmsms[k] <- match[which.min(abs(sampl_num$mz[k]-match$mz)),]$msms
            }
            
            if (nrow(match)==0){
              
              ppmtolr<- ppmtolMS1
              
              repeat{
                ppmtolr<-ppmtolr + 15
                matchtemp <- spectra_file[
                  (between(spectra_file$mz,
                           ppm(sampl_num$mz[k], "-", ppmtolr),
                           ppm(sampl_num$mz[k], "+", ppmtolr))) &
                    (between(spectra_file$rt,
                             sampl_num$rt[k] - rttolMS1 / 60,
                             sampl_num$rt[k] + rttolMS1 / 60))
                  ,]
                
                if (nrow(matchtemp)>=1 |ppmtolr > 1000) {
                  match <- matchtemp
                  break
                }
              }
                
                if (nrow(match) > 1) {
                  lst[[p]]$fullmsms[k] <- match[which.min(abs(sampl_num$mz[k] - match$mz)),]$msms
                } 
                if (nrow(match) == 1) { 
                  lst[[p]]$fullmsms[k] <- match$msms
                }
                if (nrow(match)==0) {
                  lst[[p]]$fullmsms[k] <- NA
                }
            
          }
        }
        }
        
          ## Generic with spectra exported from mzmine 
        if (results.file.type == "Generic" && spectra.file.type == "mgf_mzmine") {
          k=1
          for (k in 1:nrow(sampl_num)){
            
            setTxtProgressBar(pb, p+k/nrow(sampl_num))
            
            rttolMS1<-rttol
            
            match <- spectra_file[
                (between(spectra_file$rt,
                         sampl_num$rt[k] - rttolMS1 / 60,
                         sampl_num$rt[k] + rttolMS1 / 60))
              ,]
            
            if (nrow(match)==1) {
              lst[[p]]$fullmsms[k] <- match$msms
            }
            
            if (nrow(match)>1) {
              lst[[p]]$fullmsms[k] <- match[which.min(abs(sampl_num$rt[k]-match$rt)),]$msms
            }
            
            if (nrow(match)==0){
              
              rttolr<- rttolMS1
              
              repeat{
                rttolr<-rttolr + 2
                matchtemp <- spectra_file[
                    (between(spectra_file$rt,
                             sampl_num$rt[k] - rttolr / 60,
                             sampl_num$rt[k] + rttolr / 60))
                  ,]
                
                if (nrow(matchtemp)>=1 |rttolr > 100) {
                  match <- matchtemp
                  break
                }
              }
              
              if (nrow(match) > 1) {
                lst[[p]]$fullmsms[k] <- match[which.min(abs(sampl_num$rt[k]-match$rt)),]$msms
              } 
              if (nrow(match) == 1) { 
                lst[[p]]$fullmsms[k] <- match$msms
              }
              if (nrow(match)==0) {
                lst[[p]]$fullmsms[k] <- NA
              }
              
            }
          }
        }
        
        ## MSDIAL feature list with spectra exported from mzmine 
        if (results.file.type == "MSDIAL") {
          k=1
          for (k in 1:nrow(sampl_num)){
            setTxtProgressBar(pb, p+k/nrow(sampl_num))
            
            if(is.na(sampl_num$msms[[k]])==T){ 
              next
            }
            
            rt_dist <- abs(sampl_num$rt[k] - spectra_file$rt)
            temp_subset <- spectra_file[rt_dist < (rttol)/60,]
            
            tofind <- sampl_num$msms[[k]] %>% str_trunc(., width=100, ellipsis = "") %>% str_split(., ":|\\ ") %>% .[[1]]
            tofind <- paste0(tofind[grep("\\.", tofind)], ".*", collapse = "")
            
            subset_msms <- temp_subset$msms %>% str_trunc(., width=200, ellipsis = "")
            select_it <- grep(tofind, subset_msms)
            if(length(select_it)==0){
              subset_msms <- temp_subset$msms %>% str_trunc(., width=3000, ellipsis = "")
              select_it <- grep(tofind, subset_msms)
              if(length(select_it)==0){
                subset_msms <- temp_subset$msms
                select_it <- grep(tofind, subset_msms)
                if(length(select_it)==0){
                  next
                }
              }
            }
            
            lst[[p]]$fullmsms[k] <- temp_subset[select_it[1],]$msms
            
          }
        }
      }

      
      AllMatch <- do.call(rbind,lst) 
      AllMatch <- AllMatch[order(AllMatch$ID),]
      #check for missing msms
      MissingSpectra <- AllMatch %>%
        filter(is.na(fullmsms))
      
      if (dim(MissingSpectra)[1] > 0) {
        write.csv(
          MissingSpectra,
          paste0(
            "Troubleshooting/",
            ion.mode,
            "_AlignmentFile_MissingSpectra.csv"
          )
        )
        
        print(
          paste0(
            "Missing ",
            nrow(MissingSpectra),
            "/",
            nrow(featureDF),
            " spectra, increase ppmtol and/or rttol, see ",
            ion.mode,
            "_AlignmentFile_MissingSpectra.csv"
          )
        )
      }
      
      # featureDFFinal <-
      #   featureDFFinal[featureDFFinal$ID %in% AllMatch$ID,]
      
      # DIADataObj <-
      #   addAlignResultForProgen(DIADataObj, featureDFFinal)
      
      ##MSMS_Spectra original 
      subdf1 <- select(AllMatch, c(ID, fullmsms))
      split_list <- split(as.data.frame(subdf1), subdf1$ID)
      split_list <- lapply(split_list, function(x) x[,-1, drop = T])
      
      splits <- function(x){
        sep_peaks <- unlist(strsplit(x, "\\s+"))
        sep_heights <- unlist(strsplit(sep_peaks, ":")) %>% as.numeric(.) %>% matrix(.,ncol=2,byrow = T) %>% 
          as_tibble(sep_heights,.name_repair = "minimal")
        names(sep_heights) <- c("mz", "intensity")
        return(sep_heights)
      }
      
      splits_msp_mgf <- function(x){
        sep_peaks <- unlist(strsplit(x, "\\s+")) %>% as.numeric(.) %>% matrix(.,ncol=2,byrow = T) %>% 
          as_tibble(sep_peaks,.name_repair = "minimal")
        names(sep_peaks) <- c("mz", "intensity")
        return(sep_peaks)
      }
      
      if (spectra.file.type == "txt"){
        MSMS_Spectra <- lapply(split_list, splits) 
      }
      
      if (spectra.file.type == "msp" |spectra.file.type=="mgf"|spectra.file.type=="mgf_mzmine"){
        MSMS_Spectra <- lapply(split_list, splits_msp_mgf) 
      }
      
      MSMS_SpectraNamed <- tibble(ID = AllMatch$ID, NestedSpectra = MSMS_Spectra)
      
      MSMS_SpectraJoined <-
        left_join(
          featureDFFinal,
          MSMS_SpectraNamed,
          by = "ID"
        )
      
      DIADataObj <- addResultsNestedSpectra(DIADataObj, MSMS_SpectraJoined)
      
    }
  }