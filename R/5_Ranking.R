##### PPM window #####
ppm <- function(mz, mode, ppmtol){
  if(mode == "-"){
    mz*-ppmtol/1e6+mz
  }else{
    mz*ppmtol/1e6+mz
  }
}

##### Split MS/MS into two columns #####
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

##### Weighted dot product #####

calcProd <- function(samp, ref){
  samp_intensity <- samp$intensity
  ref_intensity <- ref$intensity
  wA <- 1/(1+(ref_intensity/(sum(ref_intensity)-0.5)))
  wB <- 1/(1+(samp_intensity/(sum(samp_intensity)-0.5)))
  dot_product <- sum(wA*wB*ref_intensity*samp_intensity)**2/( sum((wA*ref_intensity)**2) * sum((wB*samp_intensity)**2) )
  dot_product
}

calcRevDotProductNew <- function(Lipids,  max.tails=max.tails, class="TAG"){
  
  if(class!="TAG")stop("Only works for TAGs")  
  reference <- getTAGs(max.tails)
  ret <- c()
  k=3
  for (k in 1:nrow(Lipids)){
    
    dot_product <- 99
    temp_ref <- reference[reference$Name==Lipids$Name[k],]
    samp_intensity <- Lipids$MSMSref[[k]][c(1,2, 4)]
    ref_intensity <- data.frame(matrix(NA, nrow = 5, ncol = 3))
    names(ref_intensity) <- names(samp_intensity)
    
    TAG_frags <- names(temp_ref)[grep("mz", names(temp_ref))]
    ref_intensity$name <- TAG_frags[-c(2)]
    ref_intensity$mzref <- temp_ref %>% unlist(.) %>% .[ TAG_frags[-c(2)]] %>% unname() %>% as.numeric()
    ref_intensity$intensity <- c(100, 50, 999, 999, 999)
    
    samp_intensity <- samp_intensity[samp_intensity$name %in% ref_intensity$name,]
                                      
    if (any(!ref_intensity$name %in% samp_intensity$name)){
      temp <- ref_intensity[!ref_intensity$name %in% samp_intensity$name,]
      temp$intensity <- 0
      samp_intensity <- rbind(samp_intensity,temp  )
    }
    
    ref_intensity <- aggregate(ref_intensity["intensity"], by=ref_intensity["mzref"], sum)
    samp_intensity <- samp_intensity[order(samp_intensity$mz),][!duplicated(samp_intensity$mz),][,-1]
    
    ref_intensity$intensity <- ref_intensity$intensity/max(ref_intensity$intensity)*1000
    samp_intensity$intensity <- samp_intensity$intensity/max(samp_intensity$intensity)*1000
    dot_product <- calcProd(samp_intensity,ref_intensity)
    if (missing(dot_product)){
      print(k)
    }
    if (dot_product==99){
      print(k)
      stop("Nothing calculated? Inspect code.")
    }
    ret <- c(ret, dot_product)
  }
  
  Lipids$rev_dot_product <- ret
  invisible(Lipids)
}

calcRevDotProductNewAverage <- function(Lipids,  max.tails=max.tails, class="TAG"){
  
  if(class!="TAG")stop("Only works for TAGs")  
  reference <- getTAGs(max.tails)
  
  ret <- c()
  calcProd <- function(samp, ref){
    samp_intensity <- samp
    ref_intensity <- ref
    wA <- 1/(1+(ref_intensity/(sum(ref_intensity)-0.5)))
    wB <- 1/(1+(samp_intensity/(sum(samp_intensity)-0.5)))
    dot_product <- sum(wA*wB*ref_intensity*samp_intensity)**2/( sum((wA*ref_intensity)**2) * sum((wB*samp_intensity)**2) )
    dot_product
  }
  
  k=15290
  for (k in 1:nrow(Lipids)){
    
    dot_product <- 99
    temp_ref <- reference[reference$Name==Lipids$Name[k],]
    samp_intensity <- Lipids$MSMSref[[k]][c(1,2, 4)]
    ref_intensity <- data.frame(matrix(NA, nrow = 5, ncol = 3))
    names(ref_intensity) <- names(samp_intensity)
    
    TAG_frags <- names(temp_ref)[grep("mz", names(temp_ref))]
    ref_intensity$name <- TAG_frags[-c(2)]
    ref_intensity$mzref <- temp_ref %>% unlist(.) %>% .[ TAG_frags[-c(2)]] %>% unname() %>% as.numeric()
    ref_intensity$intensity <- c(100, 50, 999, 999, 999)

    samp_intensity <- Lipids$MSMSref[[k]]
    samp_intensity <- samp_intensity[samp_intensity$name %in% ref_intensity$name,]
    
    if (any(!ref_intensity$name %in% samp_intensity$name)){
      temp <- ref_intensity[!ref_intensity$name %in% samp_intensity$name,]

      temp$intensity <- 0
      temp <- unlist(temp)
      temp <- c(temp, rep(0, ncol(samp_intensity)-3))
      samp_intensity <- rbind(samp_intensity,temp )
    }
    
    ref_intensity <- aggregate(ref_intensity["intensity"], by=ref_intensity["mzref"], sum)[-1]
    samp_intensity <- samp_intensity[order(samp_intensity$mzref),][!duplicated(samp_intensity$mzref),][,-c(1,2,3,4), drop=F]
    samp_intensity[is.na(samp_intensity)] <- 0
    samp_intensity <- apply(samp_intensity, 2, as.numeric)
    samp_intensity <- apply(samp_intensity, 2, function(x) x/max(x)*1000)
    samp_intensity <- samp_intensity[,!is.na(samp_intensity[1,]), drop=F]
    
    dot_products <- apply(samp_intensity, 2, function(x) calcProd(x, ref_intensity))
    dot_product <- mean(dot_products, na.rm = T) 

    
    if (dot_product==99){
      print(k)
      stop("Nothing calculated? Inspect code.")
    }
    
    ret <- c(ret, dot_product)
  }
  
  Lipids$rev_dot_product_avg <- ret
  invisible(Lipids)
}



calcDotProductNew <- function(DIADataObj, Lipids , ppmtol=15, max.tails=max.tails, class="TAG"){
  
  if(class!="TAG")stop("Only works for TAGs")
  spectra <- DIADataObj@ResultsNestedSpectra

  
  reference <- getTAGs(max.tails)
  ret <- c()
  k=1
  for (k in 1:nrow(Lipids)){
    dot_product <- 99
    temp_ref <- reference[reference$Name==Lipids$Name[k],]
    temp_spectra <- spectra[spectra$ID==Lipids$ID.simple[k]]$NestedDeisotopedSpectra
    temp_spectra <- temp_spectra[[1]]
    

    ref_intensity <- data.frame(matrix(NA, nrow = 5, ncol = 3))
    names(ref_intensity) <- c("name", "mzref", "intensity")
    TAG_frags <- names(temp_ref)[grep("mz", names(temp_ref))]
    ref_intensity$name <- TAG_frags[-c(2)]
    ref_intensity$mzref <- temp_ref %>% unlist(.) %>% .[ TAG_frags[-c(2)]] %>% unname() %>% as.numeric()
    ref_intensity$intensity <- c(100, 50, 999, 999, 999)
    ref_intensity <- ref_intensity[order(ref_intensity$mzref),]
    
    ref_dot <- temp_spectra
    ref_dot$source <- "sample"
    ref_dot$intensity <- 0
    
    not_matched <- c()
    for (u in 1:nrow(ref_intensity)){
      match <- which(between(ref_dot$mz, 
                ppm(ref_intensity$mz[u],"-",ppmtol), 
                ppm(ref_intensity$mz[u],"+",ppmtol)))
      if (length(match>=1)){
        if (length(match)>1){
          subset <- temp_spectra[between(temp_spectra$mz, 
                          ppm(ref_intensity$mz[u],"-",ppmtol), 
                          ppm(ref_intensity$mz[u],"+",ppmtol)),]
          pick <- which.max(subset$intensity)
          match <- which(between(ref_dot$mz, 
                                 ppm(ref_intensity$mz[u],"-",ppmtol), 
                                 ppm(ref_intensity$mz[u],"+",ppmtol)))[pick]
          
        }
        ref_dot$intensity[match] <- ref_intensity$intensity[u]
      } else {
        not_matched <- c(not_matched, u)
      }
    }
    
    if(length(not_matched)!=0){
      ref_missing <- ref_intensity[not_matched, c(2,3,1)] %>% unlist()
      ref_dot <- rbind(ref_dot, ref_missing)
      ref_missing[2] <- 0
      ref_missing <- ref_missing[1:2]
      samp_intensity <- rbind(temp_spectra,ref_missing)
    } else{
      samp_intensity <- temp_spectra
    }
    
    ref_intensity <- ref_dot
    if (nrow(ref_intensity)!=nrow(samp_intensity)){
      print(k)
      stop("Different lengths? Inspect code.")
    }
    
    samp_intensity$intensity <- as.numeric(samp_intensity$intensity)
    ref_intensity$intensity <- as.numeric(ref_intensity$intensity)
  
    ref_intensity$intensity <- ref_intensity$intensity/max(ref_intensity$intensity)*1000
    samp_intensity$intensity <- samp_intensity$intensity/max(samp_intensity$intensity)*1000
    dot_product <- calcProd(samp_intensity,ref_intensity)
    if (dot_product==99){
      print(k)
      stop("Nothing calculated? Inspect code.")
    }
    ret <- c(ret, dot_product)
  }
  
  Lipids$dot_product <- ret
  invisible(Lipids)
}


calcRevDotProduct <- function(Lipids){
  
  calcProd <- function(samp, ref){
    samp_intensity <- samp
    ref_intensity <- ref
    wA <- 1/(1+(ref_intensity/(sum(ref_intensity)-0.5)))
    wB <- 1/(1+(samp_intensity/(sum(samp_intensity)-0.5)))
    dot_product <- sum(wA*wB*ref_intensity*samp_intensity)**2/( sum((wA*ref_intensity)**2) * sum((wB*samp_intensity)**2) )
    dot_product
  }
  
  
  ret <- c()
  for (k in 1:nrow(Lipids)){
    
    samp_intensity <- Lipids$MSMSref[[k]][grep("Fragment",  Lipids$MSMSref[[k]]$name),]
    ref_intensity <- samp_intensity
    ref_intensity$intensity <- 1
    ref_intensity <- aggregate(ref_intensity["intensity"], by=ref_intensity["mz"], sum)
    
    samp_intensity <- samp_intensity[order(samp_intensity$mz),][!duplicated(samp_intensity$mz),]
    if(nrow(samp_intensity)<2){
      ret <- c(ret, NA)
      next
    }
    
    ref_intensity <- ref_intensity$intensity/max(ref_intensity$intensity)
    samp_intensity <- samp_intensity[,-c(1:4), drop=F]
    samp_intensity <- apply(samp_intensity, 2, function(x) x/max(x))
    
    dot_products <- apply(samp_intensity, 2, function(x) calcProd(x, ref_intensity))
    dot_product <- mean(dot_products, na.rm = T) 
    ret <- c(ret, dot_product)
  }
  Lipids$Rev.dot.product <- ret
  invisible(Lipids)
}

calcDotProductOriginalSingle <- function(Lipids){
  
  calcProd <- function(samp, ref){
    samp_intensity <- samp
    ref_intensity <- ref
    wA <- 1/(1+(ref_intensity/(sum(ref_intensity)-0.5)))
    wB <- 1/(1+(samp_intensity/(sum(samp_intensity)-0.5)))
    dot_product <- sum(wA*wB*ref_intensity*samp_intensity)**2/( sum((wA*ref_intensity)**2) * sum((wB*samp_intensity)**2) )
    dot_product
  }
  
  ret <- c()
  k=1
  for (k in 1:nrow(Lipids)){
    
    samp_intensity <- Lipids$MSMSref[[k]][grep("Fragment",  Lipids$MSMSref[[k]]$name),]
    ref_intensity <- samp_intensity
    ref_intensity$intensity <- 1
    ref_intensity <- aggregate(ref_intensity["intensity"], by=ref_intensity["mz"], sum)
    
    samp_intensity <- samp_intensity[order(samp_intensity$mz),][!duplicated(samp_intensity$mz),]
    if(nrow(samp_intensity)<2){
      ret <- c(ret, NA)
      next
    }
    
    ref_intensity <- ref_intensity$intensity/max(ref_intensity$intensity)
    samp_intensity <- samp_intensity[,3, drop=F]
    samp_intensity <- apply(samp_intensity, 2, function(x) x/max(x))
    
    dot_product <- calcProd(samp_intensity, ref_intensity)
    ret <- c(ret, dot_product)
  }
  Lipids$dot_product_original_single <- ret
  invisible(Lipids)
}


getIntensity <- function(DIADataObj = DIA_Pos, Lipids=Identified_TAG_Combined, run.code="TJV"){

    Alignment <- DIADataObj@Results
    Lipids$MS1Intensity <- NA
    
    for (q in 1:nrow(Lipids)){
      
      if(Lipids$FeatureType[q]=="Lipid"){
        
        samples <- Lipids$MSMSref[[q]] %>% column_to_rownames(.,var="name") %>% .[, grep("sample", colnames(.)), drop=F] %>% t() %>% as.data.frame()
        fragment_index <- grep("ment", Lipids$MSMSref[[q]]$name)
        fragments <- samples[,c(fragment_index)]
        
        no_na_sample <- fragments[complete.cases(fragments),]
        if(nrow(no_na_sample)==0){
          Lipids$MS1Intensity[q] <- 0
          next
        }
        #Cannot calculate correlation if there are less than 3 samples
        
        #Matches MS1 heights in alignment file with MS2 fragments 
        #Matches by selecting sample number of each row. 
        #Ensure alignment file heights in order and complete. 
        heights <- Alignment[Alignment$ID==Lipids$ID.simple[[q]],]
        select_cols <- grep(run.code, names(heights))
        heights_sample <- heights[,..select_cols]
        sort_names <- str_sort(names(heights_sample))
        heights_sample <- heights_sample[,..sort_names]
        
        no_na_sample$PrecursorMS1 <- NA
        
        for (i in 1:nrow(no_na_sample)){
          temp <- str_extract(rownames(no_na_sample)[i], "[0-9]+") %>% as.numeric()
          no_na_sample$PrecursorMS1[i] <- heights_sample[[temp]]
        }
        no_na_sample$PrecursorMS1 <- as.numeric(no_na_sample$PrecursorMS1)

        Lipids$MS1Intensity[q] <- no_na_sample$PrecursorMS1 %>% max()
      }
      
    }
    invisible(Lipids)
}

getTails <-  function(Lipids){
  ret <- c()
  for (k in 1:nrow(Lipids)){
     samp_intensity <- Lipids$MSMSref[[k]][grep("Fragment",  Lipids$MSMSref[[k]]$name),]
     samp_intensity <- samp_intensity[order(samp_intensity$mz),][!duplicated(samp_intensity$mz),]
     tails <- nrow(samp_intensity)
     ret <- c(ret, tails)
  }
  Lipids$Tails <- ret
  invisible(Lipids)
}



## 270 seconds for 271 rows of PCS
 # spectra.file = "NewData/MSe_Pos_Data/Centroid"
 # Lipids = TAG

##### Retrieve spectra from every sample #####
retrieveSpectra <- function(Lipids, spectra.file, spectra.file.type=c("txt","mgf","msp"), ppmtol, rttolMS1=12){
  parse_msp <- function(file_path) {
    
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
    end_indices <- c(tail(indices, -1) - 5, length(data))
    fragments_temp <- sapply(seq_along(start_indices), function(i) {
      paste(msp[start_indices[i]:end_indices[i]], collapse = " ")
    })
    
    
    spectra_file <- data.frame(rt = rt_temp, mz = mz_temp,msms=fragments_temp)
    
    return(spectra_file)
  }
  parse_mgf <- function(file_path) {
    
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
  
  
  #We decrease the ppmtolMS1 because we want to narrow down the features we are matching.
  #We increase the ppmtol for MS/MS because we already know the peaks exist in one sample.
  ppmtolMS1 <- ppmtol

  spectra.file <- list.files(spectra.file, full.names = T)
  lipid_ID <- unique(Lipids$ID.simple)
  pb <- txtProgressBar(min = 1, max = length(spectra.file)+1, style = 3, width = 100)

  
  for (k in 1:length(spectra.file)){ 
    
    if (spectra.file.type=="msp"){
      sample <- parse_msp(spectra.file[k])
      colnames(sample)[colnames(sample) == "rt"] <- "RT (min)"
      colnames(sample)[colnames(sample) == "mz"] <- "Precursor m/z"
      colnames(sample)[colnames(sample) == "msms"] <- "MSMS spectrum"
    }
    if (spectra.file.type=="mgf"){
      sample <- parse_mgf(spectra.file[k])
      colnames(sample)[colnames(sample) == "rt"] <- "RT (min)"
      colnames(sample)[colnames(sample) == "mz"] <- "Precursor m/z"
      colnames(sample)[colnames(sample) == "msms"] <- "MSMS spectrum"
    }
    if (spectra.file.type=="txt"){
      sample <- read_delim(spectra.file[k],show_col_types = F) %>% as.data.frame()
    }   
    # #We load spectra samples one at a time because of file size.
    # sample <- read_delim(spectra.file[k],show_col_types = F) %>% as.data.frame()

    for (j in 1:length(lipid_ID)){
      
      #We loop through each unique species. 
      species <- Lipids[Lipids$ID.simple==lipid_ID[j],]
      setTxtProgressBar(pb, k+j/length(lipid_ID))
      
      if (species$FeatureType[1]=="Lipid"){
        
        #Filter by ppm and rt 
        match <- sample[(between(sample$`Precursor m/z`, ppm(species$mz[1],"-",ppmtolMS1), ppm(species$mz[1],"+",ppmtolMS1)) & 
                            between(sample$`RT (min)`,species$rt[1]-rttolMS1/60, species$rt[1]+rttolMS1/60)),]
  
          if (nrow(match)>1) {
            rttolr <- rttolMS1
            
            #Keep reducing retention time tolerance until there are no matches (by 4 for speed).
            repeat{
              rttolr <- rttolr-4
              matchtemp <- match[between(match$`RT (min)`,species$rt[1]-rttolr/60, species$rt[1]+rttolr/60),]
              
              if (rttolr<=4 | nrow(matchtemp)==1) {
                match <- matchtemp
                break
              } else if (nrow(matchtemp)==0){
                
                #Keep increasing retention time tolerance until there is one or more matches. 
                repeat{ 
                  rttolr <- rttolr+1
                  matchtemp <- match[between(match$`RT (min)`,species$rt[1]-rttolr/60, species$rt[1]+rttolr/60),]
                 
                  if (nrow(matchtemp)>=1) {
                    match <- matchtemp
                    break
                  }
                }
              }
              break
            }
            
            #Select larger area if there is still more than one match. 
            if (nrow(match)>1) {match <- match[which.min(abs(species$mz[1]-match$`Precursor m/z`)),]}
          }
          
          #Now we should have the MS/MS where the lipid was identified
          if (nrow(match)==1) {
            if(spectra.file.type=="txt"){
              sample_msms <- splits(match$`MSMS spectrum`)
              }
            if((spectra.file.type=="msp") |(spectra.file.type=="mgf")){
              sample_msms <- splits_msp_mgf(match$`MSMS spectrum`)
            }
              
            #Looping through every lipid identified from that feature 
            for (p in 1:nrow(species)){
            
              sample_ref <- species$MSMSref[[p]]
              
              temp <- data.frame()

              #Looping through fragments and precursor mz
              for (l in 1:nrow(sample_ref)){
                
                #Subset before finding the closest one (for speed). 
                between_msms <- sample_msms[between(sample_msms$mz, sample_ref$mzref[l]-1, sample_ref$mzref[l]+1),]
                
                #Prevents errors if the first row is empty
                if(nrow(between_msms)==0){
                  fill_it <- data.frame(mz=NA, intensity=NA)
                  temp <- rbind(temp, fill_it)
                  next
                }
                
                #Select peak that is the closest to the reference.
                mz_difference <- abs(sample_ref$mzref[l]-between_msms$mz)
                smallest_col <- which.min(mz_difference)
                temp <- rbind(temp,between_msms[smallest_col,])
              }
              
              #Set peak as NA if the closest peak is not within the ppm window
              boolean <- map2_lgl(.x=temp$mz, .y=sample_ref$mzref, ~between(.x, ppm(.y,"-",ppmtol), ppm(.y,"+",ppmtol)))
              removed_ppm_limit <- map2_dbl(.x=temp$intensity, .y=boolean, ~ifelse(.y, .x, NA))
              Lipids[Lipids$ID.simple==lipid_ID[j],]$MSMSref[[p]][paste0("sample_",k)] <-  removed_ppm_limit
            }
          } 
        }
    }
    gc()
    }
  invisible(Lipids)
}



calcCorrelationPrecursors <- function(Lipids) {
  Lipids$Correlation <- NA
  Lipids$Correlation.tan <- NA
  Lipids$Missing.fragments <- NA
  Lipids$Missing.fragments.one <- NA
  
  for (q in 1:nrow(Lipids)){
    
    if(Lipids$FeatureType[q]=="Lipid"){
      
      samples <- Lipids$MSMSref[[q]] %>% column_to_rownames(.,var="name") %>% .[, grep("sample", colnames(.)), drop=F] %>% t() %>% as.data.frame()
      precursors_index <- grep("recur", Lipids$MSMSref[[q]]$name)
      fragment_index <- grep("ment", Lipids$MSMSref[[q]]$name)
      no_na_sample <- samples[complete.cases(samples),]
      
      if (length(precursors_index)>2) { print("There is a two precursor limit for this function."); break}
      
      #Sum missing fragments (NAs)
      oneornone <- apply(samples,1, function(x) sum(is.na(x)))==1|apply(samples,1, function(x) sum(is.na(x)))==0
      Lipids$Missing.fragments.one[q]<- sum(is.na(samples[oneornone,]))
      Lipids$Missing.fragments[q] <- sum(is.na(samples))
      
      #Cannot calculate correlation if there are less than 3 samples 
      if(nrow(no_na_sample)<3){
        Lipids$Correlation[q] <- -1
        Lipids$Correlation.tan[q] <- -1
        next}  
      
      corr_mat <- cor(samples, use = "pairwise.complete.obs")
      Lipids$Correlation[q] <- mean(corr_mat[fragment_index, precursors_index])
      
      #We want to Fisher z-transform the correlation when averaging
      temp <- corr_mat[fragment_index, precursors_index]
      Lipids$Correlation.tan[q] <- atanh(temp) %>% mean %>% tanh
      
    }
  }
  invisible(Lipids)
}

calcCorrelationWithin <- function(Lipids,ms2.precursors=ms2.precursors) {
  # Lipids$Correlation.max.pf <- NA
  # Lipids$Correlation.max.ff <- NA
  
  Lipids$Correlation.tan.within <- NA
  Lipids$Missing.fragments <- NA
  # Lipids$Missing.fragments.one <- NA
  q=2
  for (q in 1:nrow(Lipids)){
    
    if(Lipids$FeatureType[q]=="Lipid"){
      
      samples <- Lipids$MSMSref[[q]] %>% column_to_rownames(.,var="name") %>% .[, grep("sample", colnames(.)), drop=F] %>% t() %>% as.data.frame()
      precursors_index <- grep("recur", Lipids$MSMSref[[q]]$name)
      fragment_index <- grep("ment", Lipids$MSMSref[[q]]$name)

      # if (length(precursors_index)>2) { print("There is a two precursor limit for this function."); break}
      
      #Sum missing fragments (NAs)
      # oneornone <- apply(samples,1, function(x) sum(is.na(x)))==1|apply(samples,1, function(x) sum(is.na(x)))==0
      # Lipids$Missing.fragments.one[q]<- sum(is.na(samples[oneornone,]))
      if (ms2.precursors == 1){
        precursors_fragments <- samples[,c(fragment_index, precursors_index)]
        Lipids$Missing.fragments[q] <- sum(is.na(precursors_fragments))
        # 
        # if (length(precursors_index)==1) {
        #   Lipids$Missing.fragments[q] <- Lipids$Missing.fragments[q] + 1
        # }
        # # if (length(precursors_index)==2) {
        # #   Lipids$Missing.fragments[q] <- Lipids$Missing.fragments[q] + 1 
        # # }
      }
      else {
        fragments <- samples[,c(fragment_index)]
        Lipids$Missing.fragments[q] <- sum(is.na(fragments))
      }
      
  
      frags_only <- samples[,fragment_index]
      # if (length(precursors_index)==1) {
      #   Lipids$Missing.fragments[q] <- Lipids$Missing.fragments[q] + 1 
      # }
      no_na_sample <- samples[complete.cases(frags_only),]
      
      #Cannot calculate correlation if there are less than 3 samples 
      if(nrow(no_na_sample)<3){
        # Lipids$Correlation.max.ff[q] <- -1
        # Lipids$Correlation.max.pf[q] <- -1
        Lipids$Correlation.tan.within[q] <- -1
        next}  
      
      corr_mat <- cor(frags_only, use = "pairwise.complete.obs")
      # Lipids$Correlation.max.pf[q] <- corr_mat[fragment_index, precursors_index][order(corr_mat[fragment_index, precursors_index], decreasing = T)][1]
      # temp.corr <- corr_mat[fragment_index, fragment_index]
      # Lipids$Correlation.max.ff[q] <- temp.corr %>% unname() %>% c() %>% .[temp.corr!=1] %>% .[order(., decreasing = T)] %>% .[1]
      #Mean of correlation matrix
      #We want to Fisher z-transform the correlation when averaging
      # temp <- corr_mat[fragment_index, precursors_index]
      
      upper_matrix <- corr_mat[upper.tri(corr_mat)] 
      upper_matrix <- upper_matrix[!is.na(upper_matrix)]
      if (any(upper_matrix==1)){
        upper_matrix[upper_matrix==1] <- 0.99999
      }
      Lipids$Correlation.tan.within[q] <-unique(upper_matrix) %>% atanh %>% mean %>% tanh
      
    }
  }
  invisible(Lipids)
}

calcCorrelation <- function(DIADataObj = DIA_Pos, Lipids, run.code="Sample") {
  # Lipids$Correlation.max.pf <- NA
  # Lipids$Correlation.max.ff <- NA
  Alignment <- DIADataObj@Results
  Lipids$Correlation.tan <- NA
  # Lipids$Missing.fragments <- NA
  q=1
  for (q in 1:nrow(Lipids)){
    
    if(Lipids$FeatureType[q]=="Lipid"){
      
      samples <- Lipids$MSMSref[[q]] %>% column_to_rownames(.,var="name") %>% .[, grep("sample", colnames(.)), drop=F] %>% t() %>% as.data.frame()
      precursors_index <- grep("recur", Lipids$MSMSref[[q]]$name)
      fragment_index <- grep("ment", Lipids$MSMSref[[q]]$name)
      
      fragments <- samples[,c(fragment_index)]
      # Lipids$Missing.fragments[q] <- sum(is.na(fragments))
      

      no_na_sample <- fragments[complete.cases(fragments),]
      
      #Cannot calculate correlation if there are less than 3 samples 
      if(nrow(no_na_sample)<3){
        # Lipids$Correlation.max.ff[q] <- -1
        # Lipids$Correlation.max.pf[q] <- -1
        Lipids$Correlation.tan[q] <- -1
        next
      }  
      
      #Matches MS1 heights in alignment file with MS2 fragments 
      #Matches by selecting sample number of each row. 
      #Ensure alignment file heights in order and complete. 
      heights <- Alignment[Alignment$ID==Lipids$ID.simple[[q]],]
      select_cols <- grep(run.code, names(heights))
      heights_sample <- heights[,..select_cols]
      sort_names <- str_sort(names(heights_sample))
      heights_sample <- heights_sample[,..sort_names]
      
      no_na_sample$PrecursorMS1 <- NA
  
      for (i in 1:nrow(no_na_sample)){
        temp <- str_extract(rownames(no_na_sample)[i], "[0-9]+") %>% as.numeric()
        no_na_sample$PrecursorMS1[i] <- heights_sample[[temp]]
      }
      no_na_sample$PrecursorMS1 <- as.numeric(no_na_sample$PrecursorMS1)
      corr_mat <- cor(no_na_sample, use = "pairwise.complete.obs") %>% as.data.frame()
      result <- corr_mat$PrecursorMS1[-length(corr_mat$PrecursorMS1)]
      result[result==1] <- 0.99999
      Lipids$Correlation.tan[q] <-  unique(result) %>% atanh %>% mean %>% tanh
      
    }
  }
  invisible(Lipids)
}



calcEntropy <- function(spectra){
  temp <- -spectra*log(spectra) 
  ret <- sum(temp)
  return(ret)
}

calcEntropyScore <- function(DIADataObj = DIA_Pos, Lipids, ppmtol=15, max.tails=max.tails, class="TAG") {
  
  if(class!="TAG") {stop("Only works for TAGs")}
  spectra <- DIADataObj@ResultsNestedSpectra
  Lipids$Entropy <- NA
  Lipids$Entropy_similarity <- NA
  
  reference <- getTAGs(max.tails)
  
  ret <- c()
  k=1
  for (k in 1:nrow(Lipids)){
    temp_spectra <- spectra[spectra$ID==Lipids$ID.simple[k]]$NestedDeisotopedSpectra
    temp_spectra <- temp_spectra[[1]]
    normalized_spectra_int <- temp_spectra$intensity/(sum(temp_spectra$intensity))
    temp_spectra$intensity <- normalized_spectra_int
    
    spectra_entropy <- calcEntropy(normalized_spectra_int)
    Lipids$Entropy[k] <- spectra_entropy
    
    temp_ref <- reference[reference$Name==Lipids$Name[k],]
    ref_intensity <- data.frame(matrix(NA, nrow = 5, ncol = 3))
    names(ref_intensity) <- c("name", "mzref", "intensity")
    TAG_frags <- names(temp_ref)[grep("mz", names(temp_ref))]
    ref_intensity$name <- TAG_frags[-c(2)]
    ref_intensity$mzref <- temp_ref %>% unlist(.) %>% .[ TAG_frags[-c(2)]] %>% unname() %>% as.numeric()
    ref_intensity$intensity <- c(100, 50, 999, 999, 999)
    normalized_spectra_int <- ref_intensity$intensity/(sum(ref_intensity$intensity))
    ref_intensity$intensity <- normalized_spectra_int
    ref_entropy <- calcEntropy(normalized_spectra_int)
    
    ref_dot <- temp_spectra
    not_matched <- c()
    for (u in 1:nrow(ref_intensity)){
      match <- which(between(ref_dot$mz, 
                             ppm(ref_intensity$mz[u],"-",ppmtol), 
                             ppm(ref_intensity$mz[u],"+",ppmtol)))
      if (length(match>=1)){
        if (length(match)>1){
          subset <- temp_spectra[between(temp_spectra$mz, 
                                         ppm(ref_intensity$mz[u],"-",ppmtol), 
                                         ppm(ref_intensity$mz[u],"+",ppmtol)),]
          pick <- which.max(subset$intensity)
          match <- which(between(ref_dot$mz, 
                                 ppm(ref_intensity$mz[u],"-",ppmtol), 
                                 ppm(ref_intensity$mz[u],"+",ppmtol)))[pick]
          
        }
        ref_dot$intensity[match] <-ref_dot$intensity[match]+ref_intensity$intensity[u]
      } else {
        not_matched <- c(not_matched, u)
      }
    }
    if(length(not_matched)!=0){
      ref_missing <- ref_intensity[not_matched, c(2,3)] %>% unlist()
      ref_dot <- rbind(ref_dot, ref_missing)
    }
    
    normalized_combined_int <- combined$intensity/(sum(combined$intensity))
    combined_entropy <- calcEntropy(normalized_combined_int)
    unweighted_entropy <- 1 - (2 * combined_entropy - ref_entropy - spectra_entropy)/log(4)
    Lipids$Entropy_similarity[k] <- unweighted_entropy
  }

  invisible(Lipids)
}


