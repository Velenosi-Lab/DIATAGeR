##### Max carbons for isotope generation #####
calcMaximumCarbon <- function(masses) { 
  
  oxygen = 15.999
  carbon = 12.011
  hydrogen  = 1.008
  CH2 = carbon + 2 * hydrogen
  
  #All lipids have at least a mass of oxygen contributing to their mass
  maximum_carbon <- floor((masses-oxygen)/CH2)
  
  #Original formula
  #maximum_carbon <- floor((masses-2*CH3-oxygen)/CH2) + 2
  
  return(maximum_carbon)
}

##### Density cutoff #####
calcDensityChange <- function(spectra, threshold=density.cutoff.value,
                              cutoff.floor=low.intensity.cutoff, cutoff.ceiling=1500) {
  
  if (sum(spectra$intensity>=1e5)>20) { 
    #If many intensities above 1e5, remove top 5 largest
    spectra_ceiling <- spectra[spectra$intensity< min(tail(sort(spectra$intensity), 5)),]
    points <- 2**20
  } else {
    #Else remove everything above 1e5 
    spectra_ceiling <- spectra[spectra$intensity<1e5,]
    points <- 2**18 
  }
  
  #Density estimate for intensity distribution
  density_calc <- density(spectra_ceiling$intensity, n=points)
  # For plotting
  # AAA<<-density_calc
  xx <- density_calc$x  
  yy <- density_calc$y 
  increment <- xx[2] - xx[1]
  
  #Range of possible intensity cutoffs by 10
  seqence <- seq(cutoff.floor, cutoff.ceiling, by=10)
  
  ret <- c()
  for (cutoff in seqence) {
    #Relative change is (y2-y1)/y1*100, where y1 the Riemann sum from 0 to the cutoff
    #and y2 the Riemann sum from 0 to the cutoff+10.
    relative_change <- ((sum(yy[xx < (cutoff+10)])*increment - sum(yy[xx < cutoff])*increment) / 
                        (sum(yy[xx < cutoff])*increment))*100
    ret <- c(ret, relative_change)
  }
  
  #If relative change greater less than 1%, we set as cutoff. 
  cutoff_df <- data.frame(change=ret, cutoff=seqence)
  thresh <- min(cutoff_df$cutoff[cutoff_df$change<threshold])
  
  #Returns intensity cutoff
  return(thresh)
}

##### Natural abundance of carbon 13 #####
carbonisotopeabundance = 0.0107 

##### Isotope distribution #####
calcIsotopes <- function(spectra){
  
  carbons <- calcMaximumCarbon(spectra$mz)
  to_normalize <- dbinom(0, carbons, carbonisotopeabundance)
  int1 <- dbinom(1, carbons, carbonisotopeabundance)/to_normalize #Isotope peak 1 %
  int2 <- dbinom(2, carbons, carbonisotopeabundance)/to_normalize #Isotope peak 2 %
  
  iso_int1 <-  spectra$intensity * int1
  iso_int2 <-  spectra$intensity * int2
  mz_1 <- spectra$mz + 1.003355
  mz_2 <-  spectra$mz + 1.003355*2
  mz_1_up <-  mz_1*ppmtol/1e6+mz_1+0.001
  mz_1_down <-  mz_1*-ppmtol/1e6+mz_1-0.001
  mz_2_up <-  mz_2*ppmtol/1e6+mz_2+0.001
  mz_2_down <-  mz_2*-ppmtol/1e6+mz_2-0.001
  
  ret <- cbind(spectra,
               iso_int1=iso_int1, iso_int2=iso_int2,
               mz_1=mz_1, mz_2=mz_2,
               mz_1_up=mz_1_up, mz_1_down=mz_1_down, 
               mz_2_up=mz_2_up, mz_2_down=mz_2_down)
  
  return(ret)
}


#' Deisotoping and denosing MS2 spectra
#' 
#' @param DIADataObj object stores feature list and MS2 spectra.
#' @param ppmtol ppm tolerance. Default 15. 
#' @param low.intensity.cutoff Anything below this value will be removed. Default 100.
#' @param density.cutoff If TRUE, a density-based cutoff will be applied. The kernel
#' density estimate of each spectra's intensity distribution will be used to implement
#' a low intensity cutoff. 
#' @param density.cutoff.value The low intensity cutoff is the intensity where the
#' percentage change in the the Riemann sum of the kernel density when the intensity
#' is increased by 10 drops below the density.cutoff.value. Higher value for a lower intensity 
#' cutoff. Default 1. Change with caution. 
#' @param lenient.removal If TRUE, removes M+1 and M+2 isotope peaks. Default TRUE
#' @param check.isotopes If TRUE, check the isotopes in the spectra. Default TRUE
#'
#' @return De-isotoped spectra
#' @export 
#'
#' @examples
#' \dontrun{
#' DIA_Pos <- KeepIsotopes(DIADataObj = DIA_Pos, 
#' ppmtol = 15,
#' low.intensity.cutoff=100,
#' density.cutoff=T,
#' density.cutoff.value=1)
#' }

RemoveIsotopes <- function(DIADataObj, ppmtol = 15, 
                         low.intensity.cutoff=100, density.cutoff=T, density.cutoff.value=1,
                         check.isotopes=T, lenient.removal=T ) {
  featureDF_Deisotoped <- DIADataObj@ResultsNestedSpectra
  p=0
  ## Selects identical spectra based on identical retention times and spectra lengths
  rows_spectra <- c()
  
  for (i in 1:nrow(featureDF_Deisotoped)){
    rows_spectra <- c(rows_spectra, nrow(featureDF_Deisotoped$NestedSpectra[[i]]))
  }
  
  unique_spectra  <- data.frame(rt=featureDF_Deisotoped$rt, length=rows_spectra,
                                unique_id=paste(featureDF_Deisotoped$rt, rows_spectra, sep = "-"))
  unique_ids <- unique(unique_spectra$unique_id)
  SpectraList <- featureDF_Deisotoped$NestedSpectra
  names(SpectraList) <- unique_spectra$unique_id
  ## 
  
  pb <- txtProgressBar(min = 1, max = length(unique_ids), style = 3)
  to_keep_list <- vector("list", length = length(SpectraList))
  environment(calcIsotopes) <- environment()
  environment(calcDensityChange) <- environment()
  

  #Example i=2875
  #Example i=4242 18.492-25463
  for(i in 1:length(unique_ids)) {
   
    ## Some test cases
    # MSMS_subset <- data.frame(mz=c(150, 151.003355, 152.00671), intensity=c(10000, 990, 150)) #above +1 and +2  
    # MSMS_subset <- data.frame(mz=c(150, 151.003355, 152.00671), intensity=c(10000, 900, 140)) #within +1 and +2
    # MSMS_subset <- data.frame(mz=c(150, 151.003355, 152.00671), intensity=c(10000, 300, 140)) #below  +1 and within +2
    ## 
    
    #Select spectra with identical IDs
    same_spectra <- which(names(SpectraList)==unique_ids[i]) 
    MSMS_subset <- SpectraList[[same_spectra[1]]] 
    
    #Internal function for setup
    MSMS_subset <- calcIsotopes(MSMS_subset)

    #We don't bother checking for isotopes if the intensity is below the low.intensity.cutoff, which also 
    #means nothing under the low.intensity.cutoff will be kept. (Speed) 
    MSMS_over_cutoff <- MSMS_subset[MSMS_subset$intensity>low.intensity.cutoff,]
    
    #Isotopes have a greater m/z at +1 and +2. Therefore, we only check a certain number of rows (add_index) beyond
    #the row being tested in MSMS_over_cutoff. (Speed) 
    add_index <- as.integer(nrow(MSMS_subset)/(max(MSMS_subset$mz)-min(MSMS_subset$mz))*10)+10

    to_keep <- NULL
    to_remove <- NULL

    if(is.na(MSMS_over_cutoff$mz[1])){
      next
    }
    
    if(check.isotopes==F){
      
      MSMS_kept <- MSMS_over_cutoff[,1:2]
    
    } else {
      
      for(j in 1:nrow(MSMS_over_cutoff)){
        
        subset_it <- as.integer(row.names(MSMS_over_cutoff)[j])
        upper <- subset_it+add_index
        
        if (upper > nrow(MSMS_subset)) {
          upper <- nrow(MSMS_subset)
        }
        
        mz1_test <- between(MSMS_subset$mz[subset_it:upper], MSMS_over_cutoff$mz_1_down[j], MSMS_over_cutoff$mz_1_up[j])
        
        if (sum(mz1_test)!=0) {
          
          #We only generate object if TRUE. (Speed) 
          MSMS_subset_temp <- MSMS_subset[subset_it:upper,]
          
          iso1_test <- MSMS_subset_temp[mz1_test,]
          iso1_bool <- iso1_test$intensity > MSMS_over_cutoff$iso_int1[j]/2.5
          
          if (any(iso1_bool)) {
  
            to_keep <- c(to_keep, rownames(MSMS_over_cutoff)[j])
            
            #Test if it is M+1 isotope, if it is between maximum and minimum intensity. 
            iso1_remove_bool <- between(iso1_test$intensity, MSMS_over_cutoff$iso_int2[j]/2.5, MSMS_over_cutoff$iso_int1[j])
            
            if(any(iso1_remove_bool)){
              to_remove <- c(to_remove, rownames(iso1_test))
            }
            
            #Test if it is M+2 isotope.
            mz2_test <- between(MSMS_subset_temp$mz, MSMS_over_cutoff$mz_2_down[j], MSMS_over_cutoff$mz_2_up[j])
         
           if (sum(mz2_test)!=0){
             
             iso2_test <- MSMS_subset_temp[mz2_test,]
             #Test if it is M+2 isotope, if it is between maximum and minimum intensity. 
             iso2_bool <- between(iso2_test$intensity, MSMS_over_cutoff$iso_int2[j]/2.5, MSMS_over_cutoff$iso_int2[j])
    
              if (any(iso2_bool)) {
               to_remove <- c(to_remove,  rownames(iso2_test))
              }
            }
          }
        }
      }
    
      to_keep <- as.integer(to_keep)
      to_remove <- as.integer(to_remove)
      
      MSMS_kept <- MSMS_subset[rownames(MSMS_subset) %in% to_keep,][,1:2]
    }
    
    if (density.cutoff) {
      intensity_cutoff <- calcDensityChange(MSMS_subset)
      MSMS_kept <- MSMS_kept[MSMS_kept$intensity>intensity_cutoff,]
    }
    
  ## lenient.removal if true, removes x if +1 is below max carbon thresh.
  if (lenient.removal){
    if(length(to_remove)!=0){
      MSMS_kept <- MSMS_kept[!(rownames(MSMS_kept) %in% to_remove),]
    } 
  }
    
  # Makes plot of density cutoff and removed spectra.
  #   plot_MSMS <- MSMS_kept
  #   plot_MSMS$intensity <-  -plot_MSMS$intensity
  #   plot_ORI <- MSMS_subset
  # 
  #   A <- ggplot()+
  #   geom_bar(data = plot_ORI, aes(x=mz, y=intensity), width = 0.001, col="#03053A", stat="identity",position = "dodge")+
  #   geom_bar(data = plot_MSMS, aes(x=mz, y=intensity), width = 0.001, col="#03053A", stat="identity",position = "dodge")+
  #   geom_text(data = plot_ORI, hjust = 0, size= 4, aes(y=50000, x=455, label=paste("Raw spectra\nPeaks:", nrow(MSMS_subset))))+
  #   geom_text(data = plot_MSMS,hjust = 0,size= 4,  aes(y=-50000, x=455, label=paste("Peaks with isotopes\nPeaks:", nrow(plot_MSMS))))+
  #   xlim(c(450, 780))+
  #   # geom_hline(yintercept=0, size = 0.45, linetype="dashed")+
  #   theme_bw()+
  #   theme(plot.title = element_text(family = 'Helvetica', color = '#666666', face = 'bold', size = 14),
  #           panel.border = element_blank(),
  #           axis.line = element_line(colour = "black", size = 0.45),
  #           panel.background = element_rect(fill = "transparent", colour = NA),
  #           plot.background = element_rect(fill = "transparent", colour = NA),
  #           panel.grid.minor = element_blank(),
  #           panel.grid.major =  element_blank(),
  #           plot.margin = unit(c(0.5,1,0,0.74), "cm"))+
  #     ylab("Intensity")+
  #     xlab("m/z")
  # 
  #   B <- ggplot()+
  #     geom_area(aes(x=AAA$x,y=AAA$y),fill="#03053A", alpha=0.7)+
  #     geom_vline(xintercept = intensity_cutoff)+
  #     theme_bw()+
  #       xlim(c(0, 1200))+
  #       theme(plot.title = element_text(family = 'Helvetica', color = '#666666', face = 'bold', size = 14),
  #           panel.border = element_blank(),
  #           axis.line = element_line(colour = "black", size = 0.45),
  #           panel.background = element_rect(fill = "transparent", colour = NA),
  #           plot.background = element_rect(fill = "transparent", colour = NA),
  #           panel.grid.minor = element_blank(),
  #           panel.grid.major =  element_blank(),
  #           plot.margin = unit(c(0,1,1,1), "cm"))+
  #     xlab("Intensity")+
  #     ylab("Density")
  # 
  # library(cowplot)
  # P <- plot_grid(A, B, labels=c("A", "B"), ncol = 1, nrow = 2, rel_heights = c(3,1))
  # 
  # ggsave(
  #   "double1.pdf",
  #   plot = P ,
  #   scale = 1,
  #   width = 9,
  #   height = 9,
  #   units = c("in")
  # )

    to_keep_list[same_spectra] <-  lapply(to_keep_list[same_spectra], function(x) x <- MSMS_kept)

    if(nrow(MSMS_kept)<10){
      print(paste("Only", nrow(MSMS_kept), "peaks kept. Try higher density.cutoff.value or alternative to keeping peaks with an isotope?"))
      p <- p+1
    }

    setTxtProgressBar(pb, i)
  } 

    MSMS_Deisotope_Nest <- tibble(ID = featureDF_Deisotoped$ID,
                                  NestedDeisotopedSpectra = to_keep_list
    )
    
    featureDF_Deisotoped <- left_join(
      featureDF_Deisotoped,
      MSMS_Deisotope_Nest,
      by = "ID"
    )
    print( paste0(p, "spectra with less than 10 spectra."))
    addResultsNestedSpectra(DIADataObj, featureDF_Deisotoped)
}
