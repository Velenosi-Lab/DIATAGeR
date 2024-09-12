### Decoys for all molecular species '

# DAGs <- readRDS("Pos/Identified_DAG_Pos_TAGsJune22_1") 
# Identified <- readRDS("Pos/Identified_TAG_Pos_Jul9_1")

#' Create decoy database
#'
#' @param Identified Identified lipids file 
#' @param version Appends this character to the output file name of identified lipids. 
#' @param max.tails list of desired number of carbons and the maximum number of double 
#' bonds for each fatty acyl chain length. 
#' @param ion.mode Pos or Neg to indicate positive or negative ion mode
#' @param exact.tails library can be customizable for specific tails 
#'
#' @return a library of decoy triacylglycerols
#' @export
#'
#' @examples
#' \dontrun{
#' Decoys <- allDecoys(version = '2024',
#'  ion.mode = 'Pos', max.tails = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1, 14.3, 15.3,
#'   16.5, 17.3, 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0", 
#'   Identified = Lipids)
#'   }

allDecoys <- function(Identified, version, max.tails, exact.tails, ion.mode="Pos"){
  
  ## Load references for class
  Identified <- Identified[Identified$FeatureType=="Lipid",]
  #Get index of short.name in identified lipids 
  identified_short_name <- grep("[sS]hort.*", colnames(Identified))
  lipid_class <- str_extract(Identified[1, identified_short_name], "[A-Z]*")
  
  if (!lipid_class %in% c("TAG","DAG", "PC", "PS", "PI", "PG","PA", "PE")){
    stop("Decoys only works for TAG, DAG, and PC.")
  }

  # references <- read_csv(paste0("inst/extdata/",lipid_class,"_", ion.mode, "_Core.csv")) 
  references <- getTAGs(tails = max.tails,exact_tails = exact.tails)
  #Get index of short.name reference list
  reference_short_name <- grep("[sS]hort.*", colnames(references))
  
  #Generates all possible FA tails to sample
  ref_max <- str_extract(references$Name, "[0-9]*:[0-9]*.*") %>% str_split(., "-") %>% unlist(., recursive = T) %>% unique(.)
  #Every unique identified lipid species
  total_identified <- unique(Identified[[identified_short_name]])

  ret <- c()
  
  for(o in 1:length(total_identified)){
    #Number of decoys to generate
    unique_tails <- sum(references[[reference_short_name]] %in% total_identified[o])
    
    #Sample with replacement with 20% extra because duplicates may be removed
    if (lipid_class=="TAG") {
      tails <- replicate(round(unique_tails*1.2,0), sample(ref_max, 3, replace = T))
    } else if (lipid_class=="DAG") {
      tails <- replicate(round(unique_tails*1.2,0), sample(ref_max, 2, replace = T))
    } else if (lipid_class=="PS" | lipid_class=="PC" |  lipid_class=="PE"| 
               lipid_class=="PG" | lipid_class=="PA" | lipid_class=="PI") {
      tails <- replicate(round(unique_tails*1.2,0), sample(ref_max, 2, replace = T))
    }
    
    #Converts sampled lipids to numbers, order smallest to largest, and revert 
    tails <- apply(tails, 2, function(x) {gsub("\\:", "\\.", x) %>% 
        as.numeric(.) %>% 
        .[order(.,decreasing = F)] %>% 
        sprintf('%.1f', .) %>%
        gsub("\\.", "\\:", .)}) 
    
    tails <- tails %>% t() %>% as.data.frame()
  
    duplicate_check <- do.call(paste, c(tails, sep = "-"))
    duplicate_check <- paste(lipid_class, duplicate_check) 
    
    #Check to see if any possible lipids were generated
    sample_tails <- Identified[Identified[[identified_short_name]]==total_identified[o],]
    is_duplicated <- which(duplicated(duplicate_check))
    is_duplicated <- c(is_duplicated, which(duplicate_check %in% sample_tails$Name))
    
    if (length(is_duplicated)!=0){
      duplicate_check <- duplicate_check[-is_duplicated]
    }
    
    #Search in reference to find fragments, formula, etc. 
    lipids_created <- data.frame(Name = duplicate_check)  
    lipids_created <- left_join(lipids_created, references, "Name")
    lipids_created[,grep("cursor", colnames(lipids_created))] <- references[references[[reference_short_name]] %in% total_identified[o],][,grep("cursor", colnames(lipids_created))][1,]
    lipids_created <- lipids_created[,match(colnames(references), colnames(lipids_created))]
    
    #For lipids with two FA tails, the references may not have lipids sorted in order smallest to largest.
    #We reverse the FA tails and search again. 
    not_matched <- lipids_created[is.na(lipids_created$Exact.mass),]
    not_matched_rev <- str_extract( not_matched$Name, "[0-9]*:[0-9]*.*") %>%
      str_split(., "-")  %>% 
      lapply(., function(x) {gsub("\\:", "\\.", x) %>% 
       as.numeric(.) %>% 
       .[order(.,decreasing = T)] %>% 
       sprintf('%.1f', .) %>%
       gsub("\\.", "\\:", .)})  %>%
      lapply(.,function(x) paste(x, collapse="-")) %>%
      unlist(.) %>% paste(lipid_class, .) %>% 
      data.frame(Name=.) %>% 
      inner_join(., references, "Name")
    
    #If still no matches, references is not comprehensive.
    if (nrow(not_matched)!=nrow(not_matched_rev)){
      stop("Is reference not comprehensive? Combination not found in reference.")
    }
    
    remove_na <- !is.na(lipids_created$Exact.mass)
    if(length(remove_na)!=0){
      lipids_created <- lipids_created[remove_na,]
      lipids_created <- rbind(lipids_created, not_matched_rev)
    }

    lipids_created$ID.Number <- total_identified[o]
    
    #Randomize order of generated lipids
    lipids_created <- lipids_created[sample(nrow(lipids_created)),]
    #Select the right number of lipids needed
    lipids_created <- lipids_created[1:unique_tails,] 
    
    ret <- rbind(ret, lipids_created)
  }
  dir.create("decoys", showWarnings = F)
  assign(paste0("Decoy_Data_",lipid_class,"_",version), ret, envir = .GlobalEnv)
  saveRDS(object = ret, file=paste0("./decoys/Decoy_Data_",lipid_class,"_",version))
  invisible(ret)
}
