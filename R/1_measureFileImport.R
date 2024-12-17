####
#'Import the feature list (MSDIAL alignment file or Generic Format)
#'
#' @param filename MSDIAL alignment file or Generic Import Format (.txt file)
#' @param ResultsFiletype tell the function if the import file is MSDIAL Format or Generic Format 
#' @param IonMode Pos or Neg to indicate positive or negative ion mode
#' @param AlignmentFileName NULL
#'
#' @return Imported feature list
#' @export
#'
#' @examples
#' \dontrun{
#' DIA_Pos <-measureFileImport(filename = "Height_0_20239151316.txt",
#' IonMode = "Pos", ResultsFiletype = "MSDIAL")
#' }

measureFileImport <-
  function(filename, ResultsFiletype = c("MSDIAL","Generic"), IonMode = c("Pos","Neg"), AlignmentFileName = NULL) {
    
    if (ResultsFiletype == "Generic") {
      Importfile_Generic <- fread(filename)
      # Add IDs and edit name of columns 
      Importfile_Generic$ID <- 0:(nrow(Importfile_Generic)-1)
      names(Importfile_Generic)[names(Importfile_Generic) == "RT [min]"] <- "rt"
      names(Importfile_Generic)[names(Importfile_Generic) == "m/z"] <- "mz"
      Importfile_Generic$rt<- as.numeric(Importfile_Generic$rt)
      Importfile_Generic$mz<-as.numeric(Importfile_Generic$mz)
      Importfile_Generic<- Importfile_Generic %>% select (ID,everything())
      
      # Importfile_Generic <- Importfile_Generic %>%
      #   rename_with(~ str_replace(., "^Area", "Sample"))
      
      # The reference spectrum is sample with highest area: 
      sample_columns_names<- names(Importfile_Generic)[4:(ncol(Importfile_Generic))]
      Importfile_Generic[, 'Spectrum reference file name' := sample_columns_names[max.col(.SD, ties.method = "first")], .SDcols = sample_columns_names]
      Importfile_Generic<- Importfile_Generic %>% select (ID,rt,mz,'Spectrum reference file name',everything())
      
      # Put "Sample" at the begining of the name of sample column: 
      sample_columns <- 5:(ncol(Importfile_Generic))
      colnames(Importfile_Generic)[sample_columns] <- paste0("Sample: ", colnames(Importfile_Generic)[sample_columns])
      }
    
    
    
    
    if (ResultsFiletype == "MSDIAL"|ResultsFiletype =="Progenesis") {
      if (ResultsFiletype == "Progenesis" & is.null(AlignmentFileName)) {
        stop("Please add AlignmentFileName")
      }
      if (IonMode == "Pos") {
        hydrogen <- 1.0078250170
      } else if (IonMode == "Neg") {
        hydrogen <- -1.0078250170
      } else {
        stop("Specify `IonMode`: Pos or Neg")
      }
      
      if (is.null(AlignmentFileName)) {
        measureFileSubStart <- fread(filename)
      } else {
        measureFileSubStart <-  fread(filename)
      }
      
      if (is.na(measureFileSubStart[1, 1])|measureFileSubStart[1, 1]=="") {
        ### Extracting MS-DIAL data
        measureFileSub1 <- measureFileSubStart %>%
          slice(., c(which(!is.na(.[, 1])& .[, 1]!="") [1]:n()) ) %>%
          setNames(., unlist(.[1, ])) %>%
          .[-1, ] %>%
          type_convert(.)
        ### extracting MSDIAL SampleTypeInfo
        measureFileSub2 <- measureFileSubStart %>%
          slice(., c(1:(which(!is.na(.[, 1])& .[, 1]!="")[1]))) %>%
          select(., !starts_with("X"))
        measureFileSub2[4, 1] <- "SampleID"
      } else {
        measureFileSub1 <- measureFileSubStart
      }
      
      measureFileSub1$`Alignment ID` <- as.integer(measureFileSub1$`Alignment ID`)
      names(measureFileSub1)[names(measureFileSub1) == "Alignment ID"] <- "ID"
      names(measureFileSub1)[names(measureFileSub1) == "Average Rt(min)"] <- "rt"
      names(measureFileSub1)[names(measureFileSub1) == "Average Mz"] <- "mz"
      names(measureFileSub1)[names(measureFileSub1) == "Average Mz"] <- "mz"
      
      # Put "Sample" at the begining of the name of sample column: 
      sample_columns <- 33:(ncol(measureFileSub1)-2)
      
      colnames(measureFileSub1)[sample_columns] <- paste0("Sample: ", colnames(measureFileSub1)[sample_columns])
      
      
      # load adducts; Add_DF returns: character (empty)
      #Add_DF <- list.files("inst/adducts", full.names = T, recursive = T)
      Add_DF <- list.files(system.file("adducts", package="DIATAGeR"), full.names = T)
      
      Add_DF <- read_csv(Add_DF[grep(IonMode, Add_DF)], col_types = cols())
      Add_sub <- measureFileSub1 %>%
        select(ID, `Adduct type`) %>%
        rename(Adduct = `Adduct type`) %>%
        left_join(., Add_DF, "Adduct")
      
      # Test All adducts in feature list exist in adduct mass spreadsheet
      AdductTest <- Add_sub %>%
        filter(is.na(`Adduct Mass`)) %>%
        select(Adduct)
      # Features missing adducts are converted to M-H+ or M-H- 
      if (dim(AdductTest)[1] > 0) {
        # warning(c("Missing adduct masses for ", paste(levels(as.factor(AdductTest$Adduct)))))
        if (IonMode == "Pos") {
          Add_sub$'Adduct Mass'[is.na(Add_sub$'Adduct Mass')] <- 1.007276
        } else if (IonMode == "Neg") {
          Add_sub$'Adduct Mass'[is.na(Add_sub$'Adduct Mass')] <- -1.007276
        }
      }
      
      Add_sub <- Add_sub$`Adduct Mass`
      keep <- which(!colnames(measureFileSub1) %>% duplicated())
      measureFileSub1 <- measureFileSub1 %>% select(all_of(keep))
      # Add adduct mass to measureFileSub1
      measureFileSub1 <- add_column(measureFileSub1, `Adduct Mass` = Add_sub, .after = 5)
      
      # measureFileSub1[is.na(measureFileSub1$`Adduct type`),]
      
      # Calculate neutral mass using adduct mass
      measureFileSub1 <- measureFileSub1 %>%
        mutate(Neutral = mz - `Adduct Mass`) %>%
        select(1:3, Neutral, everything())
    }
    if (ResultsFiletype == "Progenesis") {
      measureFileSubStartProgen <- suppressWarnings(read_csv(filename, col_types = cols()))
      
      if (is.na(measureFileSubStartProgen[1, 1])) {
        measureFileSub1Progen <- measureFileSubStartProgen %>%
          slice(., c(which(!is.na(.[, 1]))[1]:n())) %>%
          setNames(., .[1, ]) %>%
          .[-1, ] %>%
          type_convert(.)
        
        measureFileSub2Progen <- measureFileSubStartProgen %>%
          slice(., c(1:(which(!is.na(.[, 1]))[1]))) %>%
          select(., c((which(!is.na(.[1, ]))[1]):ncol(.))) %>%
          t(.) %>%
          as_tibble(., .name_repair = "minimal") %>%
          setNames(., c("Group", "SampleID")) %>%
          fill(Group, .direction = "down") %>%
          select("SampleID", "Group")
      } else {
        measureFileSub1Progen <- measureFileSubStartProgen
      }
      
      measureFileSub1Progen <- measureFileSub1Progen %>%
        rename(Neutral = `Neutral mass (Da)`, rt = `Retention time (min)`, mz = `m/z`) %>%
        add_column(., ID = c(1:nrow(.)), .before = 1) %>%
        mutate(`Adduct Mass` = Neutral - mz, .after = 5) %>%
        select(ID, rt, mz, Neutral, everything())
      
      WithNeutral <- measureFileSub1Progen %>%
        filter(!is.na(Neutral))
      WithOutNeutral <- measureFileSub1Progen %>%
        filter(is.na(Neutral)) %>%
        mutate(Neutral = mz - hydrogen, `Adduct Mass` = hydrogen)
      
      measureFileSub1Progen <- bind_rows(WithNeutral, WithOutNeutral) %>%
        arrange(ID)
    }
    
    object <- new("DIAData")
    if (ResultsFiletype == "Generic") {
      object@Results <- Importfile_Generic
      if (exists("SampleType_Generic")) {
        object@SampleType <- SampleType_Generic
      }
    }
    if (ResultsFiletype == "MSDIAL") {
      object@Results <- measureFileSub1
      if (exists("measureFileSub2")) {
        object@SampleType <- measureFileSub2
      }
    }
    if (ResultsFiletype == "Progenesis") {
      object@Results <- measureFileSub1Progen
      object@AlignResultForProgen <- measureFileSub1
      if (exists("measureFileSub2Progen")) {
        object@SampleType <- measureFileSub2Progen
      }
    }
    # object@ID <- measureFileSub1$`Alignment ID`
    # object@rt <- measureFileSub1$`Average Rt(min)`
    # object@mz <- measureFileSub1$`Average Mz`
    # object@name <- measureFileSub1$`Metabolite name`
    # object@adduct <- measureFileSub1$`Adduct type`
    # object@SpectrumRef <- measureFileSub1$`Spectrum reference file name`
    # object@MS1IsotopeSpectrum <- measureFileSub1$`MS1 isotopic spectrum`
    
    return(object)
  }
