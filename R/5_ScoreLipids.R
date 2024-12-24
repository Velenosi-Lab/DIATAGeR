
##### Generates Scores #####
getScoreNew <- function(TAGsAll, standard.file= NULL, rt.type=c("truncated", "exact"), rttol=2){

  Unscored <- TAGsAll[TAGsAll$Correlation.tan==-1,]
  TAGsAll <- TAGsAll[TAGsAll$Correlation.tan!=-1,]
  TAGsAll$Rev.dot.product[is.na(TAGsAll$Rev.dot.product)] <- median(TAGsAll$Rev.dot.product, na.rm=T)
  # TAGsAll$dot_product_original_single[is.na(TAGsAll$dot_product_original_single)] <- median(TAGsAll$dot_product_original_single, na.rm=T)


  
  # Caret functions to build model
  inTrain <- createDataPartition(TAGsAll$Class, p = 0.70, list = FALSE)
  
  TAGTRAIN <- TAGsAll[inTrain, ]
  TAGTEST <- TAGsAll[-inTrain, ]

  #Removes standards from training data if exists
  #sa <- standardsCheck(TAGTRAIN, standard.file= standard.file, rt.type=rt.type, rttol=rttol)
  #TAGTRAIN <- TAGTRAIN[-which((TAGTRAIN$ID.Number %in% sa$ID.Number) & (TAGTRAIN$rt %in% sa$rt)), ]

   if (all(TAGTRAIN$Correlation.tan == -1)) {
    stop("No correlations available.")
   }
  

  
  # Customize the train method in trainControl
  train_model <- trainControl(method = "cv", number = 5, 
                              summaryFunction = twoClassSummary, classProbs = TRUE,
                              savePredictions = "all")
  model1 <- suppressWarnings({
    train(
    Class ~   Rev.dot.product + sqrt(Missing.fragments)  + log10(FragIntensity) + Correlation.tan.within+ Correlation.tan,
# rev_dot_product + sqrt(Missing.fragments) + Correlation.tan + log10(FragIntensity),
    data = TAGTRAIN,
    method = "glm",
    family = "binomial",
    trControl = train_model
  )
  })

  # Predict scores for all lipids.
  Score_1 <- predict(model1, newdata = TAGsAll, type = "prob")$Tp
  
  # Prints AUC and creates AUC object as a metric for the user.
  lipid.pred <- predict(model1, newdata = TAGTEST, type = "prob")
  myroc1 <- roc(TAGTEST$Class, lipid.pred$Tp)
  auc1 <- round(auc(myroc1), 3)
  
  TAGsAll$Score <- Score_1
  Unscored$Score <- 0
  assign('AUC', value = auc1, envir=.GlobalEnv)
  
  TAGsAll <- rbind(TAGsAll,Unscored)
}

##### Generates Scores #####
getScoreDotProd <- function(TAGsAll, standard.file= NULL, rt.type=c("truncated", "exact"), rttol=2){

  Unscored <- TAGsAll[TAGsAll$Correlation.tan==-1,]
  TAGsAll <- TAGsAll[TAGsAll$Correlation.tan!=-1,]
  # TAGsAll$dot_product_original[is.na(TAGsAll$dot_product_original)] <- median(TAGsAll$dot_product_original[TAGsAll$Class=="Tp"], na.rm=T)
  # TAGsAll$dot_product_original_single[is.na(TAGsAll$dot_product_original_single)] <- median(TAGsAll$dot_product_original_single[TAGsAll$Class=="Tp"], na.rm=T)
  
  # Caret functions to build model
  inTrain <- createDataPartition(TAGsAll$Class, p = 0.70, list = FALSE)
  
  TAGTRAIN <- TAGsAll[inTrain, ]
  TAGTEST <- TAGsAll[-inTrain, ]
  
  #Removes standards from training data if exists
  #sa <- standardsCheck(TAGTRAIN, standard.file= standard.file, rt.type=rt.type, rttol=rttol)
  #TAGTRAIN <- TAGTRAIN[-which((TAGTRAIN$ID.Number %in% sa$ID.Number) & (TAGTRAIN$rt %in% sa$rt)), ]
  
  if (all(TAGTRAIN$Correlation.tan == -1)) {
    stop("No correlations available.")
  }
  
  # Customize the train method in trainControl
  train_model <- trainControl(method = "cv", number = 5, 
                              summaryFunction = twoClassSummary, classProbs = TRUE,
                              savePredictions = "all")
  model1 <- train(
    Class ~   rev_dot_product+ dot_product + Tails,
    
    data = TAGTRAIN,
    method = "glm",
    family = "binomial",
    trControl = train_model
  )
  
  #+ sqrt(Missing.fragments) + Correlation.tan + log10(FragIntensity)+
  # Predict scores for all lipids.
  Score_1 <- predict(model1, newdata = TAGsAll, type = "prob")$Tp
  
  # Prints AUC and creates AUC object as a metric for the user.
  lipid.pred <- predict(model1, newdata = TAGTEST, type = "prob")
  myroc1 <- roc(TAGTEST$Class, lipid.pred$Tp)
  auc1 <- round(auc(myroc1), 3)
  
  TAGsAll$DotProdScore <- Score_1
  Unscored$DotProdScore <- 0
  assign('AUC', value = auc1, envir=.GlobalEnv)
  
  TAGsAll <- rbind(TAGsAll,Unscored)
}




##### Generates Scores #####
getScore <- function(TAGsAll){
  set.seed(123)
  Unscored <- TAGsAll[TAGsAll$Correlation==-1,]
  TAGsAll <- TAGsAll[TAGsAll$Correlation!=-1,]
  TAGsAll$dot_product[is.na(TAGsAll$dot_product)] <- median(TAGsAll$dot_product[TAGsAll$Class=="Tp"], na.rm=T)
  
  
  #Caret functions to build model
  inTrain<- createDataPartition(TAGsAll$Class,p=0.700, list=FALSE)
  TAGTRAIN <- TAGsAll[ inTrain,]
  TAGTEST  <- TAGsAll[-inTrain,]
  
  if(all(TAGTRAIN$Correlation.tan==-1)){
    stop("No correlations available.")
  }
  
  train_model <- trainControl(method = "cv", number = 10,
                              summaryFunction = twoClassSummary,  classProbs = T, 
                              savePredictions = "all")
  
  model1 <- train(
    Class ~ Missing.fragments+ Correlation.tan + log10(FragIntensity)+ rev_dot_product+ Entropy+Entropy_similarity,
    
    data = TAGTRAIN,
    method = "glm",
    family = "binomial",
    trControl = train_model,
    metric="ROC")
  print(model1$finalModel)
  #Predict scores for all lipids.
  Score_1 <- predict(model1, newdata = TAGsAll, type="prob")$Tp
  
  #Prints AUC and creates AUC object as a metric for the user. 
  lipid.pred = predict(model1, newdata = TAGTEST, type="prob")
  myroc1 <- roc(TAGTEST$Class, lipid.pred$Tp)
  auc1 <- round(auc(myroc1),3)
  assign('AUC', value = auc1, envir=.GlobalEnv)
  print(paste("Area under the ROC Curve:", paste0(auc1, "."), "Any number over 0.6 is good."))
  
  TAGsAll$Score <- Score_1
  Unscored$Score <- 0
  
  TAGsAll <- rbind(TAGsAll,Unscored )
  
  return(TAGsAll)
}

#' Generate Scores for ranking lipids. 
#'
#' @param DIADataObj The object stores alignment file and spectra
#' @param format tell the function if the alignment file is MSDIAL or Generic Format
#' @param ion.mode Mass spectrometry ion mode, the same as in LipidIdentifier. 
#' @param lipid Lipid class, the same as in LipidIdentifier. 
#' @param version Version name, the same as in LipidIdentifier. Also used to pull identified lipids from the ion.mode folder. 
#' @param ppmtol PPM tolerance, the same as in LipidIdentifier. Default 15. 
#' @param rttol Retention time tolerance, the same as in LipidIdentifier. Default 5. 
#' @param spectra.file File path for retrieving spectra from all samples. 
#' @param ppmtolMS1 Defaults to ppmtol if not specified. 
#' @param run.code for calculating the correlation between fragments and precursor intensities. 
#' This helps recognize the columns of samples in the alignment file. Default "Sample"
#' @param ms1.precursors Number of precursors found in MS1 to identify compounds. Default to 1 
#' @param ms2.precursors Number of precursor found in MS2 to identify compounds. 
#' @param spectra.file.type specify txt, msp or mgf files. The user will have determined the spectra name
#' @param max.tails list of desired number of carbons and the maximum number of double 
#' bonds for each fatty acyl chain length. 
#' Default = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1, 14.3, 15.3, 16.5, 17.3, 18.5, 
#' 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0"
#' @param ... Additional parameters
#' @param exact.tails  library can be customizable for specific tails
#'
#' @return scores for targeted and decoy triacylglycerols
#' @export
#'
#' @examples
#' \dontrun{
#' ResultsTAG <- ScoreLipids(DIADataObj = DIA_Pos,
#' lipid = "TAG", 
#' ion.mode = "Pos",
#' version = "24", 
#' ppmtol= 15, 
#' rttol = 5, 
#' ms1.precursors =1, ms2.precursors =1, 
#' spectra.file.type ="txt",
#' max.tails = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1, 14.3, 15.3, 16.5, 17.3, 
#' 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0",
#' spectra.file = "Data/Centroid/")
#' }
#' 



ScoreLipids <-  function(DIADataObj, 
                         ion.mode,format = c("MSDIAL","Generic","Progenesis"), 
                         lipid, version, ppmtol= 15, ppmtolMS1, rttol = 5, 
                         run.code ="Sample",
                         ms1.precursors= c(1), ms2.precursors = c(0,1),
                         spectra.file.type=c("txt","msp","mgf","mgf_mzmine"),
                         max.tails= "8.2, 9.0, 10.2, 11.0, 12.3, 13.1,14.3, 15.3, 16.5, 17.3, 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0",
                         exact.tails = NULL,
                         spectra.file, ...){
  
  if (missing(spectra.file)){
    stop("Please provide path of the spectra files.")
  }

  files <- grep(pattern = paste0("/Identified_", lipid, "_", ion.mode, "_", version, "_"), list.files(ion.mode, full.names = T))
  
  if(length(files)>1){
    stop(paste(paste("\n",list.files(ion.mode, full.names = T)[files], collapse = ""), "exist. \n Please remove duplicates with same version!" ))
  } 
  
  if(length(files)==0){
    stop(paste(paste0("/Identified_TAG_Pos_", version, "_#"), "does not exist."))
  }
    
  Identified <- readRDS(list.files(ion.mode, full.names = T)[files])
  
  allDecoys(version = version, ion.mode = ion.mode, max.tails = max.tails, exact.tails = exact.tails, Identified = Identified)
  
  DecoyIdentifier(DIADataObj, lipid = lipid, ion.mode = ion.mode, format = format, ppmtol = ppmtol,
                  ppmtolMS1 = ppmtolMS1 , rttol = rttol,
                  ms1.precursors=ms1.precursors, ms2.precursors=ms2.precursors,
                  version = version, decoy.name = version, 
                  max.tails = max.tails,
                  write.annotations = T)
  
  #Loads the lipids identified from the decoy reference data. 
  filesDecoy <- grep(pattern = paste0("/IdentifiedDecoys_", lipid, "_", ion.mode, "_", version, "_"), list.files(ion.mode, full.names = T))
  
  if(length(filesDecoy)>1){
    stop(paste(paste("\n",list.files(ion.mode, full.names = T)[filesDecoy], collapse = ""), "exist. \n Please remove duplicates with same version!" ))
  } 
  
  IdentifiedDecoy <- readRDS(paste0(ion.mode, "/", paste0("/IdentifiedDecoys_", lipid, "_", ion.mode, "_", version, "_1")))
  if (is.null(IdentifiedDecoy)){
    stop("No decoys were identified in the data!")
  }
  #We have Identified, which is the positive spectra matches, and IdentifiedDecoy, which are the decoy spectra matches. 
  #We combined them, replace the ID.Number with the lipid species name, the name of the original species whose spectra 
  #was used to identify decoy lipids. 
  Identified$Class <- "Tp"
  Identified <- Identified[Identified$FeatureType=="Lipid",]
  IdentifiedDecoy$Class <- "Fp"
  IdentifiedDecoy <- IdentifiedDecoy[IdentifiedDecoy$FeatureType=="Lipid",]
  AllIdentified <- rbind(Identified,IdentifiedDecoy)
  AllIdentified[grep("^[0-9]", AllIdentified$ID.Number),]$ID.Number <- AllIdentified[grep("^[0-9]", AllIdentified$ID.Number),]$Short.name
  
  AllIdentified <- retrieveSpectra(Lipids=AllIdentified, ppmtol = ppmtol, spectra.file = spectra.file,spectra.file.type=spectra.file.type)
  AllIdentified <- calcRevDotProduct(Lipids=AllIdentified)
  AllIdentified <- calcCorrelation(DIADataObj = DIADataObj, run.code = run.code, Lipids=AllIdentified)
  AllIdentified <- calcCorrelationWithin( Lipids=AllIdentified,ms2.precursors=ms2.precursors)
  AllIdentified <- getScoreNew(AllIdentified)
  
  ##11th Dec: Edit format of the final results 
  AllIdentified<-AllIdentified[,-c(1,2,8)]
  AllIdentified<-AllIdentified[,c("ID.simple","Formula","Short.name","Name","rt","First_Precursor_mz","mz","MSMSref","FeatureType","Class","Rev.dot.product","Correlation.tan","Correlation.tan.within","Missing.fragments","FragIntensity","Score")]
  
  colnames(AllIdentified)<-c("Feature.ID","Formula","Species","Molecular.species","Retention.time","Precursor.mz.ref","Precursor.mz","MSMSref",
                             "FeatureType","Class","Rev.dot.product","Corr.precusor.fragments","Corr.fragments.fragments","Missing.fragments","FragIntensity","Score")
  return(AllIdentified)
}
