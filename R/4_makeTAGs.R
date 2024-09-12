

getFAs <- function(Max,exact_tails=NULL) {
  if (!is.null(exact_tails)) {
    return(exact_tails)
  }
  
  FAs <- c()
  for (i in 1:length(Max)){
    FAs <- c(FAs, seq(floor(Max[i]), Max[i], by = 0.1))
  }
  return(FAs)
}

getCarbons <- function(x) {
  y = lapply(x, function(x) lapply(x, function(x) x[[1]] ))
  lapply(y, function(x) sum(as.numeric(unlist(x)))+3)
}

getDoubleBonds <- function(x) {
  y = lapply(x, function(x) lapply(x, function(x) x[[2]] ))
  lapply(y, function(x) sum(as.numeric(unlist(x))))
}

getHydrogens <- function(Carbons, Doubles){
  (Carbons-3)*2+2-Doubles*2
}

getFormula <- function(Carbons, Hydrogens, Oxygen = 6){
  paste0("C", Carbons, "H", Hydrogens, "O", Oxygen)
}

getMass <- function(Carbons, Hydrogens, Oxygen = 6){
  12.000000*Carbons+1.007825*Hydrogens+15.994915*Oxygen
}

doNothing <- function(x) x

###Enter max double bonds here:
#Enter the desired number of carbons and the maximum number of double bonds for each fatty acyl chain length.
#Can be seperated by space or with a comma. For example, either of the following: 
#tails = "8.2 9.0 10.2 11.0 12.3 13.1 14.3 15.3 16.5 17.3 18.5 19.5 20.6 21.5 22.6 23.0 24.4 25.0 26.0"
#tails = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1, 14.3, 15.3, 16.5, 17.3, 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0"

#' Create triacylglycerols library 
#'
#' @param tails list of desired number of carbons and the maximum number of double 
#' bonds for each fatty acyl chain length. Format: number of carbons. number of
#' double bonds (Eg: tails = 8.2 creates three fatty acyl chains: 8:0, 8:1 and 8:2)
#' @param exact_tails library can be customizable for specific tails 
#' (Eg: exact_tails = c("18.1","18.2) creates library of TGs made up 18:1 and 18:2 fatty acids)
#'
#' @return library of triacylglycerols with information of precursors and fragments masses
#' @export
#'
#' @examples
#' \dontrun{
#' Library <- getTAGs(tails = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1, 14.3, 15.3, 
#' 16.5, 17.3, 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0")
#' }
#' 

getTAGs <- function(tails = "8.2, 9.0, 10.2, 11.0, 12.3, 13.1,14.3, 15.3, 16.5, 17.3, 18.5, 19.5, 20.6, 21.5, 22.6, 23.0, 24.4, 25.0, 26.0",exact_tails=NULL){
  
  if (grepl(",", tails)) { 
    Max = read.delim(sep = ",", header = F, text = tails) %>% t()
  } else {
    Max <- read.delim(sep = " ", header = F, text = tails) %>% t()
  }  

  FAs <- getFAs(Max,exact_tails)
  
  Tails <- combinations(length(FAs), 3, FAs, repeats.allowed = T)
  if (is.null(exact_tails)) {
    Tails <- apply(Tails, 2, function(x) sprintf('%.1f', x) )
  }  
  Formula <- as.data.frame(t(Tails))
  Formula <- lapply(Formula, doNothing)
  Formula <- lapply(Formula, str_split, "\\.")
  
  CarbonsNumber <- getCarbons(Formula) %>% as.numeric()
  DoubleBondsNumber <- getDoubleBonds(Formula) %>% as.numeric()
  HydrogensNumber <- getHydrogens(CarbonsNumber,DoubleBondsNumber)
  MolecFormula <- getFormula(CarbonsNumber, HydrogensNumber)
  ExactMass <- getMass(CarbonsNumber, HydrogensNumber)
  
  Name <- do.call(paste, list(Tails[,1], Tails[,2], Tails[,3], sep = "-")) 
  Name <- sapply(Name, function(x) gsub("\\.", "\\:", x))
  Name <- paste("TAG", Name) %>% as.data.frame()
  ShortName <- paste(CarbonsNumber-3, DoubleBondsNumber, sep = ":")
  ShortName <- paste("TAG", ShortName) %>% as.data.frame()
  
  Chain_1_length <- lapply(Formula, function(x) x[[1]][[1]]) %>% unlist(.) %>% unname(.) %>% as.numeric(.)
  Chain_1_bonds <-  lapply(Formula, function(x) x[[1]][[2]]) %>% unlist(.) %>% unname(.) %>% as.numeric(.)
  Chain_2_length <- lapply(Formula, function(x) x[[2]][[1]]) %>% unlist(.) %>% unname(.) %>% as.numeric(.)
  Chain_2_bonds <- lapply(Formula, function(x) x[[2]][[2]]) %>% unlist(.) %>% unname(.) %>% as.numeric(.)
  Chain_3_length <- lapply(Formula, function(x) x[[3]][[1]]) %>% unlist(.) %>% unname(.) %>% as.numeric(.)
  Chain_3_bonds <- lapply(Formula, function(x) x[[3]][[2]]) %>% unlist(.) %>% unname(.) %>% as.numeric(.)
  
  
  
  Precursor_1 <- ExactMass+18.03437
  Precursor_2 <- ExactMass+22.98977
  Precursor_3 <- ExactMass+1.00782
  Temp <- ExactMass-(Chain_1_length*12) - (((2*(Chain_1_length-1)+1)-(2*Chain_1_bonds))*1.00782) - (2*15.99492) - (1.00782)
  Fragment_1 <- Temp+1.00782
  Temp <- ExactMass-(Chain_2_length*12) - (((2*(Chain_2_length-1)+1)-(2*Chain_2_bonds))*1.00782) - (2*15.99492) - (1.00782)
  Fragment_2 <- Temp+1.00782
  Temp <- ExactMass-(Chain_3_length*12) - (((2*(Chain_3_length-1)+1)-(2*Chain_3_bonds))*1.00782) - (2*15.99492) - (1.00782)
  Fragment_3 <- Temp+1.00782
  
  All <- data.frame(ID.Number=1:length(ExactMass), 
               logP=rep(999, length(ExactMass)),
               Formula=MolecFormula, 
               Short.name=unlist(ShortName),
               Name=unlist(Name), 
               First_Precursor_mz=Precursor_1,
               Second_Precursor_mz=Precursor_2,
               Third_Precursor_mz=Precursor_3,
               First_Fragment_mz=Fragment_1, 
               Second_Fragment_mz=Fragment_2,
               Third_Fragment_mz=Fragment_3) 
  rownames(All) <- 1:nrow(All)
  
  return(All)
}


