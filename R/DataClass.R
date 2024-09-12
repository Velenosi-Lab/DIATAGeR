### DataClass
setOldClass(c("tbl_df", "tbl", "data.frame"))
DIAData<-setClass("DIAData",
                  representation(Results = "data.frame",
                                 SampleType = "data.frame",
                                 RawSpectra = "list",
                                 ResultsNestedSpectra = "data.frame",
                                 DeisotopeSpectra = "list",
                                 Identified = "list",
                                 Decoys = "list", 
                                 AlignResultForProgen = "data.frame"),

                  prototype(Results = data.frame(),
                    SampleType = data.frame(),
                            RawSpectra = list(),
                            DeisotopeSpectra = list(),
                            ResultsNestedSpectra = data.frame(),
                            Identified = list(),
                            Decoys = list()
                           ))

setGeneric("addDecoyClass", function(object, ...)
  standardGeneric("addDecoyClass"))
setMethod("addDecoyClass", "DIAData",
          function(object,Lipid,value){
            object@Decoys[[Lipid]]<-value
            if (validObject(object))
              return(object)
          })

setGeneric("addLipidClass", function(object, ...)
  standardGeneric("addLipidClass"))
setMethod("addLipidClass", "DIAData",
          function(object,Lipid,value){
            object@Identified[[Lipid]]<-value
            if (validObject(object))
              return(object)
          })

setGeneric("addLipidClassMultiCore", function(object, ...)
  standardGeneric("addLipidClassMultiCore"))
setMethod("addLipidClassMultiCore", "DIAData",
          function(object,value){
            object@Identified<-value
            if (validObject(object))
              return(object)
          })


LipidMode<-setClass("LipidMode",
                    representation (
                      Pos = "DIAData",
                      Neg = "DIAData",
                      IdentifiedCombined = "list",
                      FeaturesCombined = "data.frame"
                    ),)

#NO
setGeneric("MergeID", function(object, ...)
  standardGeneric("MergeID"))
setMethod("MergeID", "LipidMode",
          function(object){
            PosIdentified<-object@Pos@Identified
            NegIdentified<-object@Neg@Identified
            if (validObject(object))
              return(PosIdentified)
            return(NegIdentified)
          })
#NO
setGeneric("IdentifiedCombined", function(object, ...)
  standardGeneric("IdentifiedCombined"))
setMethod("IdentifiedCombined", "LipidMode",
          function(object,value){
            object@IdentifiedCombined<-value
            if (validObject(object))
              return(object)

          })

setGeneric("addAlignResultForProgen", function(object, ...)
  standardGeneric("addAlignResultForProgen"))
setMethod("addAlignResultForProgen", "DIAData",
          function(object,value){
            object@AlignResultForProgen<-value
            if (validObject(object))
              return(object)

          })
#NO
setGeneric("CombineModeClass", function(PosObject, ...)
  standardGeneric("CombineModeClass"))
setMethod("CombineModeClass", "DIAData",
          function(PosObject=PosObject,NegObject=NegObject){
            if(length(grep("+",PosObject@Results$`Adduct type`))
               <length(grep("+",NegObject@Results$`Adduct type`)))
              stop("Positive and Negative objects are reversed")

            object<-new("LipidMode")
            object@Pos<-PosObject
            object@Neg<-NegObject


            if (validObject(object))
              return(object)
          })

setGeneric("addRawSpectra", function(object, ...)
  standardGeneric("addRawSpectra"))
setMethod("addRawSpectra", "DIAData",
          function(object,value){
            object@RawSpectra<-value
            if (validObject(object))
              return(object)
          })

setGeneric("addDeisotopeSpectra", function(object, ...)
  standardGeneric("addDeisotopeSpectra"))
setMethod("addDeisotopeSpectra", "DIAData",
          function(object,value){
            object@DeisotopeSpectra<-value
            if (validObject(object))
              return(object)
          })

setGeneric("addResultsNestedSpectra", function(object, ...)
  standardGeneric("addResultsNestedSpectra"))
setMethod("addResultsNestedSpectra", "DIAData",
          function(object,value){
            object@ResultsNestedSpectra<-value
            if (validObject(object))
              return(object)
          })

####### cbind.fill function
# cbind.fill<-function(...,fill=NULL)
# {
#   inputs<-list(...)
#   inputs<-lapply(inputs,vert)
#   maxlength<-max(unlist(lapply(inputs,len)))
#   bufferedInputs<-lapply(inputs,buffer,length.out=maxlength,fill,preserveClass=FALSE)
#   return(Reduce(cbind.data.frame,bufferedInputs))
# }
# 
# vert<-function(object)
# {
#   #result<-as.data.frame(cbind(as.matrix(object)))
#   if(is.list(object))
#     object<-cbind(object)
#   object<-data.frame(object)
# 
#   return(object)
# }
# 
# len <- function(data)
# {
#   result<-ifelse(is.null(nrow(data)),length(data),nrow(data))
#   return(result)
# }
# 
# buffer<-function(x,length.out=len(x),fill=NULL,preserveClass=TRUE)
# {
#   xclass<-class(x)
#   input<-lapply(vert(x),unlist)
#   results<-as.data.frame(lapply(input,rep,length.out=length.out))
#   if(length.out>len(x) && !is.null(fill))
#   {
#     results<-t(results)
#     results[(length(unlist(x))+1):length(unlist(results))]<-fill
#     results<-t(results)
#   }
#   if(preserveClass)
#     results<-as2(results,xclass)
#   return(results)
# }


