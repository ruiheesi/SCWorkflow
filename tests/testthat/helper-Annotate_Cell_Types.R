#This helper script will return parameters for each test dependent on the input:

getParamACT <- function(data){

  if(data == "TEC"){
      object <- selectCRObject("TEC")
      species <- "Mouse"
  } else if (data == "Chariou") {
      object <- selectCRObject("Chariou")
      species <- "Mouse"
  } else if (data == "pbmc-single") {
      object <- selectCRObject("pbmc-single")
      species <- "Human"
  } else if (data == "nsclc-multi") {
      object <- selectCRObject("nsclc-multi")
      species <- "Human"
  } else if (data == "BRCA") {
      object <- selectCRObject("BRCA")
      species <- "Human"
  }
 
  return(list("object" = object, "species" = species))  
}
