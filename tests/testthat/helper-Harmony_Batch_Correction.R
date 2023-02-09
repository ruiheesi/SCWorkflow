
#This helper script will return parameters for each test dependent on the data input:

getparam_Harmony <- function(data){

  if(data == "TEC"){
    
      object = readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.vars = c("orig.ident")

  } else if (data == "Chariou") {
    
      object = readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.vars = c("orig.ident")
      
  } else if (data == "NSCLC_Single") {
    
      object = readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.vars = c("Phase")
    
  } else if (data == "NSCLC_Multi") {

      object = readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.vars = c("orig.ident")

  } else if (data == "BRCA") {

      object = readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
      genes.to.add = sample(rownames(object), 5, replace = FALSE)
      group.by.vars = c("orig.ident")
    
  }
  
  return(list("object" = object, "genes.to.add"= genes.to.add, 
              "group.by.vars" = group.by.vars))  
}
