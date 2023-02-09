
#This helper script will return parameters for each test dependent on the data input:

getparam_Cbg_auto <- function(data){

  if(data == "TEC"){
    
      object = readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
      samples.to.include = unique(object$orig.ident)
      samples.to.display = unique(object$orig.ident)
      marker.list = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      cells.of.interest = colnames(marker.list)[1:3]
      
  } else if (data == "Chariou") {
    
      object = readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
      samples.to.include = unique(object$orig.ident)
      samples.to.display = unique(object$orig.ident)
      marker.list = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      cells.of.interest = colnames(marker.list)[4:6]
      
  } else if (data == "NSCLC_Single") {
    
    object = readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
    samples.to.include = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    marker.list = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                             rand_type2 = sample(rownames(object), 5, replace = FALSE),
                             rand_type3 = sample(rownames(object), 5, replace = FALSE))
    cells.of.interest = colnames(marker.list)
    
  } else if (data == "NSCLC_Multi") {

    object = readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
    samples.to.include = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    marker.list = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                             rand_type2 = sample(rownames(object), 5, replace = FALSE),
                             rand_type3 = sample(rownames(object), 5, replace = FALSE))
    cells.of.interest = colnames(marker.list)

  } else if (data == "BRCA") {

    object = readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    samples.to.include = unique(object$orig.ident)
    samples.to.display = unique(object$orig.ident)
    marker.list = data.frame(rand_type1 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE),
                             rand_type2 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE),
                             rand_type3 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE))
    cells.of.interest = colnames(marker.list)
    
  }
  
  return(list("object" = object, "samples.to.include"= samples.to.include, 
              "samples.to.display" = samples.to.display, 
              "marker.list" = marker.list, 
              "cells.of.interest" = cells.of.interest))  
}
