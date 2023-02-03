
#This helper script will return parameters for each test dependent on the data input:

getparam_Cbg_auto <- function(data){

  if(data == "TEC"){
    
      object = readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
      samples_to_include = unique(object$orig.ident)
      samples_to_display = unique(object$orig.ident)
      marker_list = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      cells_of_interest = colnames(marker_list)[1:3]
      
  } else if (data == "Chariou") {
    
      object = readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
      samples_to_include = unique(object$orig.ident)
      samples_to_display = unique(object$orig.ident)
      marker_list = read.csv(test_path("fixtures", "Marker_Table_demo.csv"))
      cells_of_interest = colnames(marker_list)[4:6]
      
  } else if (data == "NSCLC_Single") {
    
    object = readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
    samples_to_include = unique(object$orig.ident)
    samples_to_display = unique(object$orig.ident)
    marker_list = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                             rand_type2 = sample(rownames(object), 5, replace = FALSE),
                             rand_type3 = sample(rownames(object), 5, replace = FALSE))
    cells_of_interest = colnames(marker_list)
    
  } else if (data == "NSCLC_Multi") {

    object = readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
    samples_to_include = unique(object$orig.ident)
    samples_to_display = unique(object$orig.ident)
    marker_list = data.frame(rand_type1 = sample(rownames(object), 5, replace = FALSE),
                             rand_type2 = sample(rownames(object), 5, replace = FALSE),
                             rand_type3 = sample(rownames(object), 5, replace = FALSE))
    cells_of_interest = colnames(marker_list)

  } else if (data == "BRCA") {

    object = readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    samples_to_include = unique(object$orig.ident)
    samples_to_display = unique(object$orig.ident)
    marker_list = data.frame(rand_type1 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE),
                             rand_type2 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE),
                             rand_type3 = sample(rownames(object$SCT@scale.data), 5, replace = FALSE))
    cells_of_interest = colnames(marker_list)
    
  }
  
  return(list("object" = object, "samples_to_include"= samples_to_include, 
              "samples_to_display" = samples_to_display, 
              "marker_list" = marker_list, 
              "cells_of_interest" = cells_of_interest))  
}
