
#This helper script will return parameters for each test dependent on the data input:

getparam_Pseudobulk <- function(data){

  if(data == "TEC"){
    
      object = readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
      contrasts = c("F-M")
      replicate = 'orig.ident'
      comparison_level = 'Gender'
      label = "Gender"

  } else if (data == "Chariou") {
    
      object = readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
      contrasts = c("F-M")
      replicate = 'orig.ident'
      comparison_level = 'Gender'
      label = "Gender"
      
  } else if (data == "NSCLC_Single") {
    
      object = readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
      object@meta.data$label_cluster = paste("cluster", object@meta.data$SCT_snn_res.0.2, sep = "_")
      contrasts = c("cluster_1-cluster_2")
      replicate = 'Phase'
      comparison_level = 'label_cluster'
      label = "label_cluster"
    
  } else if (data == "NSCLC_Multi") {

      object = readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
      object@meta.data$label_cluster = paste("cluster", object@meta.data$SCT_snn_res.0.2, sep = "_")
      contrasts = c("cluster_1-cluster_2")
      replicate = 'Phase'
      comparison_level = 'label_cluster'
      label = "label_cluster"

  } else if (data == "BRCA") {

      object = readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
      contrasts = c("G2M-S")
      replicate = 'orig.ident'
      comparison_level = 'Phase'
      label = "Phase"
    
  }
  
  return(list("object" = object, "contrasts"= contrasts, 
              "replicate" = replicate,
              "comparison_level" = comparison_level,
              "label" = label))  
}
