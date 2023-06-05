getPseudobulkParam <- function(data){
  
  if(data == "TEC"){
    
    object = select_crobject("TEC")
    contrasts = c("F-M")
    replicate = 'orig.ident'
    subgroup = 'Gender'
    group = "Gender"
    
  } else if (data == "Chariou") {
    
    object = select_crobject("Chariou")
    contrasts = c("F-M")
    replicate = 'orig.ident'
    subgroup = 'Gender'
    group = "Gender"
    
  } else if (data == "pbmc.single") {
    
    object = select_crobject("nsclc-single")
    object@meta.data$group_cluster = paste("cluster", 
                                           object@meta.data$SCT_snn_res.0.2, 
                                           sep = "_")
    contrasts = c("cluster_1-cluster_2")
    replicate = 'Phase'
    subgroup = 'group_cluster'
    group = "group_cluster"
    
  } else if (data == "nsclc.multi") {
    
    object = select_crobject("nsclc-multi")
    object@meta.data$group_cluster = paste("cluster", 
                                           object@meta.data$SCT_snn_res.0.2, 
                                           sep = "_")
    contrasts = c("cluster_1-cluster_2")
    replicate = 'Phase'
    subgroup = 'group_cluster'
    group = "group_cluster"
    
  } else if (data == "BRCA") {
    
    object = select_crobject("BRCA")
    contrasts = c("G2M-S")
    replicate = 'orig.ident'
    subgroup = 'Phase'
    group = "Phase"
    
  }
  
  return(list("object" = object, 
              "contrasts"= contrasts, 
              "replicate" = replicate,
              "subgroup" = subgroup,
              "group" = group))  
}
