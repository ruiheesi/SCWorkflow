getPseudobulkParam <- function(data){
  
  if(data == "TEC"){
    
    object = selectCRObject("TEC")
    contrasts = c("F-M")
    replicate = 'orig.ident'
    subgroup = 'Gender'
    group = "Gender"
    
  } else if (data == "Chariou") {
    
    object = selectCRObject("Chariou")
    contrasts = c("F-M")
    replicate = 'orig.ident'
    subgroup = 'Gender'
    group = "Gender"
    
  } else if (data == "pbmc.single") {
    
    object = selectCRObject("pbmc.single")
    object@meta.data$group_cluster = paste("cluster", 
                                           object@meta.data$SCT_snn_res.0.2, 
                                           sep = "_")
    contrasts = c("cluster_1-cluster_2")
    replicate = 'Phase'
    subgroup = 'group_cluster'
    group = "group_cluster"
    
  } else if (data == "nsclc.multi") {
    
    object = selectCRObject("nsclc.multi")
    object@meta.data$group_cluster = paste("cluster", 
                                           object@meta.data$SCT_snn_res.0.2, 
                                           sep = "_")
    contrasts = c("cluster_1-cluster_2")
    replicate = 'Phase'
    subgroup = 'group_cluster'
    group = "group_cluster"
    
  } else if (data == "BRCA") {
    
    object = selectCRObject("BRCA")
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
