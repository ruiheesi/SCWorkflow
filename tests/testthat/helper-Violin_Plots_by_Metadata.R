selectViolin <- function(dataset) {
  
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    object = selectCRObject("TEC")
    group.by = "orig_ident"
    group.subset =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
                               replace = FALSE)
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset") 
    object = selectCRObject("Chariou")
    group.by = "orig_ident"
    group.subset =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
                               replace = FALSE)
    
    
  } else if (dataset == "pbmc.single"){
    
    print("selected nsclc_single dataset") 
    object = selectCRObject("pbmc.single")
    group.by = "orig.ident"
    group.subset =  unique(object$orig.ident)
    genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
                               replace = FALSE)
    
    
  } else if (dataset == "nsclc.multi"){
    
    print("selected nsclc_multi dataset") 
    object = selectCRObject("nsclc.multi")
    group.by = "orig.ident"
    group.subset =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
                               replace = FALSE)
    
  } else if (dataset == "brca"){
    
    print("selected BRCA dataset") 
    object = selectCRObject("BRCA")
    group.by = "orig.ident"
    group.subset =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object$SCT@scale.data), 5, 
                               replace = FALSE)
    
  }
  
  return(list("object" = object, 
              "group.by" = group.by,
              "group.subset" = group.subset, 
              "genes.of.interest" = genes.of.interest))
  
}
