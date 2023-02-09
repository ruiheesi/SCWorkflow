select_dataset_SCviolin <- function(dataset) {
  
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    object <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
    ident.of.interest = "orig_ident"
    groups.of.interest =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object), 5, replace = FALSE)
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset") 
    object <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
    ident.of.interest = "orig_ident"
    groups.of.interest =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object), 5, replace = FALSE)
    
    
  } else if (dataset == "NSCLC_Single"){
    
    print("selected NSCLC_Single dataset") 
    object <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
    ident.of.interest = "orig.ident"
    groups.of.interest =  unique(object$orig.ident)
    genes.of.interest = sample(rownames(object), 5, replace = FALSE)
    
    
  } else if (dataset == "NSCLC_Multi"){
    
    print("selected NSCLC_Multi dataset") 
    object <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
    ident.of.interest = "orig.ident"
    groups.of.interest =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object), 5, replace = FALSE)
    
  } else if (dataset == "BRCA"){
    
    print("selected BRCA dataset") 
    object <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    ident.of.interest = "orig.ident"
    groups.of.interest =  unique(object$orig.ident)[1:3]
    genes.of.interest = sample(rownames(object$SCT@scale.data), 5, replace = FALSE)
    
  }
  
  return(list("object" = object, "ident.of.interest" = ident.of.interest,
              "groups.of.interest" = groups.of.interest, "genes.of.interest" = genes.of.interest))
  
}
