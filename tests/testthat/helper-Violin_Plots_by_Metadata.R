select_dataset_SCviolin <- function(dataset) {
  
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    object <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
    ident_of_interest = "orig_ident"
    groups_of_interest =  unique(object$orig.ident)[1:3]
    genes_of_interest = sample(rownames(object), 5, replace = FALSE)
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset") 
    object <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
    ident_of_interest = "orig_ident"
    groups_of_interest =  unique(object$orig.ident)[1:3]
    genes_of_interest = sample(rownames(object), 5, replace = FALSE)
    
    
  } else if (dataset == "NSCLC_Single"){
    
    print("selected NSCLC_Single dataset") 
    object <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
    ident_of_interest = "orig.ident"
    groups_of_interest =  unique(object$orig.ident)
    genes_of_interest = sample(rownames(object), 5, replace = FALSE)
    
    
  } else if (dataset == "NSCLC_Multi"){
    
    print("selected NSCLC_Multi dataset") 
    object <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
    ident_of_interest = "orig.ident"
    groups_of_interest =  unique(object$orig.ident)[1:3]
    genes_of_interest = sample(rownames(object), 5, replace = FALSE)
    
  } else if (dataset == "BRCA"){
    
    print("selected BRCA dataset") 
    object <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    ident_of_interest = "orig.ident"
    groups_of_interest =  unique(object$orig.ident)[1:3]
    genes_of_interest = sample(rownames(object$SCT@scale.data), 5, replace = FALSE)
    
  }
  
  return(list("object" = object, "ident_of_interest" = ident_of_interest,
              "groups_of_interest" = groups_of_interest, "genes_of_interest" = genes_of_interest))
  
}
