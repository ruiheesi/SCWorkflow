select_crobject <- function(dataset) {
  
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    inputObject <- readRDS(test_path("fixtures/TEC", 
        "TEC_Combine_and_Renormalize_SO_downsample.rds"))
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset")
    inputObject <- readRDS(test_path("fixtures/Chariou", 
        "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
    
  } else if (dataset == "nsclc-single"){
    
    print("selected NSCLC Single dataset")
    inputObject <- readRDS(test_path("fixtures/NSCLC_Single",
    "NSCLCsingle_Cell_Types_SingleR_SO_downsample.rds"))
    
  } else if (dataset == "nsclc-multi"){
    
    print("selected NSCLC Multi dataset")
    inputObject <- readRDS(test_path("fixtures/NSCLC_Multi",
    "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  
  } else if (dataset == "BRCA"){
    
    print("selected BRCA dataset")
    inputObject <- readRDS(test_path("fixtures/BRCA",
        "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
  }
  return(inputObject) 
  
}


select_srobject <- function(dataset){
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    inputObject <- readRDS(test_path("fixtures/TEC", 
          "TEC_CellTypesSingleR_SO_downsample.rds"))
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset")
    inputObject <- readRDS(test_path("fixtures/Chariou", 
          "Chariou_Cell_Types_SingleR_SO_downsample.rds"))
    
  } else if (dataset == "nsclc-single"){
    
    print("selected NSCLCsingle dataset")
    inputObject <- readRDS(test_path("fixtures/NSCLC_Single", 
          "NSCLCsingle_Cell_Types_SingleR_SO_downsample.rds"))
  
  } else if (dataset == "nsclc-multi"){
    
    print("selected NSCLCmulti dataset")
    inputObject <- readRDS(test_path("fixtures/NSCLC_Multi", 
       "NSCLCmulti_Cell_Types_SingleR_SO_downsample.rds"))
  
  } 
}
  
