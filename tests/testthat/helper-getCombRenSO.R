select_crobject <- function(dataset) {
  
  if (dataset == "TEC"){
    
    print("selected TEC dataset") 
    inputObject <- readRDS(test_path("fixtures/TEC", 
        "TEC_Combine_and_Renormalize_SO_downsample.rds"))
    
  } else if (dataset == "Chariou"){
    
    print("selected Chariou dataset")
    inputObject <- readRDS(test_path("fixtures/Chariou", 
        "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
    
  } 
  
  return(inputObject)
  
}