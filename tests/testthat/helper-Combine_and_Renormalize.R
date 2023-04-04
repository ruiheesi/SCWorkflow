getParamCR <- function(data){
  
  if(data == "TEC"){
    object <- readRDS(
      test_path(paste0("fixtures/",data), 
                paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    npcs = 30
    vars.to.regress = c()
    clust.res.low=0.2
    clust.res.high = 1.2
    only.var.genes = FALSE 
    
    
  } else if (data == "Chariou") {
    object <- readRDS(
      test_path(paste0("fixtures/",data), 
                paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    npcs = 21
    vars.to.regress = c("percent.mt")
    clust.res.low=2.4
    clust.res.high = 2.8
    only.var.genes = FALSE 
    
    
  } else if (data == "pbmc-single") {
    object <- readRDS(
      test_path(paste0("fixtures/",data), 
                paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    npcs = 30
    vars.to.regress = c()
    clust.res.low=0.2
    clust.res.high = 1.2
    only.var.genes = FALSE 
    
    
  } else if (data == "nsclc-multi") {
    object <- readRDS(
      test_path(paste0("fixtures/",data), 
                paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    npcs = 30
    vars.to.regress = c()
    clust.res.low=0.2
    clust.res.high = 1.2
    only.var.genes = FALSE 
    
    
  } else if (data == "BRCA") {
    object <- readRDS(
      test_path(paste0("fixtures/",data), 
                paste0(data,"_Filtered_PCA_Norm_SO_downsample.rds")))
    npcs = 30
    vars.to.regress = c()
    clust.res.low=0.2
    clust.res.high = 1.2
    only.var.genes = TRUE 
    
  }
  
  return(list("object" = object,
              "npcs" = npcs,
              "vars.to.regress" = vars.to.regress,
              "clust.res.low"=clust.res.low,
              "clust.res.high" = clust.res.high,
              "only.var.genes" = only.var.genes ))  
}