getParamPN <- function(data){

  if(data == "TEC"){
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    vars.to.regress = c("percent.mt")
    vars.to.plot = c('percent.mt')
    
  } else if (data == "Chariou") {
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    vars.to.regress = c("percent.mt")
    vars.to.plot = c('percent.mt')

  } else if (data == "pbmc-single") {
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    vars.to.regress = c("percent.mt")
    vars.to.plot = c('percent.mt')
      
  } else if (data == "nsclc-multi") {
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    vars.to.regress = c("percent.mt")
    vars.to.plot = c('percent.mt')
    
  } else if (data == "BRCA") {
    object=readRDS(test_path(
      paste0("fixtures/",data), paste0(data,"_Filtered_SO_downsample.rds")))
    vars.to.regress = c("percent.mt")
    vars.to.plot = c('percent.mt')
    
  }
  
  return(list("object" = object,
              "npcs" = npcs,
              "vars.to.regress" = vars.to.regress,
              "vars.to.plot"=vars.to.plot))    
}