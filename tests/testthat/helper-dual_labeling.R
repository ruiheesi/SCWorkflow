#This helper script will return parameters for each test dependent on the data input:
getparamdl <- function(data){

  if(data == "TEC"){
      object <- select_crobject("TEC") 
      samples <- c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")
      set.seed(15)
      marker1 <- sample(rownames(object),1, replace = FALSE)
      marker.1.type = "SCT"
      set.seed(1)
      marker2 <- sample(rownames(object),1, replace = FALSE)
      marker.2.type = "SCT"
      data.reduction = "tsne"
      density.heatmap = TRUE
      display.unscaled.values = TRUE
  } else if (data == "Chariou") {
      object <- select_crobject("Chariou")
      samples <- c("PBS","CD8dep","ENT","NHSIL12","Combo")
      varfeat <- VariableFeatures(object@assays$SCT)[1:100]
      set.seed(254)
      marker1 <- sample(varfeat,1)
      marker.1.type = "SCT"
      set.seed(35)
      marker2 <- sample(varfeat,1)
      marker.2.type = "SCT"
      data.reduction = "tsne"
      density.heatmap = TRUE
      display.unscaled.values = TRUE
  }
  return(list("object" = object,
              "samples" = samples,
              "marker1" = marker1,
              "marker.1.type" = marker.1.type,
              "marker2" <- marker2,
              "marker.2.type" = marker.2.type,
              "data.reduction" = data.reduction,
              "density.heatmap" = density.heatmap,
              "display.unscaled.values" = display.unscaled.values))  
}