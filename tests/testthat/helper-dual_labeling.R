#This helper script will return parameters for each test dependent on the data input:
getParamDL <- function(data){

  if(data == "TEC"){
      object <- selectCRObject("TEC") 
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
      object <- selectCRObject("Chariou")
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
  } else if (data == "pbmc-single") {
      object <- selectCRObject("pbmc-single")
      samples <- c("PBMC_Single")
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
  } else if (data == "nsclc-multi") {
      object <- selectCRObject("nsclc-multi")
      samples <- c("Donor_1","Donor_2","Donor_3")
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
  } else if (data == "BRCA") {
      object <- selectCRObject("BRCA")
      samples <- c("CID4471","CID4290A","CID44971","CID4040","CID4513","CID4535")
      varfeat <- VariableFeatures(object@assays$SCT)
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