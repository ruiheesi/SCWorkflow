#This helper script will return parameters for each test dependent on the data 
#input ("TEC","Chariou","pbmc-single","nsclc-multi","BRCA"):
getParamDL <- function(data){

  if(data == "TEC"){
      object <- selectCRObject("TEC") 
      samples <- c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")
      set.seed(15)
      var.feat <- VariableFeatures(object@assays$SCT)[1:100]
      marker.1 <- sample(var.feat,1)
      marker.1.type = "SCT"
      set.seed(1)
      marker.2 <- sample(var.feat,1)
      marker.2.type = "SCT"
      data.reduction = "tsne"
      density.heatmap = TRUE
      display.unscaled.values = TRUE
  } else if (data == "Chariou") {
      object <- selectCRObject("Chariou")
      samples <- c("PBS","CD8dep","ENT","NHSIL12","Combo")
      var.feat <- VariableFeatures(object@assays$SCT)[1:100]
      set.seed(254)
      marker.1 <- sample(var.feat,1)
      marker.1.type = "SCT"
      set.seed(35)
      marker.2 <- sample(var.feat,1)
      marker.2.type = "SCT"
      data.reduction = "tsne"
      density.heatmap = TRUE
      display.unscaled.values = TRUE
  } else if (data == "pbmc-single") {
      object <- selectCRObject("pbmc-single")
      samples <- c("PBMC_Single")
      var.feat <- VariableFeatures(object@assays$SCT)[1:100]
      set.seed(254)
      marker.1 <- sample(var.feat,1)
      marker.1.type = "SCT"
      set.seed(35)
      marker.2 <- sample(var.feat,1)
      marker.2.type = "SCT"
      data.reduction = "tsne"
      density.heatmap = TRUE
      display.unscaled.values = TRUE
  } else if (data == "nsclc-multi") {
      object <- selectCRObject("nsclc-multi")
      samples <- c("Donor_1","Donor_2","Donor_3")
      var.feat <- VariableFeatures(object@assays$SCT)[1:100]
      set.seed(254)
      marker.1 <- sample(var.feat,1)
      marker.1.type = "SCT"
      set.seed(35)
      marker.2 <- sample(var.feat,1)
      marker.2.type = "SCT"
      data.reduction = "tsne"
      density.heatmap = TRUE
      display.unscaled.values = TRUE
  } else if (data == "BRCA") {
      object <- selectCRObject("BRCA")
      samples <- c("CID4471","CID4290A","CID44971","CID4040","CID4513","CID4535")
      var.feat <- VariableFeatures(object@assays$SCT)
      set.seed(254)
      marker.1 <- sample(var.feat,1)
      marker.1.type = "SCT"
      set.seed(35)
      marker.2 <- sample(var.feat,1)
      marker.2.type = "SCT"
      data.reduction = "tsne"
      density.heatmap = TRUE
      display.unscaled.values = TRUE
  }
  
  return(list("object" = object,
              "samples" = samples,
              "marker.1" = marker.1,
              "marker.1.type" = marker.1.type,
              "marker.2" <- marker.2,
              "marker.2.type" = marker.2.type,
              "data.reduction" = data.reduction,
              "density.heatmap" = density.heatmap,
              "display.unscaled.values" = display.unscaled.values))  

}  

.drawdualplot <- function(x, width = 10, height = 10){
#.drawdualplot <- function(x){
    path <- tempfile(fileext = ".png")
    png(path,
        width=width,
        height=height,
        units = "in",
        res = 400
    )
    on.exit(dev.off())
    grid.draw(x)
    path
}
