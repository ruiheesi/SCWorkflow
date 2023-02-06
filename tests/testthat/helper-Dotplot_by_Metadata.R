
#This helper script will return parameters for each test dependent on the data input:

getparamdp <- function(data){

  if(data == "TEC"){
      object <- select_crobject("TEC")
      metadata <- "SCT_snn_res.0.2"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c(0,1,2,3,4,5,6,7,8,9)
      plot.reverse <- FALSE
      cell.reverse.sort <- FALSE
  } else if (data == "Chariou") {
      object <- select_crobject("Chariou")
      metadata <- "orig.ident"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("PBS","CD8dep","NHSIL12","ENT","Combo")
      plot.reverse <- TRUE
      cell.reverse.sort <- TRUE
  } else if (data == "nsclc-single") {
      object <- select_srobject("nsclc-single")
      metadata <- "BP_encode_main"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("B-cells","CD4+ T-cells","CD8+ T-cells","DC",
                 "HSC","Monocytes","NK cells")  
      plot.reverse <- TRUE
      cell.reverse.sort <- TRUE
  } else if (data == "nsclc-multi") {
      object <- select_srobject("nsclc-multi")
      metadata <- "BP_encode_main"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("B-cells","CD4+ T-cells","CD8+ T-cells","DC",
                 "HSC","Monocytes","NK cells")  
      plot.reverse <- TRUE
      cell.reverse.sort <- TRUE
  } else if (data == "BRCA") {
      object <- select_crobject("BRCA")
      metadata <- "SCT_snn_res.0.2"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("0","1","2","3","4","5","6")  
      plot.reverse <- TRUE
      cell.reverse.sort <- TRUE
  }
  
  
  
  
  return(list("object" = object, "metadata"=metadata, 
              "markers" = markers, "cells" = cells, 
              "plot.reverse" = plot.reverse,
              "cell.reverse.sort" = cell.reverse.sort))  
}