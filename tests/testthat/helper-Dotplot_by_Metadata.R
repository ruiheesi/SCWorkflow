
#This helper script will return parameters for each test dependent on the data input:

getparam <- function(data){

  if(data == "TEC"){
      object <- select_crobject("TEC")
      metadata_column <- "SCT_snn_res.0.2"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c(0,1,2,3,4,5,6,7,8,9)
      plot.reverse <- FALSE
      cell.reverse.sort <- FALSE
  } else if (data == "Chariou") {
      object <- select_crobject("Chariou")
      metadata_column <- "orig.ident"
      set.seed(15)
      markers <- sample(rownames(object),10, replace = FALSE)
      cells <- c("PBS","CD8dep","NHSIL12","ENT","Combo")
      plot.reverse <- TRUE
      cell.reverse.sort <- TRUE
  }
  return(list("object" = object, "metadata"=metadata_column, 
              "markers" = markers, "cells" = cells, 
              "plot.reverse" = plot.reverse,
              "cell.reverse.sort" = cell.reverse.sort))  
}