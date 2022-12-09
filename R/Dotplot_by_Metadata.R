# This code comes from NIDAP 'Dotplot of Gene Expression by Metadata [scRNA-Seq][CCBR] (79573b27-8a93-4f22-9863-993be1a44fc1): v16' code template

#' @title Dotplot of Gene Expression by Metadata
#' @description This Dotplot plotter plots average gene expression values for a set of genes from an input table.  Input table contains a single column for "Genes" and a single column for "Cell Type".  The values in the "Cell Type" column should match the values provided in the metadata template (Metadata Category to Plot).  The plot will order the genes (x-axis, left to right) and Cell Types (y-axis, top to bottom) in the order in which it appears in the input table.  Based on the additional column selected (Sample Column where sample or any other metadata category), it will display a contingency table for the Metadata 
#' @details This method provides a visualization plot showing the frequency of positively
#' 
#' @param input.dataset Seurat Object
#' @param metadata Metadata category to plot
#' @param sample Column that contains sample id 
#' @param cell.type Column that contains the cell type information
#' @param genes.column Column that contains gene names
#' @param dot.color Dot color (default is "dark red")

#' @import tidyverse
#' @import Seurat
#' @import cowplot

#' @export 
#' 
#' @return Dotplot with markers and cell types. 

DotplotMet <- function(object,
                       metadata, 
                       sample.column, 
                       cell.type,
                       markers, 
                       dot.color = "darkred") {
  
  #Make the metadata match:
  metadata.df <- object@meta.data
  #colnames(metadata.df) <- gsub("\\.","_",colnames(metadata.df))
  Idents(object) <- metadata.df[[metadata]]
  
  #Bring in input genes and custom names
  markers <- gsub("[[:space:]]", "", markers)
  markers <- markers[!duplicated(markers)]
  markers <- markers[markers %in% rownames(object)]
  cell.type <- cell.type[cell.type != ""]
  dp <- DotPlot(object, assay="SCT", features=markers, 
                dot.scale=4,
                cols = c("lightgrey", dot.color))
  dp
  dp$data$id <- factor(dp$data$id, levels=rev(cell.type))
  dp$data$features.plot <- factor(dp$data$features.plot, levels=markers)
  
  plot <- ggplot(data = dp$data, mapping = aes_string(x = "features.plot", y = "id")) + 
          geom_point(mapping = aes_string(size = "pct.exp", color = "avg.exp.scaled")) + 
          #scale_color_gradient2(midpoint=0, low="blue", mid="white",
          #                high="red", space ="Lab" ) +
          scale_color_gradient(low = "lightgrey", high = dot.color) +
          theme_cowplot() + 
          theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) + 
          labs(y=metadata)
  
  # Generate Contingency Table for Annotated cell types
  sample_column = sub("_",".",sample_column)
  cluster_num <- as.data.frame.matrix(table(object@meta.data[[sample_column]],object@meta.data[[metadata]]))
  cluster_num %>% rownames_to_column("Samples") -> cluster_num
  
  result.list <- list("plot" = plot, "data" = cluster_num)
  
  return(result.list)
}




