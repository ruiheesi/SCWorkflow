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
  l1 <- length(markers)
  dups <- markers[duplicated(markers)]
  markers <- markers[!duplicated(markers)]
  
  l2 <- length(markers)
  if(l1 > l2){message(paste0("\n\nThe following duplicate genes removed: ",dups))}
  missing.genes <- markers[!markers %in% rownames(object)]
  
  l3 <- length(missing.genes)
  missing.genes <- paste(shQuote(missing.genes), collapse=", ")
  if(l3 == l2){stop("No genes listed are found in dataset.")}
  if(l3 > 0){message(paste0("\n\n",l3," genes are absent from dataset:", 
                          missing.genes,
                          ". Possible reasons are that gene is not official gene symbol",
                          " or gene is not highly expressed and has been filtered.\n\n "))}
  markers <- markers[markers %in% rownames(object)]
  cell.type <- cell.type[cell.type != ""]
  cell.type <- cell.type[cell.type %in% dp$data$id]
  if(length(cell.type) < 1){
    error("None of the cell types in list match cell types in metadata")
  }
  cell.type.missing <- cell.type[!cell.type %in% dp$data$id]
  if(length(cell.type.missing) < length(cell.type)){
    message("")
  }
  
  dp <- DotPlot(object, assay="SCT", features=markers, 
                dot.scale=4,
                cols = c("lightgrey", dot.color))
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
  
  
  fs <- input_dataset$fileSystem()
  path <- fs$get_path("seurat_object.rds", 'r')
  so <- readRDS(path)
  
  metadata.df <- so@meta.data
  colnames(metadata.df) <- gsub("\\.","_",colnames(metadata.df))
  Idents(so) <- metadata.df[[metadata_category_to_plot]]
  
  #Bring in input category:
  celltype <- category_labels_and_genes_table[[category_labels_to_plot]][!is.na(category_labels_and_genes_table[[category_labels_to_plot]])]
  celltype <- celltype[celltype != "null"]
  
  #Calculate difference in input category elements from metadata column elements
  a <- length(unique(Idents(so)))
  b <- length(celltype)
  c <- a - b
  d <- b - a
  
  if(sum(celltype %in% Idents(so)) == 0){
    stop("ERROR: The category from the input table you wish to plot should match the metadata column you are plotting.")
  } else if(c>0){
    missinglab <- unique(Idents(so))[!unique(Idents(so)) %in% celltype]
    print(paste0("There are ",c," elements in the metadata table category that are missing from your input table: "))
    missinglab <- cat(paste(as.character(missinglab),collapse="\n"))
    cat("\n")
  } else if(d>0){
    missinglab2 <- celltype[!celltype %in% unique(Idents(so))]
    print(paste0("There are ",d," elements in the category from your input table that are missing from your metadata table: "))
    missinglab2 <- cat(paste(as.character(missinglab2),collapse="\n"))
    cat("\n")
  }
  
  #Bring in input gene names:
  markers <- category_labels_and_genes_table[[genes_column]][!is.na(category_labels_and_genes_table[[genes_column]])]
  markers <- markers[markers != "null"]
  
  #Deal with missing genes:
  dupmarkers <- markers[duplicated(markers)]
  if(length(dupmarkers) > 0){
    dupmarkers <- paste(as.character(dupmarkers),collapse="\n")
    print(paste0("There are duplicated markers:"))
    dupmarkers <- cat(dupmarkers)
    cat("\n")
  }
  
  markers <- markers[!duplicated(markers)]
  missingmarkers <- markers[!markers %in% rownames(so)]
  e <- length(missingmarkers)
  
  if(e>0){
    print(paste0("There are ",e," markers that are missing from your single cell dataset: "))
    missingmarkers <- cat(paste(as.character(missingmarkers),collapse="\n"))
  }
  
  #Draw Dotplot using Seurat function 
  dp <- DotPlot(so, assay="SCT", features=markers, dot.scale=4,cols = c("lightgrey", dot_color))
  dp$data$id <- factor(dp$data$id, levels=rev(celltype))
  dp$data$features.plot <- factor(dp$data$features.plot, levels=markers)
  dp$data <- na.omit(dp$data) #remove na's which show up in plot
  
  #Plot Bubbleplot using Seurat dotplot ggplot object:
  plot <- ggplot(data = dp$data, mapping = aes_string(x = "features.plot", y = "id")) + 
    geom_point(mapping = aes_string(size = "pct.exp", color = "avg.exp.scaled")) + 
    scale.func(range = c(0, 4)) + 
    theme_cowplot() + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) + labs(y=metadata_category_to_plot)
  print(plot)
  
  # Show numeric data matching image:
  if(show_percentage_cells_expressing == TRUE){
    dp$data %>% select(features.plot,pct.exp,id) %>% pivot_wider(names_from = features.plot, values_from = pct.exp) -> dp.pct.tab
    return(dp.pct.tab)
  } else {
    dp$data %>% select(features.plot,avg.exp,id) %>% pivot_wider(names_from = features.plot, values_from = avg.exp) -> dp.exp.tab
    return(dp.exp.tab)
  }
}

}




