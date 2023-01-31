# This code comes from NIDAP 'Dotplot of Gene Expression by Metadata'

#' @title Dotplot of Gene Expression by Metadata
#' @description This function uses the Dotplot function from Seurat and plots 
#'        average gene expression values and percent expressed for a set of genes.  
#' @details This method provides a dotplot showing the percent frequency of 
#'        gene-positive cells as size of dot and degree of expression as color 
#'        of dot.
#' 
#' @param object Seurat Object
#' @param metadata Metadata category to plot
#' @param cells Metadata category factors to plot
#' @param markers Column that contains gene names
#' @param plot.reverse Set Metadata categories to x-axis (default is FALSE)
#' @param cell.reverse.sort Reverse Metadata Categories to plot (default is FALSE)
#' @param dot.color Dot color (default is "dark blue")

#' @importFrom tidyr pivot_wider
#' @importFrom Seurat subset Idents DotPlot

#' @export 
#' 
#' @return Dotplot with markers and cell types. 

DotplotMet <- function(object,
                       metadata, 
                       cells,
                       markers,
                       plot.reverse = FALSE,
                       cell.reverse.sort = FALSE,
                       dot.color = "darkblue") {
  
  #Make the metadata match:
  metadata.df <- object@meta.data
  Idents(object) <- metadata.df[[metadata]]
  
  #Error messages depending on the input celltype category
  
  #Calculate difference in input category elements from metadata column elements
  a <- sum(cells %in% unique(Idents(object)))
  b <- sum(!unique(Idents(object)) %in% cells)
  c <- sum(!cells %in% unique(Idents(object)))
  
  if(a < 2){
    stop("At least 2 metadata categories you wish to plot should 
         match the metadata column you are plotting.")
  }
  
  if(b>0){
    missinglab <- unique(Idents(object))[!unique(Idents(object)) %in% cells]
    warning(paste0("There are ",b," additional element(s) in the metadata table category that were missing from
                   your input categories:",
                   paste(as.character(missinglab),collapse="\n")))
  } 
  
  if(c>0){
    missinglab2 <- celltype[!celltype %in% unique(Idents(object))]
    warning(paste0("There are ",c," additional elements in your input categories 
                 that are missing from your metadata table: "))
    missinglab2 <- cat(paste(as.character(missinglab2),collapse="\n"))
  }
  
  #Subset object by identity
  object <- subset(object, idents = cells)
  
  #Bring in input genes and custom names
  markers <- gsub("[[:space:]]", "", markers)
  l1 <- length(markers)
  dups <- markers[duplicated(markers)]
  markers <- markers[!duplicated(markers)]
  
  l2 <- length(markers)
  print(paste("There are ",l2," total unique genes in the genelist"))
  if(l1 > l2){warning(paste0("\n\nThe following duplicate genes were removed: ",
                             dups))}
  missing.genes <- markers[!markers %in% rownames(object)]
  
  l3 <- length(missing.genes)
  missing.genes <- paste(shQuote(missing.genes), collapse=", ")
  if(l3 == l2){stop("No genes listed are found in dataset.")}
  if(l3 > 0){warning(paste0("\n\n",l3," genes are absent from dataset:", 
                      missing.genes,
              ". Possible reasons are that gene is not official gene symbol",
              " or gene is not highly expressed and has been filtered.\n "))}
  markers <- markers[markers %in% rownames(object)]
  cells <- cells[cells != ""]
  
  
  dp <- DotPlot(object, assay="SCT", features=markers, 
                dot.scale=4,
                cols = c("lightgrey", dot.color))
  cells <- cells[cells %in% dp$data$id]
  if(length(cells) < 1){
    stop("None of the cell types in list match cell types in metadata")
  }
  cells.missing <- cells[!cells %in% dp$data$id]
  if(length(cells.missing) > 0){
    cells.missing <- paste(shQuote(cells.missing), collapse=", ")
    warning(paste0("Some categories are missing from your dataset: ", cells.missing))
  }
  dp$data$id <- factor(dp$data$id, levels=rev(cells))
  dp$data$features.plot <- factor(dp$data$features.plot, levels=markers)
  
  plot <- ggplot(data = dp$data, 
                mapping = aes_string(x = "features.plot", y = "id")) + 
                geom_point(mapping = aes_string(size = "pct.exp", 
                                color = "avg.exp.scaled")) + 
                scale_color_gradient(low = "lightgrey", high = dot.color) +
                theme_cowplot() + 
                theme(axis.title.x = element_blank(), 
                      axis.text.x = element_text(angle = 90)) + 
                labs(y=metadata) 
  
  if(plot.reverse == TRUE){
        plot <- plot + coord_flip() 
  }
  if(cell.reverse.sort == TRUE){
        plot <- plot + scale_y_discrete(limits = rev(levels(dp$data$id)))
  }
  
  #Tabular format of Dotplot data are provided
  dp$data %>% select(features.plot,pct.exp,id) %>% tidyr::pivot_wider(names_from = features.plot, values_from = pct.exp) -> dp.pct.tab
  dp$data %>% select(features.plot,avg.exp.scaled,id) %>% tidyr::pivot_wider(names_from = features.plot, values_from = avg.exp.scaled) -> dp.exp.tab
  
  result.list <- list("plot" = plot, "pct" = dp.pct.tab, "exp" = dp.exp.tab)
  
  return(result.list)

}





