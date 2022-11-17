# This code comes from NIDAP 'Metadata Table [scRNA-seq][CCBR]' code template
# Template Manager https://nidap.nih.gov/workspace/vector/templates/ri.vector.main.template.6a8139b7-45b4-4c6a-8648-c8ec34e6fc60
# Documentation https://nidap.nih.gov/workspace/notepad/view/ri.notepad.main.notepad.70a5a28a-fa78-4f77-bee2-30f2d4854304
# Example https://nidap.nih.gov/workspace/vector/view/ri.vector.main.workbook.d9fbfd7b-fb08-4814-ae16-7327cb7f4f37?branch=master


#' @title Extracting meta.data and reductions slots from Seurat Object
#' @description Returns the per cell metadata and dimension embeddings (optional) of the Seurat Object as a data table.
#' @details Exposes the metadata of any combined Seurat object, though it is recommended to input a SingleR annotated Seurat Object.
#' 
#' @param SO object of class Seurat
#' @param return.cell.embeddings If TRUE, cell embeddings from dimensional reductions present in the Seurat Object are added to the output metadata table
#' 
#' @importFrom Seurat AddMetaData
#' @importFrom methods slotNames
#' @importFrom tibble rownames_to_column
#'   
#' @export
#' 
#' @seealso
#' [SampleNames()] [Seurat::AddMetaData()]
#' 
#' @examples
#' \dontrun{ 
#' # Return metadata table
#' metadata <- MetadataTable(SO = seurat.object)
#' head(metadata)
#' 
#' # Return dimensional reductions
#' metadata <- MetadataTable(SO = seurat.object, return.cell.embeddings = TRUE)
#' umap_coordinates <- metadata %>% dplyr::select(Barcode, seurat_clusters, contains("UMAP"))
#' 
#' ## Plot 
#' library(ggplot2)
#' ggplot(umap_coordinates, aes(x = UMAP_1, y = UMAP_2, colour = seurat_clusters)) + 
#' geom_point(size = 0.1) + theme_bw()
#' }
#' 
#' @return a data.frame extracted from the Seurat object


MetadataTable <- function(SO, return.cell.embeddings = FALSE) {
  # get cell embeddings? ====
  
  if (return.cell.embeddings) {
    reds = names(SO@reductions)
    cat(sprintf(
      "\n Extracting cell embeddings from Seurat Object: %s\n\n",
      paste(reds, collapse = ", ")
    ))
    for (i in seq_along(reds)) {
      SO = AddMetaData(SO, as.data.frame(SO@reductions[[i]]@cell.embeddings))
    }
  } else {
    cat("\n Skipping extracting cell embeddings from Seurat Object\n\n")
  }
  
  # extract meta.data ====
  
  cat(" Extracting Metadata Table from Seurat Object\n\n")
  if ("meta.data" %in% slotNames(SO)) {
    met.df <- SO@meta.data
  } else {
    met.df <- SO$RNA@meta.data
  }
  
  # barcode column ====
  
  if (!("Barcode" %in% colnames(met.df))) {
    met.df %>% rownames_to_column("Barcode") -> met.df
  }
  
  # return output ====
  
  #print(summary(met.df))
  return(met.df)
}