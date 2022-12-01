# This code comes from NIDAP 'Metadata Table [scRNA-seq][CCBR]' code template https://nidap.nih.gov/workspace/vector/templates/ri.vector.main.template.51ea15ce-2be0-415f-902b-3c86175eb6cd
# Documentation https://nidap.nih.gov/workspace/notepad/view/ri.notepad.main.notepad.53192624-dba0-4e63-9cbb-9bec776ede47
# Example https://nidap.nih.gov/workspace/vector/view/ri.vector.main.workbook.00e7a393-4cd7-4f20-bdd2-1c1c0cfd1ea2?branch=master


#' Extracting cell identities from the meta.data slot of Seurat-class object
#' 
#' Returns tabulation of cell identities present in the orig_ident column
#' 
#' It is recommended to run this function after doing SingleR annotations and before doing any visualizations/downstream analysis that requires inputting sample names.
#' 
#' @param SO Seurat-class object
#' 
#'   
#' @export
#' @return A table (data.frame) with the number of cells per each sample name (column names)


SampleNames <- function(SO) {
  # extract metadata table and sample names
  
  cat("\n Extracting Sample Names from metadata table: ")
  if ("meta.data" %in% slotNames(SO)) {
    if ("orig.ident" %in% colnames(SO@meta.data)) {
      dim.df <-
        as.data.frame(t(as.matrix(table(
          SO@meta.data$orig.ident
        ))))
    } else {
      dim.df <-
        as.data.frame(t(as.matrix(table(
          SO@meta.data$orig_ident
        ))))
    }
    colnames(dim.df) <-
      lapply(colnames(dim.df), function(x)
        gsub(".h5", "", x))
  } else {
    if ("orig.ident" %in% colnames(SO$RNA@meta.data)) {
      dim.df <-
        as.data.frame(t(as.matrix(
          table(SO$RNA@meta.data$orig.ident)
        )))
    } else {
      dim.df <-
        as.data.frame(t(as.matrix(
          table(SO$RNA@meta.data$orig_ident)
        )))
    }
    dim.df <-
      as.data.frame(t(as.matrix(table(
        SO$RNA@meta.data$orig.ident
      ))))
    colnames(dim.df) <-
      lapply(colnames(dim.df), function(x)
        gsub(".h5", "", x))
  }
  
  # output data.frame
  cat(sprintf("%g sample names returned", ncol(dim.df)))
  return(dim.df)
  
}