#' Returns tabulation of cell identities present in the orig_ident column
#' 
#' It is recommended to run this function after doing SingleR annotations and before doing any visualizations/downstream analysis that requires inputting sample names.
#' 
#' @param SO Seurat-class object
#' 
#' @importFrom methods slotNames
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