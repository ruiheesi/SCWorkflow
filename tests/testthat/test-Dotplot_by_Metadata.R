
test_that("dotplot produced and contingency table returned", {
  library(Seurat)
  obj <- readRDS("/rstudio-files/ccbr-data/users/maggie/SCWorkflow/tests/testthat/fixtures/SO_moduleScore.rds")
  celltype_table <- read.csv2("/rstudio-files/ccbr-data/users/maggie/SCWorkflow/tests/testthat/fixtures/cell_types_genes_table.tsv", header = TRUE, sep = "\t")
  metadata_column <- "Likely_CellType"
  sample_column <- "orig.ident"
  genes_column <- "Gene_Names"
  cell_type_column <- "Cluster_Names"
  
  markers <- celltype_table[[genes_column]][!is.na(celltype_table[[genes_column]])]
  celltype <- celltype_table[[cell_type_column]][!is.na(celltype_table[[cell_type_column]])]
  
  results.list <- DotplotMet(object = obj,
                             metadata = metadata_column, 
                             sample.column = sample_column, 
                             cell.type = celltype,
                             markers = markers, 
                             dot.color = "darkred")
    
  print(results.list)
  expected.elements <- c("plot","data")
  expect_setequal(names(results.list), expected.elements)
  
})  
