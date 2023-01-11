
test_that("dotplot produced and contingency table returned", {
  library(Seurat)
  obj <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  metadata_column <- "Likely_CellType"
  sample_column <- "orig.ident"
  genes_column <- "Gene_Names"
  cell_type_column <- "Cluster_Names"
  
  markers <- c("Adgre1","Apoe","Arg1","Ccr2","Cd163","Cd38","Cd3d","Cd3e","Cd3g","Cd4","Cd40","Cd68","Cd86","Cd8a","Csf1r","Csf3r","Cx3cr1","Cxcl10","Cxcl9","F13a1","Fn1","Foxp3","G0s2","Gzmk","Il1b","Il2ra","Itgam","Lgals3","Ly6c1","Ly6g","Mafb","Mmp12","Mmp13","Mrc1","Nos2","Nr4a1","S100a8","S100a9","Sell_neg","Siglec1","Vcan")
  celltypes <- c("CD4_T","cDCs","Monocytes","M1","Neutrophils","B_cells","Tregs","pDC","CD8_T","Macrophages")
  
  results.list <- DotplotMet(object = obj,
                             metadata = metadata_column, 
                             sample.column = sample_column, 
                             cell.type = celltypes,
                             markers = markers, 
                             dot.color = "blue")
    
  print(results.list)
  expected.elements <- c("plot","data")
  expect_setequal(names(results.list), expected.elements)
  
})  
