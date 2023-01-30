
obj <- select_crobject("TEC")
metadata_column <- "Likely_CellType"
sample_column <- "orig.ident"
genes_column <- "Gene_Names"
cell_type_column <- "Cluster_Names"

test_that("dotplot produced and contingency table returned", {
  
  markers <- c("Apoe","Arg1","Cd38","Cd3d","Cd3e","Cd3g","Cd4","Cd40","Cd86","Cd8a","Csf1r","Cx3cr1","Cxcl10",
               "Cxcl9","F13a1","Fn1","G0s2","Il2ra","Itgam","Lgals3","Ly6c1","Mafb","Nos2","Nr4a1","S100a8","S100a9","Vcan")
  celltypes <- c("CD4_T","cDCs","Monocytes","M1","Neutrophils","B_cells","Tregs","pDC","CD8_T","Macrophages")
  
  results.list <- DotplotMet(object = obj,
                             metadata = metadata_column, 
                             sample.column = sample_column, 
                             cell.type = celltypes,
                             markers = markers, 
                             dot.color = "darkred")
    
  expected.elements <- c("plot","data")
  expect_setequal(names(results.list), expected.elements)
  
})  

test_that("dotplot run with message for duplicate genes", {
  
  markers <- c("Apoe","Apoe","Arg1","Arg1","Cd38","Cd3d","Cd3e","Cd3g","Cd4","Cd40","Cd86","Cd8a","Csf1r","Cx3cr1","Cxcl10",
               "Cxcl9","F13a1","Fn1","G0s2","Il2ra","Itgam","Lgals3","Ly6c1","Mafb","Nos2","Nr4a1","S100a8","S100a9","Vcan")
  celltypes <- c("CD4_T","cDCs","Monocytes","M1","Neutrophils","B_cells","Tregs","pDC","CD8_T","Macrophages")
  
  expect_message(results.list <- DotplotMet(object = obj,
                             metadata = metadata_column, 
                             sample.column = sample_column, 
                             cell.type = celltypes,
                             markers = markers, 
                             dot.color = "darkred"), "The following duplicate genes removed: Apoe")
  

})  

test_that("dotplot run with message for missing genes", {
  
  obj <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  metadata_column <- "Likely_CellType"
  sample_column <- "orig.ident"
  genes_column <- "Gene_Names"
  cell_type_column <- "Cluster_Names"
  markers <- c("Apoe","Arg1","Cd38","Cd3d","Cd3e","Cd3g","Cd4","Cd40","Cd86","Cd8a","Csf1r","Cx3cr1","Cxcl10",
               "Cxcl9","F13a1","Fn1","G0s2","Il2ra","Itgam","Lgals3","Ly6c1","Mafb","Nos2","Nr4a1","S100a8","S100a9","Vcan",
               "Adgre1","Ccr2")
  celltypes <- c("CD4_T","cDCs","Monocytes","M1","Neutrophils","B_cells","Tregs","pDC","CD8_T","Macrophages")
  
  expect_message(results.list <- DotplotMet(object = obj,
                                            metadata = metadata_column, 
                                            sample.column = sample_column, 
                                            cell.type = celltypes,
                                            markers = markers, 
                                            dot.color = "darkred"), 
                 "2 genes are absent from dataset:'Adgre1', 'Ccr2'. Possible reasons are that gene is not official gene symbol or gene is not highly expressed and has been filtered.")
  
}) 

test_that("dotplot run with message for missing all genes, eg. human genes instead of mouse", {
  
  obj <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  metadata_column <- "Likely_CellType"
  sample_column <- "orig.ident"
  genes_column <- "Gene_Names"
  cell_type_column <- "Cluster_Names"
  markers <- c("APOE","ARG1","CD38","CD3D","CD3E","CD3G")
  celltypes <- c("CD4_T","cDCs","Monocytes","M1","Neutrophils","B_cells","Tregs","pDC","CD8_T","Macrophages")
  
  expect_error(results.list <- DotplotMet(object = obj,
                                            metadata = metadata_column, 
                                            sample.column = sample_column, 
                                            cell.type = celltypes,
                                            markers = markers, 
                                            dot.color = "darkred"), 
                                            "No genes listed are found in dataset.")
  
}) 



