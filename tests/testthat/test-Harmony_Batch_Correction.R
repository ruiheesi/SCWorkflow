test_that("Harmony returns seurat object with adjusted embeddings for TEC data", {
  
  TEC <- getparam_Harmony("TEC")
  
  so_harmonized <- harmonyBatchcorrect(so = TEC$object,
                                         variable.features = 200,
                                         genes.to.add = TEC$genes.to.add,
                                         group.by.vars = TEC$group.by.vars)
  
  expected_elements = c("pca","umap", "tsne", "harmony")
  expect_setequal(names(so_harmonized@reductions), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for Chariou data", {
  
  Chariou <- getparam_Harmony("Chariou")
  
  so_harmonized <- harmonyBatchcorrect(so = Chariou$object,
                                         variable.features = 200,
                                         genes.to.add = Chariou$genes.to.add,
                                         group.by.vars = Chariou$group.by.vars)
  
  expected_elements = c("pca","umap", "tsne", "harmony")
  expect_setequal(names(so_harmonized@reductions), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for NSCLC_Single data", {
  
  NSCLC_Single <- getparam_Harmony("NSCLC_Single")
  
  so_harmonized <- harmonyBatchcorrect(so = NSCLC_Single$object,
                                         variable.features = 200,
                                         genes.to.add = NSCLC_Single$genes.to.add,
                                         group.by.vars = NSCLC_Single$group.by.vars)
  
  expected_elements = c("pca","umap", "tsne", "harmony")
  expect_setequal(names(so_harmonized@reductions), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for NSCLC_Multi data", {
  
  NSCLC_Multi <- getparam_Harmony("NSCLC_Multi")
  
  so_harmonized <- harmonyBatchcorrect(so = NSCLC_Multi$object,
                                         variable.features = 200,
                                         genes.to.add = NSCLC_Multi$genes.to.add,
                                         group.by.vars = NSCLC_Multi$group.by.vars)
  
  expected_elements = c("pca","umap","tsne","protein_umap","protein_tsne","harmony")
  expect_setequal(names(so_harmonized@reductions), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for BRCA data", {
  
  BRCA <- getparam_Harmony("BRCA")
  
  so_harmonized <- harmonyBatchcorrect(so = BRCA$object,
                                         variable.features = 100,
                                         genes.to.add = BRCA$genes.to.add,
                                         group.by.vars = BRCA$group.by.vars)
  
  expected_elements = c("pca","umap", "tsne", "harmony")
  expect_setequal(names(so_harmonized@reductions), expected_elements)
})

## Error Testing ##

test_that("Harmony provides warning when genes are not found in the data", {
  
  NSCLC_Single <- getparam_Harmony("NSCLC_Single")
  
  expect_warning(harmonyBatchcorrect(so = NSCLC_Single$object,
                                         variable.features = 200,
                                         genes.to.add = c("wrong_gene","wrong_gene2"),
                                         group.by.vars = NSCLC_Single$group.by.vars), "specified genes were not found and therefore cannot be added")
  
  
})

test_that("Harmony stops when variable features to subset by exceeds number of genes in the data", {
  
  NSCLC_Single <- getparam_Harmony("NSCLC_Single")
  
  expect_error(harmonyBatchcorrect(so = NSCLC_Single$object,
                                       variable.features = 4000,
                                       genes.to.add = NULL,
                                       group.by.vars = NSCLC_Single$group.by.vars), "Number of variable features to subset by cannot exceed the total number of variable genes in the data")
  
  
})