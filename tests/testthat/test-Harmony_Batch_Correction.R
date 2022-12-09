test_that("Harmony returns seurat object with adjusted embeddings", {
  
  so_demo <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  so_harmonized <- harmony_batch_correct(so = so_demo,
                        genes_to_add = c("Epcam"),
                        group.by.vars = c("orig.ident"))
  
  expected_elements = c("pca","umap", "tsne", "harmony")
  expect_setequal(names(so_harmonized@reductions), expected_elements)
})
