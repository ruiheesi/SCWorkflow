# load data
seurat.object <-readRDS(test_path("fixtures", "SO_moduleScore.rds"))

# run tests
test_that("reductions are returned", {
  
  metadata <-
    MetadataTable(SO = seurat.object,
                  return.cell.embeddings = TRUE)
  
  expect_true(sum(grepl("PC_1$|PC_2$", colnames(metadata))) == 2)
  expect_true(sum(grepl("UMAP_1|UMAP_2", colnames(metadata))) == 2)
  expect_true(sum(grepl("tSNE_1|tSNE_2", colnames(metadata))) == 2)
  
})

test_that("function returns correct class", {
  
  metadata <-
    MetadataTable(SO = seurat.object,
                  return.cell.embeddings = TRUE)
  expect_s3_class(metadata, "data.frame")

})