

test_that("Test Filter Seurat Object by Metadata using Downsampled TEC (Mouse) data", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  FilterByMeta.result <- filterSeuratObjectByMetadata(object = obj,
                                                      samples.to.include = 'c("1_Embryo_13_5","3_Newborn")',
                                                      sample.name = "orig.ident",
                                                      category.to.filter = "seurat_clusters",
                                                      values.to.filter = "3")
  
  
#  ggsave(file="BeforeFiltering.pdf", FilterByMeta.result$plot1)
#  ggsave(file="AfterFiltering.pdf", FilterByMeta.result$plot2)
  expect_type(FilterByMeta.result,"list")
  expected.elements = c("object", "plot1", "plot2")
  expect_setequal(names(FilterByMeta.result), expected.elements)
  
})




test_that("Test Filter Seurat Object by Metadata using Chariou (Mouse) dataset", {
  obj <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
  FilterByMeta.result <- filterSeuratObjectByMetadata(object = obj,
                                                      samples.to.include = 'c("CD8dep","NHSIL12")',
                                                      sample.name = "orig.ident",
                                                      category.to.filter = "seurat_clusters",
                                                      values.to.filter = "2")
  
  
#  ggsave(file="BeforeFiltering.pdf", FilterByMeta.result$plot1)
#  ggsave(file="AfterFiltering.pdf", FilterByMeta.result$plot2)
  expect_type(FilterByMeta.result,"list")
  expected.elements = c("object", "plot1", "plot2")
  expect_setequal(names(FilterByMeta.result), expected.elements)
  
})


test_that("Test Filter Seurat Object by Metadata using BRCA (Human) dataset", {
  obj <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
  FilterByMeta.result <- filterSeuratObjectByMetadata(object = obj,
                                                      samples.to.include = 'c("CID3586","CID3946", "CID4513" )',
                                                      sample.name = "orig.ident",
                                                      category.to.filter = "seurat_clusters",
                                                      values.to.filter = "1")
  
  
#  ggsave(file="BeforeFiltering.pdf", FilterByMeta.result$plot1)
#  ggsave(file="AfterFiltering.pdf", FilterByMeta.result$plot2)
  expect_type(FilterByMeta.result,"list")
  expected.elements = c("object", "plot1", "plot2")
  expect_setequal(names(FilterByMeta.result), expected.elements)
  
})


test_that("Test Filter Seurat Object by Metadata using NSCLCmulti (Human) dataset", {
  obj <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  FilterByMeta.result <- filterSeuratObjectByMetadata(object = obj,
                                                      samples.to.include = 'c("Donor_1","Donor_4")',
                                                      sample.name = "orig.ident",
                                                      category.to.filter = "seurat_clusters",
                                                      values.to.filter = "16")
  
  
#  ggsave(file="BeforeFiltering.pdf", FilterByMeta.result$plot1)
#  ggsave(file="AfterFiltering.pdf", FilterByMeta.result$plot2)
  expect_type(FilterByMeta.result,"list")
  expected.elements = c("object", "plot1", "plot2")
  expect_setequal(names(FilterByMeta.result), expected.elements)
  
})




test_that("Test Filter Seurat Object by Metadata using NSCLCsingle (Human) dataset", {
  obj <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
  FilterByMeta.result <- filterSeuratObjectByMetadata(object = obj,
                                                      samples.to.include = 'c("NSCLC_Single")',
                                                      sample.name = "orig.ident",
                                                      category.to.filter = "seurat_clusters",
                                                      values.to.filter = "5")
  
  
#  ggsave(file="BeforeFiltering.pdf", FilterByMeta.result$plot1)
#  ggsave(file="AfterFiltering.pdf", FilterByMeta.result$plot2)
  expect_type(FilterByMeta.result,"list")
  expected.elements = c("object", "plot1", "plot2")
  expect_setequal(names(FilterByMeta.result), expected.elements)
  
})

