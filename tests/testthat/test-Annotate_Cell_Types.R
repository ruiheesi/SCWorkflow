test_that("Test Annotate Cell Types for UMAP plots (Mouse TEC dataset)", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  Anno.result <- annotateCellTypes(object = obj)
  ggsave("output/TEC_annotateCellTypes.p1.png",Anno.result$p1, width = 10, height = 10)
  expect_snapshot_file("output","TEC_annotateCellTypes.p1.png")
  ggsave("output/TEC_annotateCellTypes.p2.png",Anno.result$p2, width = 10, height = 10)
  expect_snapshot_file("output","TEC_annotateCellTypes.p2.png")
  expect_type(Anno.result,"list")
  expected.elements = c("object", "p1", "p2")
  expect_setequal(names(Anno.result), expected.elements) 
})


test_that("Test Annotate Cell Types for TSNE plots (Mouse TEC dataset)", {
    obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
    Anno.result <- annotateCellTypes(object = obj, reduction.type = "tsne")
    ggsave("output/TEC_annotateCellTypes_tsne.p1.png",Anno.result$p1, width = 10, height = 10)
    expect_snapshot_file("output","TEC_annotateCellTypes_tsne.p1.png")
    ggsave("output/TEC_annotateCellTypes_tsne.p2.png",Anno.result$p2, width = 10, height = 10)
    expect_snapshot_file("output","TEC_annotateCellTypes_tsne.p2.png")
    expect_type(Anno.result,"list")
    expected.elements = c("object", "p1", "p2")
    expect_setequal(names(Anno.result), expected.elements) 
  })


test_that("Test Annotate Cell Types with FineTuning (Mouse TEC dataset)", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  Anno.result <- annotateCellTypes(object = obj, do.finetuning = TRUE)
  ggsave("output/TEC_annotateCellTypes_finetune.p1.png",Anno.result$p1, width = 10, height = 10)
  expect_snapshot_file("output","TEC_annotateCellTypes_finetune.p1.png")
  ggsave("output/TEC_annotateCellTypes_finetune.p2.png",Anno.result$p2, width = 10, height = 10)
  expect_snapshot_file("output","TEC_annotateCellTypes_finetune.p2.png")
  expect_type(Anno.result,"list")
  expected.elements = c("object", "p1", "p2")
  expect_setequal(names(Anno.result), expected.elements) 
})
    

# test_that("Test Annotate Cell Types using legend.dot.size (Mouse TEC dataset)", {
#   obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
#   Anno.result <- annotateCellTypes(object = obj, legend.dot.size = 1)
#   expected.elements = c("object", "p1", "p2")
#   expect_setequal(names(Anno.result), expected.elements) 
# })



test_that("Test Annotate Cell Types using Chariou (Mouse) dataset", {    
  obj <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
  Anno.result <- annotateCellTypes(object = obj, species = "Mouse")
  ggsave("output/Chariou_annotateCellTypes.p1.png",Anno.result$p1, width = 10, height = 10)
  expect_snapshot_file("output","Chariou_annotateCellTypes.p1.png")
  ggsave("output/Chariou_annotateCellTypes.p2.png",Anno.result$p2, width = 10, height = 10)
  expect_snapshot_file("output","Chariou_annotateCellTypes.p2.png")
  expect_type(Anno.result,"list")
  expected.elements = c("object", "p1", "p2")
  expect_setequal(names(Anno.result), expected.elements) 
})




test_that("Test Annotate Cell Types using BRCA (Human) dataset", {    
    obj <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
    Anno.result <- annotateCellTypes(object = obj, species = "Human")
    ggsave("output/BRCA_annotateCellTypes.p1.png",Anno.result$p1, width = 10, height = 10)
    expect_snapshot_file("output","BRCA_annotateCellTypes.p1.png")
    ggsave("output/BRCA_annotateCellTypes.p2.png",Anno.result$p2, width = 10, height = 10)
    expect_snapshot_file("output","BRCA_annotateCellTypes.p2.png")
    expect_type(Anno.result,"list")
    expected.elements = c("object", "p1", "p2")
    expect_setequal(names(Anno.result), expected.elements) 
  })


test_that("Test Annotate Cell Types using NSCLCmulti (Human) dataset", {    
  obj <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  Anno.result <- annotateCellTypes(object = obj, species = "Human")
  ggsave("output/NSCLCmulti_annotateCellTypes.p1.png",Anno.result$p1, width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_annotateCellTypes.p1.png")
  ggsave("output/NSCLCmulti_annotateCellTypes.p2.png",Anno.result$p2, width = 10, height = 10)
  expect_snapshot_file("output","NSCLCmulti_annotateCellTypes.p2.png")
  expect_type(Anno.result,"list")
  expected.elements = c("object", "p1", "p2")
  expect_setequal(names(Anno.result), expected.elements) 
})


test_that("Test Annotate Cell Types using NSCLCsingle (Human) dataset", {    
  obj <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
  Anno.result <- annotateCellTypes(object = obj, species = "Human")
  ggsave("output/NSCLCsingle_annotateCellTypes.p1.png",Anno.result$p1, width = 10, height = 10)
  expect_snapshot_file("output","NSCLCsingle_annotateCellTypes.p1.png")
  ggsave("output/NSCLCsingle_annotateCellTypes.p2.png",Anno.result$p2, width = 10, height = 10)
  expect_snapshot_file("output","NSCLCsingle_annotateCellTypes.p2.png")
  expect_type(Anno.result,"list")
  expected.elements = c("object", "p1", "p2")
  expect_setequal(names(Anno.result), expected.elements) 
})