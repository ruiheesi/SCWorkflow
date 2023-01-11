test_that("correct number of Sample Names returned", {
  
  # load data
  seurat.object <-
    readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  # run function
  samples <- SampleNames(SO = seurat.object)
  
  # test
  expect_true(ncol(samples) == length(unique(seurat.object$orig.ident)))
  
})

test_that("all Sample Names returned", {
  
  # load data
  seurat.object <-
    readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  # run function
  samples <- SampleNames(SO = seurat.object)
  
  # test
  expect_equal(sort(colnames(samples)), sort(unique(seurat.object$orig.ident)))
  
})