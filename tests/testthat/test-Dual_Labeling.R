test_that("Test Dual labeling TEC Data", {
  CRObject <- getParamDL("TEC")
  output <- do.call(dualLabeling,CRObject)  
  
  expect_snapshot_file(
    ggsave("TEC_duallabel.png",
         output$plot, 
         width = 10, 
         height = 10))
  expect_snapshot_file("TEC_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Dual labeling with missing gene", {
  
  CRObject <- getParamDL("TEC")
  CRObject$marker1 <- "Cd8"
  expect_error(do.call(dualLabeling,CRObject),
               "is not found in dataset$")
})

test_that("Dual labeling with umap", {
  
  CRObject <- getParamDL("TEC")
  CRObject$data.reduction = "umap"
  output <- do.call(dualLabeling,CRObject)
  
  ggsave(file="umap_duallabel.png",
         output$plot, 
         width = 10, 
         height = 10)
  expect_snapshot_file("umap_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
})

test_that("Dual labeling without density heatmap", {
  
  CRObject <- getParamDL("TEC")
  CRObject$density.heatmap = FALSE
  output <- do.call(dualLabeling,CRObject)
  
  expect_snapshot_file(ggsave(file="nodensity_duallabel.png",
         output$plot, 
         width = 10, 
         height = 10))
  expect_snapshot_file("nodensity_duallabel.png")
  expect_length(output$plot$grobs,6)

})

test_that("Test Dual labeling Chariou Data", {
  
  CRObject <- getParamDL("Chariou")
  output <- do.call(dualLabeling,CRObject)  

  expect_snapshot_file(ggsave("chariou_duallabel.png",
         output$plot, 
         width = 10, 
         height = 10))
  expect_snapshot_file("chariou_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Test Dual labeling PBMC-single Data", {
  
  CRObject <- getParamDL("pbmc-single")
  output <- do.call(dualLabeling,CRObject)  
  
  expect_snapshot_file(ggsave("pbmc-single_duallabel.png",
         output$plot, 
         width = 10, 
         height = 10))
  expect_snapshot_file("pbmc-single_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Test Dual labeling NSCLC-multi Data", {
  
  CRObject <- getParamDL("nsclc-multi")
  output <- do.call(dualLabeling,CRObject)  
  
  expect_snapshot_file(ggsave("nsclc-multi_duallabel.png",
         output$plot, 
         width = 10, 
         height = 10))
  expect_snapshot_file("nsclc-multi_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Test Dual labeling BRCA Data", {
  
  CRObject <- getParamDL("BRCA")
  output <- do.call(dualLabeling,CRObject)  
  
  expect_snapshot_file(ggsave("BRCA_duallabel.png",
         output$plot, 
         width = 10, 
         height = 10))
  expect_snapshot_file("BRCA_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})
