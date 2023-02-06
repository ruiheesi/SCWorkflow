test_that("Test Dual labeling TEC Data", {
  CRObject <- getparamdl("TEC")
  output <- do.call(dualLabeling,CRObject)  
  
  ggsave("output/TEC_duallabel.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","TEC_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Dual labeling with missing gene", {
  
  CRObject <- getparamdl("TEC")
  CRObject$marker1 <- "Cd8"
  expect_error(do.call(dualLabeling,CRObject),
               "is not found in dataset$")
})

test_that("Dual labeling with umap", {
  
  CRObject <- getparamdl("TEC")
  CRObject$data.reduction = "umap"
  output <- do.call(dualLabeling,CRObject)
  
  ggsave(file="output/umap_duallabel.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","umap_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
})

test_that("Dual labeling without density heatmap", {
  
  CRObject <- getparamdl("TEC")
  CRObject$density.heatmap = FALSE
  output <- do.call(dualLabeling,CRObject)
  
  ggsave(file="output/nodensity_duallabel.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nodensity_duallabel.png")
  expect_length(output$plot$grobs,6)

})

test_that("Test Dual labeling Chariou Data", {
  
  CRObject <- getparamdl("Chariou")
  output <- do.call(dualLabeling,CRObject)  

  ggsave("output/Chariou_duallabel.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","chariou_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Test Dual labeling NSCLC-single Data", {
  
  CRObject <- getparamdl("nsclc-single")
  output <- do.call(dualLabeling,CRObject)  
  
  ggsave("output/nsclc-single_duallabel.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nsclc-single_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Test Dual labeling NSCLC-multi Data", {
  
  CRObject <- getparamdl("nsclc-multi")
  output <- do.call(dualLabeling,CRObject)  
  
  ggsave("output/nsclc-multi_duallabel.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","nsclc-multi_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})

test_that("Test Dual labeling BRCA Data", {
  
  CRObject <- getparamdl("BRCA")
  output <- do.call(dualLabeling,CRObject)  
  
  ggsave("output/BRCA_duallabel.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","BRCA_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs,7)
  
})
