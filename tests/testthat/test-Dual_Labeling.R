CRObject <- getparamdl("TEC")

test_that("Test Dual labeling TEC Data", {
  
  dual.label.result <- do.call(DualLabeling,CRObject)  
  #grid.arrange(dual.label.result$plot)
  
  tpath <- paste0(test_path(),"/output")
  ggsave(file=paste0(tpath,"/tsne_duallabel.png"),dual.label.result$plot)
         
  expect_snapshot_file(tpath,"tsne_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(dual.label.result), expected.elements)
  expect_length(dual.label.result$plot$grobs,7)
  
})

test_that("Dual labeling with missing gene", {
  
  CRObject_Test <- CRObject
  CRObject_Test$marker1 <- "Cd8"
  expect_error(do.call(DualLabeling,CRObject_Test),
               "is not found in dataset$")
})

test_that("Dual labeling with umap", {
  
  CRObject_Test <- CRObject
  CRObject_Test$data.reduction = "umap"
  dual.label.result <- do.call(DualLabeling,CRObject_Test)
  #grid.arrange(dual.label.result$plot)
  
  tpath <- paste0(test_path(),"/output")
  ggsave(file=paste0(tpath,"/umap_duallabel.png"),dual.label.result$plot)
  
  expect_snapshot_file(tpath,"umap_duallabel.png")
  expected.elements = c("so","plot")
  expect_setequal(names(dual.label.result), expected.elements)
})

test_that("Dual labeling without density heatmap", {
  
  CRObject_Test <- CRObject
  CRObject_Test$density.heatmap = FALSE
  dual.label.result <- do.call(DualLabeling,CRObject_Test)
  
  tpath <- paste0(test_path(),"/output")
  ggsave(file=paste0(tpath,"/nodensity_duallabel.png"),
         dual.label.result$plot)
  
  expect_snapshot_file(tpath,"nodensity_duallabel.png")
  expect_length(dual.label.result$plot$grobs,6)

})

test_that("Test Dual labeling Chariou Data", {
  
  CRObject <- getparamdl("Chariou")
  dual.label.result <- do.call(DualLabeling,CRObject)  

  tpath <- paste0(test_path(),"/output")
  ggsave(file=paste0(tpath,"/chariou_duallabel.png"),
         dual.label.result$plot)
  
  expect_snapshot_file(tpath,"chariou_duallabel.png")
  
  expected.elements = c("so","plot")
  expect_setequal(names(dual.label.result), expected.elements)
  expect_length(dual.label.result$plot$grobs,7)
  
})


        