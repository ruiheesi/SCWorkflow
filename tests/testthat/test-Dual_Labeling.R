test_that("Test Dual labeling TEC Data", {
  CRObject <- getParamDL("TEC")
  output <- do.call(dualLabeling, CRObject)
  
  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs, 7)
  
  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plot),
        "TEC_duallabel.png"
    )
})

test_that("Dual labeling with missing gene", {
  CRObject <- getParamDL("TEC")
  CRObject$marker.1 <- "Cd8"
  expect_error(do.call(dualLabeling, CRObject),
               "is not found in dataset$")
})

test_that("Dual labeling with umap", {
  CRObject <- getParamDL("TEC")
  CRObject$data.reduction = "umap"
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs, 7)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plot),
             "umap_duallabel.png"
    )
})

test_that("Dual labeling without density heatmap", {
  CRObject <- getParamDL("TEC")
  CRObject$density.heatmap = FALSE
  output <- do.call(dualLabeling, CRObject)
  expect_length(output$plot$grobs, 6)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plot),
             "nodensity_duallabel.png"
    )
})

test_that("Test Dual labeling Chariou Data", {
  CRObject <- getParamDL("Chariou")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs, 7)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plot),
       "chariou_duallabel.png"
    )
})

test_that("Test Dual labeling PBMC-single Data", {
  CRObject <- getParamDL("pbmc-single")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs, 7)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plot),
             "pbmc-single_duallabel.png"
    )
})

test_that("Test Dual labeling NSCLC-multi Data", {
  CRObject <- getParamDL("nsclc-multi")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs, 7)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plot),
             "nsclc-multi_duallabel.png"
    )
})

test_that("Test Dual labeling BRCA Data", {
  CRObject <- getParamDL("BRCA")
  output <- do.call(dualLabeling, CRObject)

  expected.elements = c("object", "plot")
  expect_setequal(names(output), expected.elements)
  expect_length(output$plot$grobs, 7)

  skip_on_ci()
  expect_snapshot_file(
    .drawdualplot(output$plot),
             "BRCA_duallabel.png"
    )
})
