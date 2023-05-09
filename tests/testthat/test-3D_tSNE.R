#These lines are necessary to launch Kaleido properly (temporary patch):
reticulate::py_run_string(
  "import sys;print(sys.version); sys.path.append('/rstudio-files/R_environments/single-cell-rna-seq-r4'); print(sys.path)"
)

test_that("Produce 3D tsne plot and return tsne coordinates - TEC Data", {
  cr.object <- getParam3D("TEC")
  output <- do.call(tSNE3D, cr.object)
  
  expect_snapshot_file(
    plotly::save_image(output$plot,
                       file = "output/TEC_plotly.png"),
    "TEC_plotly.png"
  )
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Run 3DTSNE with error for color selection - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$color.variable <- "Likely_CellType"
  expect_error(output <- do.call(tSNE3D, cr.object),
               "^The metadata variable selected for color")
  
})

test_that("Run 3DTSNE with error for color selection - TEC Data", {
  cr.object <- getParam3D("TEC")
  cr.object$label.variable <- "Likely_CellType"
  expect_error(output <- do.call(tSNE3D, cr.object),
               "^The metadata variable selected for labeling")
  
})

test_that("Produce 3D tsne plot and return tsne coordinates - Chariou Data", {
  cr.object <- getParam3D("Chariou")
  output <- do.call(tSNE3D, cr.object)
  
  expect_snapshot_file(
    plotly::save_image(output$plot,
                       file = "output/Chariou_plotly.png"),
    "Chariou_plotly.png"
  )
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Produce 3D tsne plot and return tsne - PBMC-single Data", {
  cr.object <- getParam3D("pbmc-single")
  output <- do.call(tSNE3D, cr.object)
  
  expect_snapshot_file(
    plotly::save_image(output$plot,
                       file = "output/PBMC_single_plotly.png"),
    "PBMC_single_plotly.png"
  )
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Produce 3D tsne plot and return tsne - NSCLC-multi Data", {
  cr.object <- getParam3D("nsclc-multi")
  output <- do.call(tSNE3D, cr.object)
  
  expect_snapshot_file(
    plotly::save_image(output$plot,
                       file = "output/NSCLC_multi_plotly.png"),
    "NSCLC_multi_plotly.png"
  )
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Produce 3D tsne plot and return tsne - BRCA Data", {
  cr.object <- getParam3D("BRCA")
  output <- do.call(tSNE3D, cr.object)
  
  expect_snapshot_file(
    plotly::save_image(output$plot,
                       file = "output/BRCA_plotly.png"),
    "BRCA_plotly.png"
  )
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
  
})
