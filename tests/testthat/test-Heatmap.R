test_that("Produce heatmap and return plot and filtered dataframe: TEC data", {
    cr.object <- getParamHM("TEC")
    output <- do.call(heatmapSC, cr.object)
            
    expect_snapshot_file(
      ggsave(
         "output/TEC_heatmap.png",
          output$plot,
          width = 10,
          height = 10
        ),
        "TEC_heatmap.png"
    )
            
    expect_type(output, "list")
    expected.elements = c("plot", "data")
    expect_setequal(names(output), expected.elements)
})

test_that("Heatmap scaled vs unscaled", {
  cr.object <- getParamHM("TEC")
  cr.object$scale.data <- FALSE
  output2 <- do.call(heatmapSC, cr.object)
  
  expect_snapshot_file(
    ggsave(
      "output/TEC_heatmap_unscaled.png",
      output2$plot,
      width = 10,
      height = 10
    ),
    "TEC_heatmap_unscaled.png"
  )
  
  #compare scaled (a) vs. nonscaled data (b) to be different
  a <- rowMeans(as.data.frame.matrix(output$data)[, -1])
  b <- rowMeans(as.data.frame.matrix(output2$data)[, -1])
  
  expect_false(isTRUE(all.equal(a, b)))
  
})

test_that("Heatmap run with bad gene name", {
  cr.object <- getParamHM("TEC")
  
  #Spike in with bad gene name
  cr.object$transcripts <- c("badgene", cr.object$transcripts)
  
  expect_warning(do.call(heatmapSC, cr.object), "^There are")
  
  #Actual warning:
  #There are 1 gene(s) absent from dataset:'badgene'.
  #Possible reasons are that gene is not official gene symbol
  #or gene is not highly expressed and has been filtered.
  
})

test_that("Heatmap run with duplicate gene names", {
  cr.object <- getParamHM("TEC")
  
  #Spike in with bad gene name
  cr.object$transcripts <- c("Map7", cr.object$transcripts)
  
  expect_warning(do.call(heatmapSC, cr.object),
                 "The following duplicate genes were removed: Map7")
  
})

test_that("Heatmap run with error for missing all genes", {
  cr.object <- getParamHM("TEC")
  cr.object$transcripts <-
    c("APOE", "ARG1", "CD38", "CD3D", "CD3E", "CD3G")
  expect_error(do.call(heatmapSC, cr.object),
               "No genes listed are found in dataset.")
  
})

test_that("Produce heatmap - Chariou data", {
  cr.object <- getParamHM("Chariou")
  output <- do.call(heatmapSC, cr.object)
  
  expect_snapshot_file(
    ggsave(
      "output/Chariou_heatmap.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "Chariou_heatmap.png"
  )
  expect_type(output, "list")
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
})

test_that("Produce heatmap - PBMC single data", {
  cr.object <- getParamHM("pbmc-single")
  output <- do.call(heatmapSC, cr.object)
  
  expect_snapshot_file(
    ggsave(
      "output/pbmc-single_heatmap.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "pbmc-single_heatmap.png"
  )
  
  expect_type(output, "list")
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
})

test_that("Produce heatmap with protein using NSCLC multi data", {
  cr.object <- getParamHM("nsclc-multi")
  cr.object$legend.font = 3
  output <- do.call(heatmapSC, cr.object)
  
  expect_snapshot_file(
    ggsave(
      "output/nsclc-multi_heatmap.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "nsclc-multi_heatmap.png"
  )
  
  expect_type(output, "list")
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
})

test_that("Produce heatmap with filtered dataframe - BRCA data", {
  cr.object <- getParamHM("BRCA")
  output <- do.call(heatmapSC, cr.object)
  
  expect_snapshot_file(
    ggsave(
      "output/BRCA_heatmap.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "BRCA_heatmap.png"
  )
  
  expect_type(output, "list")
  expected.elements = c("plot", "data")
  expect_setequal(names(output), expected.elements)
})  