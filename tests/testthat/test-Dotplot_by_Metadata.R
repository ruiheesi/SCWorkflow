test_that("dotplot run with normal parameters - TEC Data", {
  tec.data <- getParamDP("TEC")
  output <- do.call(dotPlotMet, tec.data)
  expect_snapshot_file(ggsave(
    "output/TEC_dotplot.png",
    output$plot,
    width = 10,
    height = 10
  ), "TEC_dotplot.png")
  expect_type(output, "list")
  expected.elements <- c("plot", "pct", "exp")
  expect_setequal(names(output), expected.elements)
})

test_that("dotplot run with reversed axes - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data$plot.reverse <- TRUE
  output <- do.call(dotPlotMet, tec.data)
  expect_snapshot_file(ggsave(
    "output/TEC_dotplot_reversed.png",
    output$plot,
    width = 10,
    height = 10
  ), "TEC_dotplot_reversed.png")
  expect_type(output, "list")
  expected.elements <- c("plot", "pct", "exp")
  expect_setequal(names(output), expected.elements)
})

test_that("dotplot run with different color - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data$dot.color <- "red"
  output <- do.call(dotPlotMet, tec.data)
  expect_snapshot_file(ggsave(
    "output/TEC_dotplot_red.png",
    output$plot,
    width = 10,
    height = 10
  ), "TEC_dotplot_red.png")
  expect_type(output, "list")
  expected.elements <- c("plot", "pct", "exp")
  expect_setequal(names(output), expected.elements)
})

test_that("dotplot run with warning for fewer categories - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data$cells <- c(0, 1, 2, 3, 4, 5, 6)
  expect_warning(do.call(dotPlotMet, tec.data),
                 "^There are")
})

test_that("dotplot run with error for non-matched category - TEC Data", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$metadata <- "orig.ident"
  expect_error(
    do.call(dotPlotMet, tec.data),
    "^At least 2 metadata categories you wish to plot should"
  )
})

test_that("dotplot run with warning for duplicate genes", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$markers <- c("Map7", tec.data$markers)
  expect_warning(
    do.call(dotPlotMet, tec.data),
    fixed = TRUE,
    "The following duplicate genes were removed: Map7"
  )
  
})

test_that("dotplot run with warning for missing genes", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$markers <- c("Adgre1", "Ccr2", tec.data$markers)
  expect_warning(
    do.call(dotPlotMet, tec.data), 
    "^There are "
  )
  
  #Full message:
  #There are 2 gene(s) absent from dataset:'Adgre1', 'Ccr2'.
  #Possible reasons are that gene is not official gene symbol or
  #gene is not highly expressed and has been filtered.
})

test_that("dotplot run with error for missing all genes", {
  tec.data <- getParamDP("TEC")
  tec.data <- tec.data
  tec.data$markers <- c("APOE", "ARG1", "CD38", "CD3D", "CD3E", "CD3G")
  expect_error(do.call(dotPlotMet, tec.data),
               "No genes listed are found in dataset.")
  
})

test_that("dotplot produced and contingency table returned - Chariou Data", {
  chariou.data <- getParamDP("Chariou")
  output <- do.call(dotPlotMet, chariou.data)
  
  expect_snapshot_file(
    ggsave(
      "output/Chariou_dotplot.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "Chariou_dotplot.png"
  )
  
  expect_type(output, "list")
  expected.elements <- c("plot", "pct", "exp")
  expect_setequal(names(output), expected.elements)
  
})

test_that("dotplot produced and contingency table returned - PBMC-single Data",
          {
            pbmc.single.data <- getParamDP("pbmc-single")
            output <- do.call(dotPlotMet, pbmc.single.data)
            
            expect_snapshot_file(
              ggsave(
                "output/pbmc-single_dotplot.png",
                output$plot,
                width = 10,
                height = 10
              ),
              "pbmc-single_dotplot.png"
            )
            
            expect_type(output, "list")
            expected.elements <- c("plot", "pct", "exp")
            expect_setequal(names(output), expected.elements)
            
          })

test_that("dotplot produced and contingency table returned - NSCLC-multi Data",
          {
            nsclc.multi.data <- getParamDP("nsclc-multi")
            expect_warning(output <- do.call(dotPlotMet, nsclc.multi.data),
                           "^There are")
            
            expect_snapshot_file(
              ggsave(
                "output/nsclc-multi_dotplot.png",
                output$plot,
                width = 10,
                height = 10
              ),
              "nsclc-multi_dotplot.png"
            )
            
            expect_type(output, "list")
            expected.elements <- c("plot", "pct", "exp")
            expect_setequal(names(output), expected.elements)
            
          })

test_that("dotplot produced and contingency table returned - BRCA Data", {
  brca.data <- getParamDP("BRCA")
  expect_warning(output <- do.call(dotPlotMet, brca.data),
                 "^There are")
  
  expect_snapshot_file(ggsave(
    "output/brca_dotplot.png",
    output$plot,
    width = 10,
    height = 10
  ),
  "brca_dotplot.png")
  
  expect_type(output, "list")
  expected.elements <- c("plot", "pct", "exp")
  expect_setequal(names(output), expected.elements)
  
})
