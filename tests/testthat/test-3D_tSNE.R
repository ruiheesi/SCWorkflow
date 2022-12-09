test_that("Produce 3D tsne plot and return tsne coordinates", {
  
    obj <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))

    tsne.res <- tSNE3D(object = obj,
                       color.variable = "orig_ident",
                       label.variable = "Likely_CellType",
                       fileName = "tsneplot.html",
                       save.plot = TRUE)
    
    expected.elements = c("plot","data")
    expect_setequal(names(tsne.res), expected.elements)
    
})