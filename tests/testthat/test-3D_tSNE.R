test_that("Produce 3D tsne plot and return tsne coordinates - TEC Data", {
  
    CRObject <- getparam3d("TEC")
    output <- do.call(tSNE3D,CRObject)
    
    expect_snapshot_file("output","TEC_tsneplot.html")
    expected.elements = c("plot","data")
    expect_setequal(names(output), expected.elements)
    
})

test_that("Run 3DTSNE with error for color selection - TEC Data", {
  
  CRObject <- getparam3d("TEC")
  CRObject$color.variable <- "Likely_CellType"
  expect_error(output <- do.call(tSNE3D,CRObject),
            "^The metadata variable selected for color")
  
})

test_that("Run 3DTSNE with error for color selection - TEC Data", {
  
  CRObject <- getparam3d("TEC")
  CRObject$label.variable <- "Likely_CellType"
  expect_error(output <- do.call(tSNE3D,CRObject),
            "^The metadata variable selected for labeling")
  
})

test_that("Produce 3D tsne plot and return tsne coordinates - Chariou Data", {
  
  CRObject <- getparam3d("Chariou")
  output <- do.call(tSNE3D,CRObject)
  
  expect_snapshot_file("output","Chariou_tsneplot.html")
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Produce 3D tsne plot and return tsne coordinates - NSCLC-single Data", {
  
  CRObject <- getparam3d("nsclc-single")
  output <- do.call(tSNE3D,CRObject)
  expect_snapshot_file("output","nsclc-single_tsneplot.html")
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Produce 3D tsne plot and return tsne coordinates - NSCLC-multi Data", {
  
  CRObject <- getparam3d("nsclc-multi")
  output <- do.call(tSNE3D,CRObject)
  expect_snapshot_file("output","nsclc-multi_tsneplot.html")
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
  
})

test_that("Produce 3D tsne plot and return tsne coordinates - BRCA Data", {
  
  CRObject <- getparam3d("BRCA")
  output <- do.call(tSNE3D,CRObject)
  expect_snapshot_file("output","BRCA_tsneplot.html")
  expected.elements = c("plot","data")
  expect_setequal(names(output), expected.elements)
  
})


