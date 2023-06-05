test_that("Violin plot works for TEC data", {
  tec.data = selectViolin("TEC")
  
  violin_test = do.call(violinPlot, tec.data)
  
  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for Chariou data", {
  chariou.data = selectViolin("Chariou")
  
  violin_test = do.call(violinPlot, chariou.data)
  
  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for pbmc.single data", {
  pbmc.single = selectViolin("pbmc.single")
  
  violin_test = do.call(violinPlot, pbmc.single)
  
  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for nsclc.multi data", {
  nsclc.multi = selectViolin("nsclc.multi")
  
  violin_test = do.call(violinPlot, nsclc.multi)
  
  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

test_that("Violin plot works for brca data", {
  brca = selectViolin("brca")
  
  violin_test = do.call(violinPlot, brca)
  
  expected_elements = c("gg", "ggplot")
  expect_setequal(class(violin_test), expected_elements)
  
})

## Check code detects warnings and errors ##

test_that("Violin plot stops when no query genes are found in the data", {
  pbmc.single <- selectViolin("pbmc.single")
  
  expect_error(
    violinPlot(
      object = pbmc.single$object,
      group.by = pbmc.single$group.by,
      group.subset = pbmc.single$group.subset,
      genes.of.interest = paste("jibberish", 1:5, sep =
                                  "_")
    ),
    "No query genes were found in the dataset."
  )
  
})

test_that("Violin plot stops when ident of interest is not found in seurat", {
  pbmc.single <- selectViolin("pbmc.single")
  
  expect_error(
    violinPlot(
      object = pbmc.single$object,
      group.by = "jibberish",
      group.subset = pbmc.single$group.subset,
      genes.of.interest = pbmc.single$genes.of.interest
    ),
    "Unable to find ident of interest in metadata."
  )
  
})

test_that("Violin plot stops when group of interest is empty", {
  pbmc.single <- selectViolin("pbmc.single")
  
  expect_error(
    violinPlot(
      object = pbmc.single$object,
      group.by = pbmc.single$group.by,
      group.subset = paste("jibberish", 1:5, sep = "_"),
      genes.of.interest = pbmc.single$genes.of.interest
    ),
    "No groups were found in the selected ident."
  )
  
})

test_that("Violin plot stops when user attempts to rename group.by as
          Gene, Expression, or Scaled",
          {
            pbmc.single <- selectViolin("pbmc.single")
            
            expect_error(
              violinPlot(
                object = pbmc.single$object,
                group.by = pbmc.single$group.by,
                group.subset = pbmc.single$group.subset,
                genes.of.interest = pbmc.single$genes.of.interest,
                rename.ident = "Gene"
              ),
              "New ident name cannot be one of Gene, Expression, or scaled."
            )
            
          })
