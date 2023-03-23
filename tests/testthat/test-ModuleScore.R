# Dataset Testing
test_that("modScore returns metadata with scores and cell calls for tec", {
  
  tec <- getModuleScoreParam("tec")
  
  suppressWarnings(modscore.demo <- do.call(modScore, tec))
  
  expected_elements <- c("Likely_CellType",tec$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$object@meta.data)), 
                          failure_message = "modscore results not found")
  
})

test_that("modScore returns metadata with scores and cell calls for chariou", {
  
  chariou <- getModuleScoreParam("chariou")
  
  suppressWarnings(modscore.demo <- do.call(modScore, chariou))
  
  expected_elements <- c("Likely_CellType",chariou$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("modScore returns metadata with scores and cell calls for 
          nsclc.single", {
  
  nsclc.single <- getModuleScoreParam("nsclc.single")
  
  suppressWarnings(modscore.demo <- do.call(modScore, nsclc.single))
  
  expected_elements <- c("Likely_CellType",nsclc.single$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("modScore returns metadata with scores and cell calls for 
          nsclc.multi", {
  
  nsclc.multi <- getModuleScoreParam("nsclc.multi")
  
  suppressWarnings(modscore.demo <- do.call(modScore, nsclc.multi))
  
  expected_elements <- c("Likely_CellType",nsclc.multi$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$object@meta.data)), 
         failure_message = "modscore results not found")
  
})

test_that("modScore returns metadata with scores and cell calls for brca", {
  
  brca <- getModuleScoreParam("brca")
  
  suppressWarnings(modscore.demo <- do.call(modScore, brca))
  
  expected_elements <- c("Likely_CellType",brca$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$object@meta.data)), 
         failure_message = "modscore results not found")
  
})

# Error Testings
test_that("modScore detects when no genes are found in the data", {
  
  chariou <- getModuleScoreParam("chariou")
  
  expect_error(suppressWarnings(modscore.demo <- modScore(
    object = chariou$object, 
    samples.subset = chariou$samples.subset,
    sample.to.display = chariou$sample.to.display,
    marker.table = apply(chariou$marker.table,2, function(x) toupper(x)),
    celltypes = chariou$celltypes,
    general.class = chariou$general.class,
    lvl.df = chariou$levels_df_demo,
    nbins = 10), "No genes from list was found in data"))
  
})

test_that("modScore detects when threshold number does not match number of 
          cells to analyze", {
  
  chariou <- getModuleScoreParam("chariou")
  
  expect_error(suppressWarnings(modscore.demo <- modScore(
    object = chariou$object, 
    samples.subset = chariou$samples.subset,
    sample.to.display = chariou$sample.to.display,
    marker.table = chariou$marker.table,
    celltypes = chariou$celltypes,
    general.class = chariou$general.class,
    lvl.df = chariou$levels_df_demo,
    threshold = rep(0.1,5),
    nbins = 10), 
    "Threshold length does not match # celltypes to analyze"))
  
})
