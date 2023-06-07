# Dataset Testing
test_that("modScore returns metadata with scores and cell calls for tec", {
  
  tec <- getModuleScoreParam("tec")
  
  suppressWarnings(modscore.demo <- do.call(modScore, tec))
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "tec_MS_1.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "tec_MS_2.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "tec_MS_3.png"
  )
  
  expect_equal(mean(modscore.demo$ms.object$CD8_T),
               0.044506, tolerance = 1e-1)
  
  expect_equal(mean(modscore.demo$ms.object$CD4_T),
               0.044506, tolerance = 1e-1)
  
  expect_equal(mean(modscore.demo$ms.object$Tregs),
               0.0913447, tolerance = 1e-1)
  
  expected_elements <- c("Likely_CellType",tec$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)), 
                          failure_message = "modscore results not found")
  
})

test_that("modScore returns metadata with scores and cell calls for chariou", {

  chariou <- getModuleScoreParam("chariou")

  suppressWarnings(modscore.demo <- do.call(modScore, chariou))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "chariou_MS_1.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "chariou_MS_2.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "chariou_MS_3.png"
  )

  expect_equal(mean(modscore.demo$ms.object$CD8_T),
               0.05458893, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$CD4_T),
               0.05458893, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$Tregs),
               0.07577882, tolerance = 1e-1)

  expected_elements <- c("Likely_CellType",chariou$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
         failure_message = "modscore results not found")

})

test_that("modScore returns metadata with scores and cell calls for
          pbmc.single", {

  pbmc.single <- getModuleScoreParam("pbmc.single")

  suppressWarnings(modscore.demo <- do.call(modScore, pbmc.single))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "pbmc_single_MS_1.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "pbmc_single_MS_2.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "pbmc_single_MS_3.png"
  )

  expect_equal(mean(modscore.demo$ms.object$rand_type1),
               0.2129345, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$rand_type2),
               0.1917001, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$rand_type3),
               0.1798533, tolerance = 1e-1)

  expected_elements <- c("Likely_CellType",pbmc.single$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
         failure_message = "modscore results not found")

})

test_that("modScore returns metadata with scores and cell calls for
          nsclc.multi", {

  nsclc.multi <- getModuleScoreParam("nsclc.multi")

  suppressWarnings(modscore.demo <- do.call(modScore, nsclc.multi))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "nsclc_multi_MS_1.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "nsclc_multi_MS_2.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "nsclc_multi_MS_3.png"
  )

  expect_equal(mean(modscore.demo$ms.object$rand_type1),
               0.1266389, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$rand_type2),
               0.1794167, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$rand_type3),
               0.2627949, tolerance = 1e-1)

  expected_elements <- c("Likely_CellType",nsclc.multi$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
         failure_message = "modscore results not found")

})

test_that("modScore returns metadata with scores and cell calls for brca", {

  brca <- getModuleScoreParam("brca")

  suppressWarnings(modscore.demo <- do.call(modScore, brca))

  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[1]]),
    "brca_multi_MS_1.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[2]]),
    "brca_multi_MS_2.png"
  )
  
  skip_on_ci()
  expect_snapshot_file(
    .drawMSfig(modscore.demo$ms.figures[[3]]),
    "brca_multi_MS_3.png"
  )

  expect_equal(mean(modscore.demo$ms.object$rand_type1),
               0.1958383, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$rand_type2),
               0.1079014, tolerance = 1e-1)

  expect_equal(mean(modscore.demo$ms.object$rand_type3),
               0.03652896, tolerance = 1e-1)

  expected_elements <- c("Likely_CellType",brca$celltypes)
  expect(all(expected_elements %in% colnames(modscore.demo$ms.object@meta.data)),
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
