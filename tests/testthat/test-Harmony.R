test_that("Harmony returns seurat object with adjusted embeddings for
          TEC data", {

  tec = getHarmonyParam("TEC")

  object.harmonized = do.call(harmonyBatchCorrect, tec)
   
  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$adj.tsne),
    "tec_harm.png"
  )
  
  expect_equal(mean(object.harmonized$adj.object[["harmSCT"]]["Lum",]),
               0.06279209, tolerance = 1e-1)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for Chariou
          data", {

  chariou = getHarmonyParam("Chariou")

  object.harmonized = do.call(harmonyBatchCorrect, chariou)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$adj.tsne),
    "char_harm.png"
  )

  expect_equal(mean(object.harmonized$adj.object[["harmSCT"]]["Ccl8",]),
               0.2046126, tolerance = 1e-1)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          pbmc_single data", {

  pbmc.single = getHarmonyParam("pbmc_single")

  object.harmonized = do.call(harmonyBatchCorrect, pbmc.single)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$adj.tsne),
    "pbmc_single_harm.png"
  )

  expect_equal(mean(object.harmonized$adj.object[["harmSCT"]]["CLU",]),
               -0.002671558, tolerance = 1e-1)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          nsclc_multi data", {

  nsclc.multi = getHarmonyParam("nsclc_multi")

  object.harmonized = do.call(harmonyBatchCorrect, nsclc.multi)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$adj.tsne),
    "nsclc_multi_harm.png"
  )

  expect_equal(mean(object.harmonized$adj.object[["harmSCT"]]["LCN2",]),
               0.1265227, tolerance = 1e-1)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          BRCA data", {

  brca = getHarmonyParam("BRCA")

  object.harmonized = do.call(harmonyBatchCorrect, brca)

  skip_on_ci()
  expect_snapshot_file(
    .drawHarmonyFig(object.harmonized$adj.tsne),
    "BRCA_harm.png"
  )

  expect_equal(mean(object.harmonized$adj.object[["harmSCT"]]["SCGB2A2",]),
               -0.1078249, tolerance = 1e-1)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony provides warning when genes are not found in the data", {

  tec = getHarmonyParam("TEC")

  expect_warning(harmonyBatchCorrect(
    object = tec$object,
    nvar = tec$nvar,
    genes.to.add = c("wrong_gene","wrong_gene2"),
    group.by.var = tec$group.by.var),
    "specified genes were not found and therefore cannot be added")

})

test_that("Harmony stops when variable features to subset by exceeds number of
          genes in the data", {

  tec = getHarmonyParam("TEC")

  expect_error(harmonyBatchCorrect(
    object = tec$object,
    nvar = 100000,
    genes.to.add = NULL,
    group.by.var = tec$group.by.var),
    "nvar exceed total number of variable genes in the data")

})