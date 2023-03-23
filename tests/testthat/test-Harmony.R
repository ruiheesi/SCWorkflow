test_that("Harmony returns seurat object with adjusted embeddings for
          TEC data", {

  tec = getHarmonyParam("TEC")

  object.harmonized = do.call(harmonyBatchCorrect, tec)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for Chariou
          data", {

  chariou = getHarmonyParam("Chariou")

  object.harmonized = do.call(harmonyBatchCorrect, chariou)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          nsclc_single data", {

  nsclc.single = getHarmonyParam("nsclc_single")

  object.harmonized = do.call(harmonyBatchCorrect, nsclc.single)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          nsclc_multi data", {

  nsclc.multi = getHarmonyParam("nsclc_multi")

  object.harmonized = do.call(harmonyBatchCorrect, nsclc.multi)

  expected_elements = c("adj.object","adj.tsne")
  expect_setequal(names(object.harmonized), expected_elements)
})

test_that("Harmony returns seurat object with adjusted embeddings for
          BRCA data", {

  brca = getHarmonyParam("BRCA")

  object.harmonized = do.call(harmonyBatchCorrect, brca)

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