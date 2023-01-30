

# tests
test_that("Function returns a list with specified names", {
  
  # load data
  seurat.object <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  identities_match <- read.csv(test_path("fixtures", "ClusterNames_match.csv"))
  
  output <-
    NameClusters(
      SO = seurat.object,
      cluster.identities.table = identities_match,
      cluster.column.from.SO = "SCT_snn_res_0_2",
      cluster.names = "Cluster_Names",
      cluster.numbers = "Cluster_Numbers"
    )
  
  expect_type(output,"list")
  
  expected.elements = c("output", "plot")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})


test_that("Function returns correct class", {
  
  # load data
  seurat.object <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  identities_match <- read.csv(test_path("fixtures", "ClusterNames_match.csv"))
  identities_diff <- read.csv(test_path("fixtures", "ClusterNames_diff.csv"))
  identities_oneMore <- read.csv(test_path("fixtures", "ClusterNames_oneMore.csv"))
  
  output <-
    NameClusters(
      SO = seurat.object,
      #metadata = metadata,
      cluster.identities.table = identities_diff,
      cluster.column.from.SO = "SCT_snn_res_0_2",
      cluster.names = "Cluster_Names",
      cluster.numbers = "Cluster_Numbers"
    )
  expect_s3_class(output$output, "data.frame")
  
  output <-
    NameClusters(
      SO = seurat.object,
      cluster.identities.table = identities_oneMore,
      cluster.column.from.SO = "SCT_snn_res_0_2",
      cluster.names = "Cluster_Names",
      cluster.numbers = "Cluster_Numbers"
    )
  expect_s3_class(output$output, "data.frame")
  
  output <-
    NameClusters(
      SO = seurat.object,
      cluster.identities.table = identities_match,
      cluster.column.from.SO = "SCT_snn_res_0_2",
      cluster.names = "Cluster_Names",
      cluster.numbers = "Cluster_Numbers"
    )
  expect_s4_class(output$output, "Seurat")
  expect_equal(class(output$plot), c("plotly", "htmlwidget"))
  
})

test_that("Function returns Clusternames column", {
  
  # load data
  seurat.object <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  identities_match <- read.csv(test_path("fixtures", "ClusterNames_match.csv"))
 
  colname = "Clusternames"
  
  output <-
    NameClusters(
      SO = seurat.object,
      cluster.identities.table = identities_match,
      cluster.column.from.SO = "SCT_snn_res_0_2",
      cluster.names = "Cluster_Names",
      cluster.numbers = "Cluster_Numbers"
    )
  
  expect_named(output$output@meta.data[, intersect(colname, colnames(output$output@meta.data)), drop=FALSE])
  
})


