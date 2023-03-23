test_that("Run Name clusters with default parameters - TEC data", {
  input <- getParamsNameClus("TEC")
  output <- do.call(nameClusters, input)
  
  #Test output values and plot:
  newclus <- output$object@meta.data$clusternames
  expect_snapshot_file(ggsave(
    "TEC_clusters.png",
    output$plot,
    width = 10,
    height = 10
  ),
  "TEC_clusters.png")
  
  expect_type(output, "list")
  expected.elements = c("object", "table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})

test_that("Run Name clusters with interactive plot", {
  # load data
  input <- getParamsNameClus("TEC")
  input$interactive = TRUE
  output <- do.call(nameClusters, input)
  
  expect_equal(class(output$plot), c("plotly", "htmlwidget"))
  #expect_snapshot_file(plotly::save_image(output$plot,
  #                                        file = "TEC_clusters2.png"),
  #                     "TEC_clusters2.png")
})

test_that("Run Name clusters with ordering celltypes", {
  input <- getParamsNameClus("TEC")
  input$order.celltypes.by = c(
    "B cells",
    "Dendritic cells",
    "Endothelial cells",
    "Macrophages",
    "Monocytes",
    "Epithelial cells",
    "Erythrocytes",
    "Fibroblasts",
    "Hepatocytes",
    "Neurons",
    "T cells"
  )
  
  output <- do.call(nameClusters, input)
  expect_snapshot_file(
    ggsave(
      "TEC_clusters_ordered.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "TEC_clusters_ordered.png"
  )
  
})

test_that("Run Name clusters with ordering warning missing some celltypes", {
  input <- getParamsNameClus("TEC")
  input$order.celltypes.by = c(
    "B cells",
    "Dendritic cells",
    "Endothelial cells",
    "Macrophages",
    "Monocytes",
    "Epithelial cells",
    #"Erythrocytes","Fibroblasts",
    "Hepatocytes",
    "Neurons",
    "T cells"
  )
  
  expect_warning(output <- do.call(nameClusters, input),
                 "^Some factors were not included in the list")
  expect_snapshot_file(
    ggsave(
      "TEC_clusters_missing.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "TEC_clusters_missing.png"
  )
  
})

test_that("Run Name clusters with ordering warning adding unknown celltypes",
          {
            input <- getParamsNameClus("TEC")
            input$order.celltypes.by = c(
              "B cells",
              "Dendritic cells",
              "Endothelial cells",
              "Macrophages",
              "Monocytes",
              "Epithelial cells",
              "celltype1",
              "celltype2",
              "Erythrocytes",
              "Fibroblasts",
              "Hepatocytes",
              "Neurons",
              "T cells"
            )
            expect_warning(output <- do.call(nameClusters, input),
                           "^Some factors are not in data")
            expect_snapshot_file(
              ggsave(
                "TEC_clusters_missing2.png",
                output$plot,
                width = 10,
                height = 10
              ),
              "TEC_clusters_missing2.png"
            )
          })

test_that("Run Name clusters with default parameters - Chariou", {
  input <- getParamsNameClus("Chariou")
  output <- do.call(nameClusters, input)
  
  #Test output values and plot:
  newclus <- output$object@meta.data$Clusternames
  expect_equal(sort(unique(newclus)), sort(input$cluster.names))
  
  expect_snapshot_file(
    ggsave(
      "Chariou_clusters.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "Chariou_clusters.png"
  )
  
  expect_type(output, "list")
  expected.elements = c("object", "table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})

test_that("Run Name clusters with default parameters - PBMC single", {
  input <- getParamsNameClus("pbmc-single")
  output <- do.call(nameClusters, input)
  
  #Test output values and plot:
  newclus <- output$object@meta.data$Clusternames
  expect_equal(sort(unique(newclus)), sort(input$cluster.names))
  
  expect_snapshot_file(
    ggsave(
      "NSCLC_single_clusters.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "NSCLC_single_clusters.png"
  )
  
  expect_type(output, "list")
  expected.elements = c("object", "table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})

test_that("Run Name clusters with default parameters - NSCLC multi", {
  input <- getParamsNameClus("nsclc-multi")
  output <- do.call(nameClusters, input)
  
  #Test output values and plot:
  newclus <- output$object@meta.data$Clusternames
  expect_equal(sort(unique(newclus)), sort(input$cluster.names))
  
  expect_snapshot_file(
    ggsave(
      "NSCLC_multi_clusters.png",
      output$plot,
      width = 10,
      height = 10
    ),
    "NSCLC_multi_clusters.png"
  )
  
  expect_type(output, "list")
  expected.elements = c("object", "table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})
