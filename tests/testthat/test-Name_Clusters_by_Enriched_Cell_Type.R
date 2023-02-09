test_that("Run Name clusters with default parameters - TEC data", {
  
  # load data
  input <- getparams_nameclus("TEC") 
  output <- do.call(nameClusters,input)
  
  #Test output values and plot:
  newclus <- output$object@meta.data$clusternames
  expect_equal(sort(unique(newclus)),sort(input$cluster.names))
  
  ggsave("output/TEC_clusters.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","TEC_clusters.png")
  
  expect_type(output,"list")
  expected.elements = c("object","table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_s3_class(output$table, "gtable")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})

test_that("Run Name clusters with interactive plot", {
  
  # load data
  input <- getparams_nameclus("TEC") 
  input$interactive = TRUE
  output <- do.call(nameClusters,input)

  expect_equal(class(output$plot), c("plotly", "htmlwidget"))
  saveWidget(as_widget(output$plot),"output/TEC_clusters_plotly.html")
  expect_snapshot_file("output","TEC_clusters_plotly.html")
})

test_that("Run Name clusters with ordering celltypes", {
  
  input <- getparams_nameclus("TEC") 
  input$order.celltypes.by = c("B cells","Dendritic cells",
              "Endothelial cells","Macrophages","Monocytes",
              "Epithelial cells","Erythrocytes","Fibroblasts",
              "Hepatocytes","Neurons","T cells")
  
  output <- do.call(nameClusters,input)
  ggsave("output/TEC_clusters_ordered.png",output$plot,width = 10, height = 10)
  expect_snapshot_file("output","TEC_clusters_ordered.png")

})

test_that("Run Name clusters with ordering warning missing some celltypes", {
  
  input <- getparams_nameclus("TEC") 
  input$order.celltypes.by = c("B cells","Dendritic cells",
                               "Endothelial cells","Macrophages",
                               "Monocytes","Epithelial cells",
                               #"Erythrocytes","Fibroblasts",
                               "Hepatocytes","Neurons","T cells")
  
  expect_warning(output <- do.call(nameClusters,input),
      "^Some factors were not included in the list")
  #output <- do.call(nameClusters,input)
  ggsave("output/TEC_clusters_missing.png",output$plot,width = 10, height = 10)
  expect_snapshot_file("output","TEC_clusters_missing.png")
  
})

test_that("Run Name clusters with ordering warning adding some unknown celltypes", {
  
  input <- getparams_nameclus("TEC") 
  input$order.celltypes.by = c("B cells","Dendritic cells",
                               "Endothelial cells","Macrophages",
                               "Monocytes","Epithelial cells",
                               "celltype1","celltype2",
                               "Erythrocytes","Fibroblasts",
                               "Hepatocytes","Neurons","T cells")
  
  expect_warning(output <- do.call(nameClusters,input),
                 "^Some factors are not in data")
  #output <- do.call(nameClusters,input)
  ggsave("output/TEC_clusters_missing2.png",output$plot,width = 10, height = 10)
  expect_snapshot_file("output","TEC_clusters_missing2.png")
  
})

test_that("Run Name clusters with default parameters - Chariou", {
  
  # load data
  input <- getparams_nameclus("Chariou") 
  expect_warning(output <- do.call(nameClusters,input),
      "^Some clusters had no detected cell types")
  
  #Test output values and plot:
  newclus <- output$object@meta.data$clusternames
  expect_equal(sort(unique(newclus)),sort(input$cluster.names))
  
  ggsave("output/Chariou_clusters.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","Chariou_clusters.png")
  
  expect_type(output,"list")
  expected.elements = c("object","table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_s3_class(output$table, "gtable")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})

test_that("Run Name clusters with default parameters - NSCLC single", {
  
  # load data
  input <- getparams_nameclus("nsclc-single") 
  output <- do.call(nameClusters,input)
  
  #Test output values and plot:
  newclus <- output$object@meta.data$clusternames
  expect_equal(sort(unique(newclus)),sort(input$cluster.names))
  
  ggsave("output/NSCLC_single_clusters.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","NSCLC_single_clusters.png")
  
  expect_type(output,"list")
  expected.elements = c("object","table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_s3_class(output$table, "gtable")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
})

test_that("Run Name clusters with default parameters - NSCLC multi", {
  
  # load data
  input <- getparams_nameclus("nsclc-multi") 
  output <- do.call(nameClusters,input)
  
  #Test output values and plot:
  newclus <- output$object@meta.data$clusternames
  expect_equal(sort(unique(newclus)),sort(input$cluster.names))
  
  ggsave("output/NSCLC_multi_clusters.png",output$plot, width = 10, height = 10)
  expect_snapshot_file("output","NSCLC_multi_clusters.png")
  
  expect_type(output,"list")
  expected.elements = c("object","table", "plot")
  expect_s4_class(output$object, "Seurat")
  expect_s3_class(output$table, "gtable")
  expect_equal(length(setdiff(expected.elements, names(output))), 0)
  
}

)

