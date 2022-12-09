test_that("Test Dual labeling", {
  
  obj <- readRDS("/rstudio-files/ccbr-data/users/maggie/SCWorkflow/tests/testthat/fixtures/SO_moduleScore.rds")
  markertypes <- c("SCT","protein")
  reductions <- c("tsne","umap")
  genenames <- list(geneset2 <- c("badgene","Cd8a"),gene1 = c("Cd8a","Cd4"))
  for(g in genenames){
    for(m in markertypes){
      for(r in reductions){
        print(m)
        duallabel.result <- DualLabeling(object = obj,
                              samples <- c("1_E13","2_E15","3_Newborn","4_Adult"),
                              marker1 <- g[1],
                              marker_1_type = m,
                              marker2 <- g[2],
                              marker_2_type = m,
                              data_reduction = r,
                              density_heatmap = TRUE,
                              display_unscaled_values = TRUE)
        ggsave(file=paste0(g[1],"_",g[2],"_",m,"_",r,"_duallabel.pdf"),
                           duallabel.result$plot)
      }
    }
  }
  
  expected.elements = c("so","plot")
  expect_setequal(names(duallabel.result), expected.elements)
  
})
        