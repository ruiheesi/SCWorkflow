test_that("Test Plot Metadata using TEC (Mouse) dataset", {
  
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  PlotMeta.result <- plotMetadata(object = obj,
                                  samples.to.include = 'c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")',
                                  metadata.to.plot = 'c("SCT_snn_res_0_2","SCT_snn_res_0_4","SCT_snn_res_0_6","SCT_snn_res_0_8","SCT_snn_res_1","SCT_snn_res_1_2")',
                                  reduction.type = "umap",
                                  use.cite.seq = FALSE,
                                  number.of.columns.for.final.image <- 0,
                                  show.labels = FALSE,
                                  legend.text.size = 1,
                                  legend.position = "right",
                                  columns.to.summarize = "c()",
                                  summarization.cut.off = 5,
                                  dot.size = 0.01)
  
  
#  ggsave("output/TEC_PlotMetadata_umap.pdf",PlotMeta.result$plot, width = 10, height = 10)
#  expect_snapshot_file("output","TEC_PlotMetadata_umap.pdf")
  expect_type(PlotMeta.result,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(PlotMeta.result), expected.elements) 
  
})





test_that("Test Plot Metadata using  using TSNE (TEC Mouse dataset)", {
  
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  PlotMeta.result <- plotMetadata(object = obj,
                                  samples.to.include = 'c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")',
                                  metadata.to.plot = 'c("SCT_snn_res_0_2","SCT_snn_res_0_4","SCT_snn_res_0_6","SCT_snn_res_0_8","SCT_snn_res_1","SCT_snn_res_1_2")',
                                  reduction.type = "tsne",
                                  use.cite.seq = FALSE,
                                  number.of.columns.for.final.image <- 0,
                                  show.labels = FALSE,
                                  legend.text.size = 1,
                                  legend.position = "right",
                                  columns.to.summarize = "c()",
                                  summarization.cut.off = 5,
                                  dot.size = 0.01)

#  ggsave("output/TEC_PlotMetadata_tsne.pdf",PlotMeta.result$plot, width = 10, height = 10)
#  expect_snapshot_file("output","TEC_PlotMetadata_tsne.pdf")
  expect_type(PlotMeta.result,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(PlotMeta.result), expected.elements) 
  
})













test_that("Test Plot Metadata using Chariou (Mouse) dataset", {
  
  obj <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
  PlotMeta.result <- plotMetadata(object = obj,
                                  samples.to.include = 'c("CD8dep","Combo","ENT","NHSIL12","PBS")',
                                  metadata.to.plot = 'c("SCT_snn_res_2_4","SCT_snn_res_2_6","SCT_snn_res_2_8")',
                                  reduction.type = "umap",
                                  use.cite.seq = FALSE,
                                  number.of.columns.for.final.image <- 0,
                                  show.labels = FALSE,
                                  legend.text.size = 1,
                                  legend.position = "right",
                                  columns.to.summarize = "c()",
                                  summarization.cut.off = 5,
                                  dot.size = 0.01)
  
  
#  ggsave("output/Chariou_PlotMetadata_umap.pdf",PlotMeta.result$plot, width = 10, height = 10)
#  expect_snapshot_file("output","Chariou_PlotMetadata_umap.pdf")
  expect_type(PlotMeta.result,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(PlotMeta.result), expected.elements) 
  
})






test_that("Test Plot Metadata using BRCA (Human) dataset", {
  
  obj <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
  PlotMeta.result <- plotMetadata(object = obj,
                                  samples.to.include = 'c("CID3586","CID3838","CID3921","CID3941","CID3946","CID3948","CID3963","CID4040","CID4066","CID4067","CID4290A","CID4398","CID44041","CID4461","CID4463","CID4465","CID4471","CID4495","CID44971","CID44991","CID4513","CID4515","CID45171","CID4523","CID4530N","CID4535")',
                                  metadata.to.plot = 'c("SCT_snn_res_0_2","SCT_snn_res_0_4","SCT_snn_res_0_6","SCT_snn_res_0_8","SCT_snn_res_1","SCT_snn_res_1_2")',
                                  reduction.type = "tsne",
                                  use.cite.seq = FALSE,
                                  number.of.columns.for.final.image <- 0,
                                  show.labels = FALSE,
                                  legend.text.size = 1,
                                  legend.position = "right",
                                  columns.to.summarize = "c()",
                                  summarization.cut.off = 5,
                                  dot.size = 0.01)
  
  
#  ggsave("output/BRCA_PlotMetadata_tsne.pdf",PlotMeta.result$plot, width = 10, height = 10)
#  expect_snapshot_file("output","BRCA_PlotMetadata_tsne.pdf")
  expect_type(PlotMeta.result,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(PlotMeta.result), expected.elements) 
  
})




test_that("Test Plot Metadata using NSCLCmulti (Human) dataset", {
  
  obj <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  PlotMeta.result <- plotMetadata(object = obj,
                                  samples.to.include = 'c("Donor_1","Donor_2","Donor_3","Donor_4","Donor_5","Donor_6","Donor_7")',
                                  metadata.to.plot = 'c("SCT_snn_res_0_2","SCT_snn_res_0_4","SCT_snn_res_0_6","SCT_snn_res_0_8","SCT_snn_res_1","SCT_snn_res_1_2")',
                                  reduction.type = "tsne",
                                  use.cite.seq = FALSE,
                                  number.of.columns.for.final.image <- 0,
                                  show.labels = FALSE,
                                  legend.text.size = 1,
                                  legend.position = "right",
                                  columns.to.summarize = "c()",
                                  summarization.cut.off = 5,
                                  dot.size = 0.01)
  
  
#  ggsave("output/NSCLCmulti_PlotMetadata_tsne.pdf",PlotMeta.result$plot, width = 10, height = 10)
#  expect_snapshot_file("output","NSCLCmulti_PlotMetadata_tsne.pdf")
  expect_type(PlotMeta.result,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(PlotMeta.result), expected.elements) 
  
})



test_that("Test Plot Metadata using NSCLCsingle (Human) dataset", {
  
  obj <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
  PlotMeta.result <- plotMetadata(object = obj,
                                  samples.to.include = 'c("NSCLC_Single")',
                                  metadata.to.plot = 'c("SCT_snn_res_0_2","SCT_snn_res_0_4","SCT_snn_res_0_6","SCT_snn_res_0_8","SCT_snn_res_1","SCT_snn_res_1_2")',
                                  reduction.type = "tsne",
                                  use.cite.seq = FALSE,
                                  number.of.columns.for.final.image <- 0,
                                  show.labels = FALSE,
                                  legend.text.size = 1,
                                  legend.position = "right",
                                  columns.to.summarize = "c()",
                                  summarization.cut.off = 5,
                                  dot.size = 0.01)
  
  
#  ggsave("output/NSCLCsingle_PlotMetadata_tsne.pdf",PlotMeta.result$plot, width = 10, height = 10)
#  expect_snapshot_file("output","NSCLCsingle_PlotMetadata_tsne.pdf")
  expect_type(PlotMeta.result,"list")
  expected.elements = c("object", "plot")
  expect_setequal(names(PlotMeta.result), expected.elements) 
  
})