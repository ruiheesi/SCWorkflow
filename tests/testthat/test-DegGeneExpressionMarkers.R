test_that("Test DEG Gene Expression Markers using TEC (Mouse) dataset", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  DEG.result <- degGeneExpressionMarkers(object = obj,
                                         samples = 'c("1_Embryo_13_5","2_Embryo_15")',
                                         parameter.to.test = "SCT_snn_res_0_2",
                                         contrasts =  c("1-2","2-all"),
                                         test.to.use = "MAST",
                                         log.fc.threshold = 0,
                                         use.spark = FALSE,
                                         assay.to.use = "SCT",
                                         latent.vars = c())
  expected.elements = c("out_df")
  #    write.table(DEG.result[1],file="DEG.txt", quote=FALSE, sep="\t", row.names=FALSE)
  expect_setequal(names(DEG.result), expected.elements)
})





test_that("Test DEG Gene Expression Markers using negbinom (TEC Mouse dataset)", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  DEG.result <- degGeneExpressionMarkers(object = obj,
                                         samples = 'c("1_Embryo_13_5","2_Embryo_15")',
                                         parameter.to.test = "SCT_snn_res_0_2",
                                         contrasts =  c("1-2","2-all"),
                                         test.to.use = "negbinom",
                                         log.fc.threshold = 0,
                                         use.spark = FALSE,
                                         assay.to.use = "SCT",
                                         latent.vars = c())
  expected.elements = c("out_df")
  #    write.table(DEG.result[1],file="DEG.txt", quote=FALSE, sep="\t", row.names=FALSE)
  expect_setequal(names(DEG.result), expected.elements)
})



test_that("Test DEG Gene Expression Markers using Chariou (Mouse) dataset", {
  obj <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
  DEG.result <- degGeneExpressionMarkers(object = obj,
                                         samples = 'c("CD8dep","Combo","ENT","NHSIL12","PBS")',
                                         parameter.to.test = "SCT_snn_res_2_8",
                                         contrasts =  c("0-1","0-all"),
                                         test.to.use = "MAST",
                                         log.fc.threshold = 0,
                                         use.spark = FALSE,
                                         assay.to.use = "SCT",
                                         latent.vars = c())
  expected.elements = c("out_df")
  #    write.table(DEG.result[1],file="DEG.txt", quote=FALSE, sep="\t", row.names=FALSE)
  expect_setequal(names(DEG.result), expected.elements)
})




test_that("Test DEG Gene Expression Markers using BRCA (Human) dataset", {
  obj <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
  DEG.result <- degGeneExpressionMarkers(object = obj,
                                         samples = 'c("CID3586","CID3838","CID3921","CID3941","CID3946","CID3948","CID3963","CID4040","CID4066","CID4067","CID4290A","CID4398","CID44041","CID4461","CID4463","CID4465","CID4471","CID4495","CID44971","CID44991","CID4513","CID4515","CID45171","CID4523","CID4530N","CID4535")',
                                         parameter.to.test = "SCT_snn_res_0_2",
                                         contrasts =  c("0-1","0-all"),
                                         test.to.use = "MAST",
                                         log.fc.threshold = 0,
                                         use.spark = FALSE,
                                         assay.to.use = "SCT",
                                         latent.vars = c())
  expected.elements = c("out_df")
  #    write.table(DEG.result[1],file="DEG.txt", quote=FALSE, sep="\t", row.names=FALSE)
  expect_setequal(names(DEG.result), expected.elements)
})



test_that("Test DEG Gene Expression Markers using NSCLCmulti (Human) dataset", {
  obj <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  DEG.result <- degGeneExpressionMarkers(object = obj,
                                         samples = 'c("Donor_1","Donor_2","Donor_3","Donor_4","Donor_5","Donor_6","Donor_7")',
                                         parameter.to.test = "SCT_snn_res_0_2",
                                         contrasts =  c("0-1","0-all"),
                                         test.to.use = "MAST",
                                         log.fc.threshold = 0,
                                         use.spark = FALSE,
                                         assay.to.use = "SCT",
                                         latent.vars = c())
  expected.elements = c("out_df")
  #    write.table(DEG.result[1],file="DEG.txt", quote=FALSE, sep="\t", row.names=FALSE)
  expect_setequal(names(DEG.result), expected.elements)
})



test_that("Test DEG Gene Expression Markers using NSCLCsingle (Human) dataset", {
  obj <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
  DEG.result <- degGeneExpressionMarkers(object = obj,
                                         samples = 'c("NSCLC_Single")',
                                         parameter.to.test = "SCT_snn_res_0_2",
                                         contrasts =  c("0-1","0-all"),
                                         test.to.use = "MAST",
                                         ##                       test_to_use = "negbinom",
                                         log.fc.threshold = 0,
                                         use.spark = FALSE,
                                         assay.to.use = "SCT",
##                                         latent.vars = c("orig_ident"))
                                         latent.vars = c())
  expected.elements = c("out_df")
  #    write.table(DEG.result[1],file="DEG.txt", quote=FALSE, sep="\t", row.names=FALSE)
  expect_setequal(names(DEG.result), expected.elements)
  
})