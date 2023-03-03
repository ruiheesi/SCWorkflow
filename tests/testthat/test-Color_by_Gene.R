test_that("Color by Gene using TEC (Mouse) dataset", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  ColorByGene.result <- colorByGene(object = obj,
                                    samples.to.include = 'c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")',
                                    gene = c("Gapdh","GAPDH","Foxp3","Il2ra","Cd4","Cd8a"))
  
  expected.elements = c("object", "plot")
  expect_setequal(names(ColorByGene.result), expected.elements) 
})


test_that("Color by Gene for TSNE plot (Mouse TEC dataset)", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  ColorByGene.result <- colorByGene(object = obj,
                             samples.to.include = 'c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")',
                             gene = c("Gapdh","GAPDH","Foxp3","Il2ra","Cd4","Cd8a"),
                             reduction.type = "tsne")

  expected.elements = c("object", "plot")
  expect_setequal(names(ColorByGene.result), expected.elements) 
})

test_that("Color by Gene with blue color (Mouse TEC dataset)", {
  obj <- readRDS(test_path("fixtures/TEC", "TEC_Combine_and_Renormalize_SO_downsample.rds"))
  ColorByGene.result <- colorByGene(object = obj,
                                    samples.to.include = 'c("1_Embryo_13_5","2_Embryo_15","3_Newborn","4_Adult")',
                                    gene = c("Gapdh","GAPDH","Foxp3","Il2ra","Cd4","Cd8a"),
                                    reduction.type = "tsne",
                                    color = "blue")
  
  expected.elements = c("object", "plot")
  expect_setequal(names(ColorByGene.result), expected.elements) 
})



test_that("Color by Gene using Chariou (Mouse) dataset", {    
  obj <- readRDS(test_path("fixtures/Chariou", "Chariou_Combine_and_Renormalize_SO_downsample.rds"))
  ColorByGene.result <- colorByGene(object = obj,
                                    samples.to.include = 'c("CD8dep","Combo","ENT","NHSIL12","PBS")',
                                    gene = c("Gapdh","GAPDH","Foxp3","Il2ra","Cd4","Cd8a"))
  
  expected.elements = c("object", "plot")
  expect_setequal(names(ColorByGene.result), expected.elements) 
})




test_that("Color by Gene using BRCA (Human) dataset", {    
  obj <- readRDS(test_path("fixtures/BRCA", "BRCA_Combine_and_Renormalize_SO_downsample.rds"))
  ColorByGene.result <- colorByGene(object = obj,
                                    samples.to.include = 'c("CID3586","CID3838","CID3921","CID3941","CID3946","CID3948","CID3963","CID4040","CID4066","CID4067","CID4290A","CID4398","CID44041","CID4461","CID4463","CID4465","CID4471","CID4495","CID44971","CID44991","CID4513","CID4515","CID45171","CID4523","CID4530N","CID4535")',
                                    gene = c("Gapdh","GAPDH","FOXP3","IL2RA","CD4","CD8A","MYC", "PARP8"),
                                    reduction.type = "tsne")
  
  expected.elements = c("object", "plot")
  expect_setequal(names(ColorByGene.result), expected.elements) 
})




test_that("Test Annotate Cell Types using NSCLCmulti (Human) dataset", {    
  obj <- readRDS(test_path("fixtures/NSCLC_Multi", "NSCLCmulti_Combine_and_Renormalize_SO_downsample.rds"))
  ColorByGene.result <- colorByGene(object = obj,
                                    samples.to.include = 'c("Donor_1","Donor_2","Donor_3","Donor_4","Donor_5","Donor_6","Donor_7")',
                                    gene = c("Gapdh","GAPDH","FOXP3","IL2RA","CD4","CD8A","MYC","AMER2","CD9","CEBPD"),
                                    reduction.type = "tsne")
  
  expected.elements = c("object", "plot")
  expect_setequal(names(ColorByGene.result), expected.elements) 
})



test_that("Test Annotate Cell Types using NSCLCsingle (Human) dataset", {    
  obj <- readRDS(test_path("fixtures/NSCLC_Single", "NSCLCsingle_Combine_and_Renormalize_SO_downsample.rds"))
  ColorByGene.result <- colorByGene(object = obj,
                                    samples.to.include = 'c("NSCLC_Single")',
                                    gene = c("Gapdh","GAPDH","FOXP3","IL2RA","CD4","CD8A","MYC"),
                                    reduction.type = "tsne")
  
  expected.elements = c("object", "plot")
  expect_setequal(names(ColorByGene.result), expected.elements) 
})

