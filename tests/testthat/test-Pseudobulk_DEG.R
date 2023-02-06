test_that("Pseudobulk DEG returns table of differentially expressed genes for TEC data", {
  
  TEC <- getparam_Pseudobulk("TEC")
  
  pseudo_res <- Pseudobulk_DEG(so = TEC$object, contrasts = TEC$contrasts, replicate = TEC$replicate, 
                               comparison_level = TEC$comparison_level,
                               label = TEC$label)
  
  expected_elements = c("Gene","cell_type","groupF-groupM_FC","groupF-groupM_logFC",
                        "groupF-groupM_tstat","groupF-groupM_pval",
                        "groupF-groupM_adjpval","de_family","de_method","de_type")
  expect_setequal(colnames(pseudo_res), expected_elements)
  
})

test_that("Pseudobulk DEG returns table of differentially expressed genes for Chariou data", {
  
  Chariou <- getparam_Pseudobulk("Chariou")
  
  pseudo_res <- Pseudobulk_DEG(so = Chariou$object, contrasts = Chariou$contrasts, replicate = Chariou$replicate, 
                               comparison_level = Chariou$comparison_level,
                               label = Chariou$label)
  
  expected_elements = c("Gene","cell_type","groupF-groupM_FC","groupF-groupM_logFC",
                        "groupF-groupM_tstat","groupF-groupM_pval",
                        "groupF-groupM_adjpval","de_family","de_method","de_type")
  expect_setequal(colnames(pseudo_res), expected_elements)
  
})

test_that("Pseudobulk DEG returns table of differentially expressed genes for NSCLC_Single data", {
  
  NSCLC_Single <- getparam_Pseudobulk("NSCLC_Single")
  
  pseudo_res <- Pseudobulk_DEG(so = NSCLC_Single$object, contrasts = NSCLC_Single$contrasts, replicate = NSCLC_Single$replicate, 
                               comparison_level = NSCLC_Single$comparison_level,
                               label = NSCLC_Single$label)
  
  expected_elements = c("Gene","cell_type","groupcluster_1-groupcluster_2_FC",
                        "groupcluster_1-groupcluster_2_logFC","groupcluster_1-groupcluster_2_tstat","groupcluster_1-groupcluster_2_pval",
                        "groupcluster_1-groupcluster_2_adjpval","de_family","de_method","de_type")
  expect_setequal(colnames(pseudo_res), expected_elements)
  
})

test_that("Pseudobulk DEG returns table of differentially expressed genes for NSCLC_Multi data", {
  
  NSCLC_Multi <- getparam_Pseudobulk("NSCLC_Multi")
  
  pseudo_res <- Pseudobulk_DEG(so = NSCLC_Multi$object, contrasts = NSCLC_Multi$contrasts, replicate = NSCLC_Multi$replicate, 
                               comparison_level = NSCLC_Multi$comparison_level,
                               label = NSCLC_Multi$label)
  
  expected_elements = c("Gene","cell_type","groupcluster_1-groupcluster_2_FC",
                        "groupcluster_1-groupcluster_2_logFC","groupcluster_1-groupcluster_2_tstat","groupcluster_1-groupcluster_2_pval",
                        "groupcluster_1-groupcluster_2_adjpval","de_family","de_method","de_type")
  expect_setequal(colnames(pseudo_res), expected_elements)
  
})

test_that("Pseudobulk DEG returns table of differentially expressed genes for BRCA data", {
  
  BRCA <- getparam_Pseudobulk("BRCA")
  
  pseudo_res <- Pseudobulk_DEG(so = BRCA$object, contrasts = BRCA$contrasts, replicate = BRCA$replicate, 
                               comparison_level = BRCA$comparison_level,
                               label = BRCA$label)
  
  expected_elements = c("Gene","cell_type","groupG2M-groupS_FC","groupG2M-groupS_logFC","groupG2M-groupS_tstat",
                        "groupG2M-groupS_pval","groupG2M-groupS_adjpval","de_family","de_method","de_type")
  expect_setequal(colnames(pseudo_res), expected_elements)
  
})

## Error Testing ##

test_that("Pseudobulk DEG stops when invalid contrast is selected", {
  
  NSCLC_Single <- getparam_Pseudobulk("NSCLC_Single")
  
  expect_error(Pseudobulk_DEG(so = NSCLC_Single$object, contrasts = c("wrong_contrast1-wrong_contrast2"), 
                              replicate = NSCLC_Single$replicate, 
                               comparison_level = NSCLC_Single$comparison_level,
                               label = NSCLC_Single$label), "contrasts not found amongst <label> metadata column")
  
})
