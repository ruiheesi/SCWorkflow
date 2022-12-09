test_that("Pseudobulk DEG returns table of differentially expressed genes", {
  so_demo <- readRDS(test_path("fixtures", "SO_moduleScore.rds"))
  
  pseudo_res <- Pseudobulk_DEG(so = so_demo, contrasts = c("CD4_T-Tregs"), replicate = 'orig.ident', comparison_level = 'Likely_CellType',
                 label = "Likely_CellType")
  
  expected_elements = c("Gene","cell_type","groupCD4_T-groupTregs_FC","groupCD4_T-groupTregs_logFC",
                        "groupCD4_T-groupTregs_tstat","groupCD4_T-groupTregs_pval",
                        "groupCD4_T-groupTregs_adjpval","de_family","de_method","de_type")
  expect_setequal(colnames(pseudo_res), expected_elements)
})
