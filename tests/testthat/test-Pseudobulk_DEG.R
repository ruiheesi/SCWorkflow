# Data Testing
test_that("Pseudobulk DEG returns table of differentially expressed genes for 
          TEC data", {
            
            tec = getPseudobulkParam("TEC")
            
            pseudo.res = do.call(pseudobulkDEG, tec)
            
            expected_elements = c("Gene",
                                  "subgroup",
                                  "groupF-groupM_FC",
                                  "groupF-groupM_logFC",
                                  "groupF-groupM_tstat",
                                  "groupF-groupM_pval",
                                  "groupF-groupM_adjpval",
                                  "de.family",
                                  "de.method",
                                  "de.type")
            
            expect_setequal(colnames(pseudo.res), expected_elements)
            
          })

test_that("Pseudobulk DEG returns table of differentially expressed genes for
          Chariou data", {
            
            chariou = getPseudobulkParam("Chariou")
            
            pseudo.res = do.call(pseudobulkDEG, chariou)
            
            expected_elements = c("Gene",
                                  "subgroup",
                                  "groupF-groupM_FC",
                                  "groupF-groupM_logFC",
                                  "groupF-groupM_tstat",
                                  "groupF-groupM_pval",
                                  "groupF-groupM_adjpval",
                                  "de.family",
                                  "de.method",
                                  "de.type")
            
            expect_setequal(colnames(pseudo.res), expected_elements)
            
          })

test_that("Pseudobulk DEG returns table of differentially expressed genes for
          pbmc.single data", {
            
            pbmc.single = getPseudobulkParam("pbmc.single")
            
            pseudo.res = do.call(pseudobulkDEG, pbmc.single)
            
            expected_elements = c("Gene",
                                  "subgroup",
                                  "groupcluster_1-groupcluster_2_FC",
                                  "groupcluster_1-groupcluster_2_logFC",
                                  "groupcluster_1-groupcluster_2_tstat",
                                  "groupcluster_1-groupcluster_2_pval",
                                  "groupcluster_1-groupcluster_2_adjpval",
                                  "de.family",
                                  "de.method",
                                  "de.type")
            
            expect_setequal(colnames(pseudo.res), expected_elements)
            
          })

test_that("Pseudobulk DEG returns table of differentially expressed genes for
          nsclc.multi data", {
            
            nsclc.multi = getPseudobulkParam("nsclc.multi")
            
            pseudo.res = do.call(pseudobulkDEG, nsclc.multi)
            
            expected_elements = c("Gene",
                                  "subgroup",
                                  "groupcluster_1-groupcluster_2_FC",
                                  "groupcluster_1-groupcluster_2_logFC",
                                  "groupcluster_1-groupcluster_2_tstat",
                                  "groupcluster_1-groupcluster_2_pval",
                                  "groupcluster_1-groupcluster_2_adjpval",
                                  "de.family",
                                  "de.method",
                                  "de.type")
            
            expect_setequal(colnames(pseudo.res), expected_elements)
            
          })

test_that("Pseudobulk DEG returns table of differentially expressed genes for
          BRCA data", {
            
            BRCA = getPseudobulkParam("BRCA")
            
            pseudo.res = do.call(pseudobulkDEG, BRCA)
            
            expected_elements = c("Gene",
                                  "subgroup",
                                  "groupG2M-groupS_FC",
                                  "groupG2M-groupS_logFC",
                                  "groupG2M-groupS_tstat",
                                  "groupG2M-groupS_pval",
                                  "groupG2M-groupS_adjpval",
                                  "de.family",
                                  "de.method",
                                  "de.type")
            
            expect_setequal(colnames(pseudo.res), expected_elements)
            
          })

# Error Testing
test_that("Pseudobulk DEG stops when invalid contrast is selected", {
  
  pbmc.single = getPseudobulkParam("pbmc.single")
  
  expect_error(pseudobulkDEG(object = pbmc.single$object,
                             contrasts = c("wrong_contrast1-wrong_contrast2"),
                             replicate = pbmc.single$replicate,
                             subgroup = pbmc.single$subgroup,
                             group = pbmc.single$group),
               "contrasts not found amongst <group> metadata column")
  
})
