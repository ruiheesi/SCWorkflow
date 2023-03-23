# Original Code comes from https://github.com/neurorestore/Libra/blob/main/R/pseudobulk_de.R

#' @title Pseudobulk Differential Gene Expression Analysis
#' @description Performs Linear Models for Microarray Data on single-cell data 
#' @details Average gene expression of replicates before performing 
#'          limma regression
#' 
#' @param object Seurat-class object
#' @param contrasts Choose the groups (group1-group2) to compare transcription 
#'                  profiles against, positive Fold Change values indicate 
#'                  upregulation in group1
#' @param replicate Column of metadata containing replicate information, 
#'                  gene expression will first be aggregated based on this 
#'                  column before differential analysis is ran
#' @param subgroup Name of metadata column containing information on the 
#'                 sub-levels you want to divide your data into (optional, 
#'                 same as group if you don't want to subdivide)
#' @param group Name of metadata column containing information on the groups 
#'              (contrasts) you want to compare
#' @param min.num.cells Filters a subgroup level based on the minimum number of 
#'                      cells at that level (Default: 3)
#' @param min.num.reps Filters a subgroup level based on the minimum number of
#'                     replicates at that level (Default: 1)
#' @param min.num.features Filters out genes if they are expressed in fewer 
#'                         than X number of cells (Default: 0)
#' @param subdivide Toggle to TRUE if you want to make comparisons across 
#'                  different levels within your data (Default: FALSE)
#' @param external Analyze external datasets (Default: TRUE)
#' @param de.family General approach of running DEG (Default: pseudobulk)
#' @param de.method Method of running DEG (Default: limma)
#' @param de.type Method of creating design for DEG (Default: LRT)
#' 
#' @import Seurat
#' @import edgeR
#' @import purrr
#' @import magrittr
#' @import tidyverse
#' @import dplyr
#' @import limma
#' @import statmod
#'   
#' @export
#' @example Do not run: pseudobulkDEG(object = seurat,
#'                                    contrasts = "Drug-Vehicle",
#'                                    replicate = "Mouse_Treatment_Replicates",
#'                                    subgroup = "Mouse_Treatment_Group",
#'                                    group = "Mouse_Treatment_Group",
#'                                    subdivide = FALSE)
#'                                    
#' @example Do not run: pseudobulkDEG(object = seurat,
#'                                    contrasts = "Drug-Vehicle",
#'                                    replicate = "Mouse_Treatment_Replicates",
#'                                    subgroup = "Species",
#'                                    group = "Mouse_Treatment_Group",
#'                                    subdivide = TRUE)

#' @return data.frame of differentially expressed genes across groups and 
#'         (if specified) subgroups

pseudobulkDEG <- function(object,
                          contrasts,
                          replicate,
                          subgroup,
                          group,
                          min.num.cells = 3,
                          min.num.reps = 1,
                          min.num.features = 0,
                          subdivide = FALSE,
                          external = T,
                          de.family = "pseudobulk",
                          de.method = "limma",
                          de.type = "LRT") {
  
  ## Check Inputs ##
  # Error Messages
  check.contrasts <- unlist(strsplit(contrasts, "-"))
  if(sum(check.contrasts %in% object@meta.data[,eval(parse(text = "group"))]) 
     < 2){
    stop("contrasts not found amongst <group> metadata column")
  }
    
  # extract cell types and group from metadata
  if ("Seurat" %in% class(object)) {
    meta = object@meta.data %>%
      droplevels()
    if (!is.null(replicate))
      replicates = as.character(meta[[replicate]])
    if (!is.factor(meta[[group]])) {
      labels = meta[[group]]
    } else {
      labels = as.character(meta[[group]])
    }
    celltypes = as.character(meta[[subgroup]])
    expr = Seurat::GetAssayData(object, slot = 'counts')
  } else {
    # check if input is sparse matrix or numberic matrix/df
    valid.input = is(object, 'sparseMatrix') ||
      is_numeric_matrix(object) ||
      is_numeric_dataframe(object)
    if (!valid.input)
      stop("input must be Seurat, monocle, sparse matrix, numeric matrix, or ",
           "numeric data frame")
    if (is.null(meta))
      stop("input matrix must be accompanied by a metadata table")
    expr = object
    if (!is.null(replicate))
      replicates = as.character(meta[[replicate]])
    labels = as.character(meta[[group]])
    celltypes = as.character(meta[[subgroup]])
  }
  
  # check dimensions are non-zero
  if (length(dim(expr)) != 2 || !all(dim(expr) > 0)) {
    stop("expression matrix has at least one dimension of size zero")
  }
  
  # check dimensions match
  ncells1 = nrow(meta)
  ncells2 = ncol(expr)
  if (ncells1 != ncells2) {
    stop("number of cells in metadata (", ncells1, ") does not match number ",
         "of cells in expression (", ncells2, ")")
  }
  
  # check at least two labels
  if (n_distinct(labels) == 1) {
    stop("only one group provided: ", unique(labels))
  }
  
  # check for missing labels or cell types
  if (any(is.na(labels))) {
    stop("labels contain ", sum(is.na(labels)), "missing values")
  }
  if (any(is.na(celltypes))) {
    stop("cell types contain ", sum(is.na(celltypes)), "missing values")
  }
  if (!is.null(replicate) && any(is.na(replicates))) {
    stop("replicates contain ", sum(is.na(replicates)), "missing values")
  }
  
  # check for missing replicates
  if (!is.null(replicate) && is.null(replicates)) {
    stop("metadata does not contain replicate information")
  }
  
  # remove missing values
  missing = is.na(expr)
  if (any(missing)) {
    stop("matrix contains ", sum(missing), "missing values")
  }
  
  # clean up the meta data
  if (!is.null(replicate)) {
    meta <- meta %>% as.data.frame() %>%
      mutate(cell.barcode = rownames(meta),
             replicate = meta[[replicate]],
             subgroup = meta[[subgroup]],
             group = meta[[group]]) %>%
      mutate_at(vars(replicate, subgroup, group), as.factor)
  } else {
    meta <- meta %>% as.data.frame() %>%
      mutate(cell.barcode = rownames(meta),
             subgroup = meta[[subgroup]],
             group = meta[[group]]) %>%
      mutate_at(vars(subgroup, group), as.factor)
  }
  
  # make sure meta contains row names and is a data frame
  rownames(meta) = colnames(expr)
  meta = as.data.frame(meta)
  inputs = list(
    expr = expr,
    meta = meta
  )

## Create a pseudobulk matrix ##
# Convert a single-cell expression matrix (i.e., genes by cells)
#   to a pseudobulk matrix by summarizing counts within biological replicates
    expr = inputs$expr
    meta = inputs$meta
  
  # convert to characters
  meta$replicate <- as.character(meta[,eval(parse(text = "replicate"))])
  meta$subgroup <- as.character(meta[,eval(parse(text = "subgroup"))])
  meta$group <- as.character(meta[,eval(parse(text = "group"))])
  
  # keep only cell types with enough cells
  keep = meta %>%
    dplyr::count(subgroup, group) %>%
    group_by(subgroup) %>%
    dplyr::filter(all(n >= min.num.cells)) %>%
    pull(subgroup) %>%
    unique()
  
  # Need to make rownames = barcode for matrix multiplication
  rownames(meta) <- meta$Barcode
  
  ##JB set default to subgroup to NULL, meaning not to partition analysis ====
  #JB if (subgroup != NULL), then partition analysis ====
  if (subdivide = TRUE){
    # process data into gene x replicate x subgroup matrices
    pseudobulks = keep %>%
      map( ~ {
        print(.)
        subgroup = .
        meta0 = meta %>% filter(subgroup == !!subgroup)
        expr0 = expr %>% magrittr::extract(, meta0$cell.barcode)
        # catch cell types without replicates or conditions
        if (n_distinct(meta0$group) < 2)
          return(NA)
        replicate.counts = distinct(meta0, group, replicate) %>%
          group_by(group) %>%
          summarise(replicates = n_distinct(replicate)) %>%
          pull(replicates)
        if (any(replicate.counts < min.num.reps))
          return(NA)

        # mm: cell x replicate membership matrix (matrix of 1s and 0s)
        mm = model.matrix(~ 0 + replicate:group, data = meta0)
        
        # mat.mm 
        mat.mm = expr0 %*% mm
        
        keep.genes = rowSums(as.array(mat.mm > 0)) > min.num.features

        mat.mm = mat.mm[keep.genes, ] %>% as.data.frame()
        mat.mm %<>% as.data.frame()
        colnames(mat.mm) = gsub("replicate|group", "", colnames(mat.mm))
        # drop empty columns
        keep.samples = colSums(mat.mm) > 0
        mat.mm %<>% magrittr::extract(, keep.samples)
        return(mat.mm)
      }) %>% setNames(keep)
    } else { keep = "all"
    pseudobulks = keep %>%
      map( ~ {
        print(.)
        subgroup = .
        meta0 = meta
        expr0 = expr
        # catch cell types without replicates or conditions
        if (n_distinct(meta0$group) < 2)
          return(NA)
        replicate.counts = distinct(meta0, group, replicate) %>%
          group_by(group) %>%
          summarise(replicates = n_distinct(replicate)) %>%
          pull(replicates)
        if (any(replicate.counts < min.num.reps))
          return(NA)

        # process data into gene X replicate X subgroup matrice
        mm = model.matrix(~ 0 + replicate:group, data = meta0)

        mat.mm = expr0 %*% mm

        keep.genes = rowSums(as.array(mat.mm > 0)) > min.num.features

        mat.mm = mat.mm[keep.genes, ] %>% as.data.frame()
        mat.mm %<>% as.data.frame()
        colnames(mat.mm) = gsub("replicate|group", "", colnames(mat.mm))
        # drop empty columns
        keep.samples = colSums(mat.mm) > 0
        mat.mm %<>% magrittr::extract(, keep.samples)
        return(mat.mm)
      }) %>%
      setNames(keep)
  }
  
  # drop NAs
  pseudobulks %<>% magrittr::extract(!is.na(.))
  
  # also filter out cell types with no retained genes
  min.dim = map(pseudobulks, as.data.frame) %>% map(nrow)
  pseudobulks %<>% magrittr::extract(min.dim > 1)
  
  # also filter out types without replicates
  min.repl = map_int(pseudobulks, ~ {
    # make sure we have a data frame a not a vector
    tmp = as.data.frame(.)
    targets = data.frame(group.sample = colnames(tmp)) %>%
      mutate(group = gsub(".*\\:", "", group.sample))
    if (n_distinct(targets$group) == 1)
      return(as.integer(0))
    min(table(targets$group))
  })
  pseudobulks %<>% magrittr::extract(min.repl >= min.num.reps)

# Run a DE analysis within each cell type using Seurat
  # check args
  if (de.method == 'limma') {
    if (de.type != 'voom') {
      # change default type to use
      de.type = 'trend'  
    }
  }
  
  results = map(pseudobulks, function(x) {
    # create targets matrix
    targets = data.frame(group.sample = colnames(x)) %>%
      mutate(group = gsub(".*\\:", "", group.sample))
    ## optionally, carry over factor levels from entire dataset
    if (is.factor(meta$group)) {
      targets$group %<>% factor(levels = levels(meta$group))
    }
    # Need to set to higher values if more than 2 contrasts variables
    if (n_distinct(targets$group) > 100000)
      return(NULL)
    
    # create design
    dm.formula <- as.formula(paste("~0", paste(c("group")), sep="+", 
                                   collapse="+"))
    design = model.matrix(dm.formula, data = targets)
    
    DE = switch(de.method,
                edgeR = {
                  tryCatch({
                    y = DGEList(counts = x, group = targets$group) %>%
                      calcNormFactors(method = 'TMM') %>%
                      estimateDisp(design)
                    test = switch(de.type,
                                  QLF = {
                                    fit = glmQLFit(y, design)
                                    test = glmQLFTest(fit, coef = -1)
                                  },
                                  LRT = {
                                    fit = glmFit(y, design = design)
                                    test = glmLRT(fit)
                                  })
                    res = topTags(test, n = Inf) %>%
                      as.data.frame() %>%
                      rownames_to_column('gene') %>%
                      # flag metrics in results
                      mutate(de.family = 'pseudobulk',
                             de.method = de.method,
                             de.type = de.type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                DESeq2 = {
                  tryCatch({
                    dds = DESeqDataSetFromMatrix(countData = x,
                                                 colData = targets,
                                                 design = ~ group)
                    dds = switch(de.type,
                                 Wald = {
                                   dds = try(DESeq(dds,
                                                   test = 'Wald',
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 },
                                 LRT = {
                                   dds = try(DESeq(dds,
                                                   test = 'LRT',
                                                   reduced = ~ 1,
                                                   fitType = 'parametric',
                                                   sfType = 'poscounts',
                                                   betaPrior = F))
                                 }
                    )
                    res = results(dds)
                    # write
                    res = as.data.frame(res) %>%
                      mutate(gene = rownames(x)) %>%
                      # flag metrics in results
                      mutate(de.family = 'pseudobulk',
                             de.method = de.method,
                             de.type = de.type)
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                },
                limma = {
                  tryCatch({
                    x = switch(de.type,
                               trend = {
                                 trend_bool = T
                                 dge = DGEList(as.matrix(x), 
                                               group = targets$group)
                                 dge = calcNormFactors(dge)
                                 x = new("EList")
                                 x$E = cpm(dge, log = TRUE, prior.count = 3)
                                 x
                               },
                               voom = {
                                 counts = all(as.matrix(x) %% 1 == 0)
                                 if (counts) {
                                   trend_bool = F
                                   x = voom(as.matrix(x), design)
                                   x
                                 }
                               })
                    
                    # run contrasts                    
                    contrasts <- unlist(contrasts %>% strsplit("-") %>% 
                                          lapply(function (x) 
                                            paste("group",x,sep="")) %>% 
                                          lapply(paste, collapse = "-"))
                    
                    cm <- makeContrasts(contrasts = contrasts, levels=design)
                    
                    fit2 = lmFit(x, design) %>% contrasts.fit(cm) %>%
                      eBayes(trend = trend_bool, robust = trend_bool)
                    
                    logFC = fit2$coefficients
                    colnames(logFC)=paste(colnames(logFC),"logFC",sep="_")
                    tstat = fit2$t
                    colnames(tstat)=paste(colnames(tstat),"tstat",sep="_")
                    FC = 2^fit2$coefficients
                    FC = ifelse(FC<1,-1/FC,FC)
                    colnames(FC)=paste(colnames(FC),"FC",sep="_")
                    pvalall=fit2$p.value
                    colnames(pvalall)=paste(colnames(pvalall),"pval",sep="_")
                    pvaladjall=apply(pvalall,2,function(x) p.adjust(x,"BH"))
                    colnames(pvaladjall)=paste(colnames(fit2$coefficients),
                                               "adjpval",sep="_")
                    
                    res=as.data.frame(cbind(FC, logFC, tstat, pvalall, 
                                            pvaladjall))
                    
                    res$Gene <- rownames(res)
                    
                    res <- res %>% mutate(
                      de.family = 'pseudobulk',
                      de.method = de.method,
                      de.type = de.type)
                    
                  }, error = function(e) {
                    message(e)
                    data.frame()
                  })
                }
    )
  })
  results %<>% bind_rows(.id = 'subgroup')
  results <- results %>% dplyr::select("Gene", everything())

return(results)
}
