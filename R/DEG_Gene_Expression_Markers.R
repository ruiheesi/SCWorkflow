# DEG (Gene Expression Markers) [scRNA-seq] [CCBR] (54b6dd44-e233-4fb8-9bae-6ab0cb46e399): v81

#' @title DEG (Gene Expression Markers)
#' @description This function performs a DEG (differential expression of genes) analysis on a merged Seurat object to identify expression markers between different groups of cells (contrasts).
#' @details The recommended input is a merged Seurat object with SingleR annotations, along with its associated sample names and metadata.
#' 
#' @param object Seurat-class object.
#' @param samples Samples to be included in the analysis. 
#' @param parameter.to.test Select the metadata column that you would like to use to perform your DEG analysis and construct your contrasts from.
#' @param contrasts Contrasts in the "A-B" format.
#' @param test.to.use The kind of algorithm you would like to use to perform your DEG analysis. Default is the MAST algorithm (wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2).
#' @param log.fc.threshold The minimum log fold-change between contrasts that you would like to analyze.
#' @param use.spark Opt to use Spark to parallelize computations. 
#' @param assay.to.use The assay to use for your DEG analysis. Default is SCT, but can use linearly scaled data by selecting RNA instead.
#' @param use.log.2 Set to FALSE if you want to use default Seurat3 settings for Fold Change (natural logarithm)
#' @param latent.vars Select possible confounding variable within metadata to account for when running DEG (Works only when test.use is one of 'LR', 'negbinom', 'poisson', or 'MAST'). Leave blank if you don't have any confounding variables you want to take into account.

#' @import Seurat
#' @import ggplot2 
#' @import RColorBrewer 
#' @import scales 
#' @import tidyverse 
#' @import ggrepel 
#' @import gdata 
#' @import reshape2 
#' @import tools
#' @import grid
#' @import gridBase
#' @import gridExtra
#' @import parallel
#' @import MAST
#' 
#'   
#' @export 
#' 
#' @return a dataframe with DEG.



degGeneExpressionMarkers <- function(object,
  samples,
  parameter.to.test = "orig_ident",
  contrasts,
  test.to.use = "MAST",
#  test.to.use = "negbinom",
  log.fc.threshold = 0.25,
  use.spark = FALSE,
  assay.to.use = "SCT",
#use.log.2 = FALSE,
  latent.vars = c()
  ){

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  metadata_table <- object@meta.data
  samples = eval(parse(text=gsub('\\[\\]','c()',samples)))
  
  if (length(samples) == 0) {
    samples = unique(object@meta.data$sample_name)
  }
  
  colnames(object@meta.data) <- gsub("orig_ident","orig.ident",colnames(object@meta.data))
  if("active.ident" %in% slotNames(object)){
    sample_name = as.factor(object@meta.data$orig.ident)
    names(sample_name)=names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample_name
    object.sub = subset(object, ident = samples)
  } else {
    sample_name = as.factor(object@meta.data$orig.ident)
    names(sample_name)=names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample_name
    object.sub = subset(object, ident = samples)
  }
  
  print("selected samples:")
  print(object.sub)
  
  colnames(object.sub@meta.data) = gsub("\\.","_",colnames(object.sub@meta.data))
  
  #define contrasts
  newcont <- list()
  for (i in 1:length(contrasts)){
    newcont[[i]] <- c(paste(unlist(strsplit(contrasts[i],"-"))))
  }
  contrasts <- newcont
  
  #ERROR CATCHING
  #collect valid names of valid columns
  validColumns <- character()
  for (i in colnames(metadata_table)) {
    if (!any(is.na(metadata_table[[i]]))) {
      validColumns <-c(validColumns,i)
    }
  }
  
  param2test <- parameter.to.test
  
  if (param2test =="") {
    mcols = colnames(object.sub@meta.data)
    param2test <-mcols[grepl("RNA_snn",mcols)][[1]]
    print(paste("No parameter selected, defaulting to",param2test))
  }
  
  contrastTarget <- object.sub@meta.data[[param2test]]
  contrastType <- param2test
  contrastCounts = as.data.frame(table(contrastTarget))
  validContrasts = subset(contrastCounts, Freq>2)[[1]]
  
  #catch malformed contrasts
  for (i in contrasts) {
    if (!(i[[1]] %in% contrastTarget)) {
      print(paste(i[[1]],"is not a valid contrast for contrast type:", contrastType))
      print("Please see below for an example of valid contrasts for your selected contrast type.")
      print(validContrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[2]] %in% contrastTarget) & (i[[2]] != "all")) {
      print(paste(i[[2]],"is not a valid contrast for contrast type:", contrastType))
      print("Please see below for an example of valid contrasts for your selected contrast type.")
      print(validContrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (length(i)>2) {
      print("Contrasts are as follows..")
      print(i)
      stop("The console says there are too many inputs in your contrasts. A contrast should only contain Group1-Group2, but the console thinks you have inputed Group1-Group2-Group3")
    } else if (!(i[[2]] %in% validContrasts) & (i[[2]] != "all")) {
      print(paste(i[[2]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[1]] %in% validContrasts)) {
      print(paste(i[[1]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
      stop("You have entered an invalid group to contrast against.")
    }
  }
  
  #print out contrast cell contrastCounts
  for (i in seq_along(contrasts)) {
    firstGroup <- contrasts[[i]][[1]]
    firstGroupCount <- subset(contrastCounts, contrastTarget == firstGroup)$Freq
    if  (contrasts[[i]][[2]]!= "all") {
      secondGroup <-contrasts[[i]][[2]]
      secondGroupCount <-subset(contrastCounts, contrastTarget == secondGroup)$Freq      
      print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. cluster",secondGroup,"with",secondGroupCount,"cells."))
    } else {
      secondGroupCount <-ncol(object.sub)-firstGroupCount
      print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. all other clusters, totalling",secondGroupCount,"cells."))
    } 
  }
  
  #define and call function for running DEG
  get_deg_table <- function(n) {
##    library(Seurat)
    
    firstCluster <-unlist(n)[1]
    secondCluster <- unlist(n)[2]
    
    

    if (secondCluster == "all") { secondCluster <- NULL }
    
    Idents(object.sub) <- param2test
    
    ### REMOVED for Seurat4!!!
    #workaround for log2/ln changes:
##    if (use.log.2) { log.fc.threshold <- log.fc.threshold/log2(exp(1)) }
    
###    markers = FindMarkers(object.sub, ident.1 = firstCluster, ident.2 = secondCluster, test.use = test.to.use, logfc.threshold = log.fc.threshold, verbose=FALSE, assay = assay.to.use, latent.vars = eval(parse(text = "latent.vars")))
    markers = FindMarkers(object.sub, ident.1 = firstCluster, ident.2 = secondCluster, test.use = test.to.use, logfc.threshold = log.fc.threshold, verbose=FALSE, assay = assay.to.use, slot="counts")

    
    ## It works, but not with MAST:
#    markers = FindMarkers(object.sub,
#                          ident.1 = firstCluster,
#                          ident.2 = secondCluster,
#                          test.use = "negbinom",
#                          logfc.threshold = log.fc.threshold,
#                          verbose=FALSE,
#                          assay = assay.to.use,
#                          latent.vars = eval(parse(text = "latent.vars")))
    
    colnames(markers) <- chartr(old=" ",new="_",paste(colnames(markers), firstCluster,"vs",secondCluster,sep = "_"))
    
    ### REMOVED for Seurat4!!!
##        if (use.log.2) { markers[grep("avg_logFC_", colnames(markers))] <- markers[grep("avg_logFC_", colnames(markers))]*log2(exp(1)) }
    
    return(markers)
  }
  
  if (use.spark) {
    deg_tables <- spark.lapply(contrasts, get_deg_table) 
  } else {
    deg_tables <- lapply(contrasts, get_deg_table) 
  }
  
  for(i in seq_along(deg_tables)){
    degtab <- deg_tables[[i]]
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] > 0) %>% dim() -> pos 
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] < 0) %>% dim() -> neg
    print(paste0("The number of upregulated genes at p<0.05 in contrast number ", i, " is:"))
    print(pos[1])
    print(paste0("The number of downregulated genes at p<0.05 in contrast number ", i, " is:"))
    print(neg[1]) 
  }
  
  #Merge the deg tables together
  out_df <- NULL
  for (i in deg_tables) {
    if (is.null(out_df)) {
      out_df <- deg_tables[1]
      out_df <- as.data.frame(out_df)
    } else {
      out_df <- merge(out_df, i, by="row.names", all=TRUE)
      rownames(out_df) <- out_df$Row.names #set the rownames
      out_df$Row.names <- NULL #drop the row.names columns which we no longer need
    }
  }
  
  out_df$Gene <- rownames(out_df)
  out_df$Row.names <- NULL
  out_df <- out_df %>% dplyr::select(Gene, everything())
#  return(out_df)
  
  result.list <- list("out_df" = out_df)
  return(result.list)
  
}
