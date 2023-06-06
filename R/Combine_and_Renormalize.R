#' Combine_and_Renormalize template in single-cell-rna-seq-r4 NIDAP environment
#' from v178
#' 
#' 
#' @title Combine & Renormalize 
#' @description Combines samples, rescales and renormalizes, runs Dimensional Reduction, and returns a combined Seurat Object. This template will summarize the multi-dimensionality of your data into a set of "principal components" to allow for easier analysis. There is an option to use this template to perform Integration, as well. This is Step 4 in the canonical Single Cell pipeline after elbow plots and regression.   
#' @details
#' 
#' 
#' @param Seurat_Object description.
#' @param npcs Select the number of principal components for your analysis. Please see the elbow plot in the previous template to figure out what number of PCs explains your variance cut-off. For example, if the elbow plot has point at (15,0.02), it means that 15 PCs encapsulate 98% of the variance in your data.
#' @param vars_to_regress Subtract (‘regress out’) this source of heterogeneity from the data. For example, to regress out mitochondrial effects, input "percent.mt." Options: percent.mt, nCount_RNA, S.Score, G2M.Score, CC.Difference
#' @param integratedata Perform integration of cells across conditions using the most variant genes to identify cells most similar to each other.  N.B. Always look at cells before deciding whether to perform integration, and seek advice from bioinformatician.
#' @param clust_res_low Select minimum resolution for clustering plots. The lower you set this, the FEWER clusters will be generated.
#' @param clust_res_high Select the maximum resolution for clustering. The higher you set this number, the MORE clusters you will produced.
#' @param clust_res_bin Select the bins for your cluster plots. For example, if you input 0.2 as your bin, and have low/high resolution ranges of 0.2 and 0.6, then the template will produce cluster plots at resolutions of 0.2, 0.4 and 0.6.
#' @param only_var_genes If dataset is larger than ~40k filtered cells, toggle to TRUE. If TRUE, only variable genes will be available for downstream analysis. Default is FALSE.
#' @param draw_umap If TRUE, draw UMAP plot.
#' @param draw_tsne If TRUE, draw TSNE plot.
#' @param imageType Select output image type. Options: png, svg
#' @param nfeatures Number of variable features.
#' @param low_cut description.
#' @param high_cut description.
#' @param low_cut_disp description.
#' @param high_cut_disp description.
#' @param selection_method Method to choose top variable features. Options: vst, mean.var.plot, dispersion
#' @param cell_hashing_data Toggle "true" if you are using cell-hashed data.
#' @param project_name description.
#' @param doMergeData Toggle FALSE to stop metadata from merging. Not recommended for standard pipeline.
#' @param seed_for_PCA description.
#' @param seed_for_TSNE description.
#' @param seed_for_UMAP description.
#' @param Do_SCTransform Set to TRUE to run SCTransform (recommended current v3 Seurat default). Set to FALSE to run ScaleData (previous v2 Seurat default) instead.
#' @param Exclude_sample Exclude unwanted samples from the merge step. The number will correspond to the order in which they appear. Leave as 0 if you want to use all samples. If you want to exclude one or several samples, separate each sample number by comma (e.g. 1,2,3,4).
#' 
#' 
#' @import Seurat
#' @import ggplot2
#' @import gridExtra 
#' @import RColorBrewer
#' @import svglite
#' @importFrom svglite svglite
#' @importFrom digest digest
#' 
#' @export
#' 
#' @return Seurat Objects and QC plots. Combines samples, rescales and renormalizes, runs Dimensional Reduction, and returns a combined Seurat Object. This template will summarize the multi-dimensionality of your data into a set of "principal components" to allow for easier analysis. There is an option to use this template to perform Integration, as well. This is Step 4 in the canonical Single Cell pipeline after elbow plots and regression.


Combine_and_Renormalize <- function(Seurat_Object,
                                    npcs = 15,
                                    vars_to_regress = NULL,
                                    integratedata = FALSE,
                                    clust_res_low = 0.2,
                                    clust_res_high = 1.2,
                                    clust_res_bin = 0.2,
                                    only_var_genes = FALSE,
                                    draw_umap = TRUE,
                                    draw_tsne = TRUE,
                                    imageType = 'png',
                                    nfeatures = 2000,
                                    low_cut = 0.1,
                                    high_cut = 8,
                                    low_cut_disp = 1,
                                    high_cut_disp = 100000,
                                    selection_method = 'vst',
                                    cell_hashing_data = FALSE,
                                    project_name = 'scRNAProject',
                                    doMergeData = TRUE,
                                    seed_for_PCA = 42,
                                    seed_for_TSNE = 1,
                                    seed_for_UMAP = 42,
                                    Do_SCTransform = TRUE,
                                    Exclude_sample = 0
                          ){
  
  
  
  ##--------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # # load data
  # object.class <- getClass(class(SO))
  # 
  # if(object.class@className == "RFoundryObject") {
  #   cat("1. Reading Seurat Object from dataset: RObjectdata.rds\n\n")
  #   
  #   SO = SO$value
  # }else{
  #   cat("1. Reading Seurat Object from dataset: SO.rds\n\n")
  #   
  #   fs <- SO$fileSystem()
  #   path <- fs$get_path("seurat_object.rds", 'r')
  #   SO <- readRDS(path)
  #   
  # }
  
  SO <- Seurat_Object
  
  #in case you want to redo this on a merged SO
  if (class(SO) =="Seurat") {
    x =list()
    x[[1]] <- SO
    SO <- x
  }
  
  # If exclude option is TRUE, filter out undesirable sample
  if (c(Exclude_sample) == 0){
    SO <- SO
  } else {
    SO <- SO[-c(Exclude_sample)]
  }
  
  
  ###############################
  ## original settings v170
  conserve_memory=only_var_genes
  
  #Auto detect number of cells and turn on Conserve memory v173
  #    ## setting a limit for cell numbers
  #    if (dim(SO)[2]<35000) {
  #        too_many_cells <- FALSE
  #    }else{
  #        too_many_cells <- TRUE
  #    #cat(“There are too many cells, only 2000 variable genes would be #reported\n”)
  #    }
  #    
  #    if (too_many_cells || only_var_genes) {
  #         conserve_memory <- TRUE
  #    }else{
  #        conserve_memory <- FALSE
  #    }
  ################################
  
  
  
  ############################
  ## Detect Citeseq
  
  #initialize Citeseq functionality as false, 
  #later the template will check for a Protein assay and run if it finds it
  doCiteSeq <- FALSE
  
  dat = vector()
  
  if (length(SO) > 1) {
    for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
    SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = project_name, merge.data = TRUE)
    allgenes <- rownames(SO_merge)
    #SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
  } else {
    SO_merge <- SO[[1]]
    allgenes <- rownames(SO_merge)
    #SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
  }
  
  if (!("orig.ident" %in% colnames(SO_merge@meta.data))) {
    SO_merge@meta.data$orig.ident <- SO_merge@meta.data$orig_ident
  }
  
  if ("Protein" %in% names(SO_merge@assays)){
    doCiteSeq <-TRUE
  }
  
  if(cell_hashing_data){
    #SO_merge <- ScaleData(SO_merge, assay = "HTO")
  }
  
  ############################
  ## SCTransform
  
  if (Do_SCTransform){
    if(is.null(vars_to_regress)){
      SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, conserve.memory = conserve_memory, return.only.var.genes = FALSE)}
    else{       
      SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE,vars.to.regress=vars_to_regress, conserve.memory = conserve_memory, return.only.var.genes = FALSE) 
    }
  }
  else{
    all.genes <- rownames(SO_merge)
    if(is.null(vars_to_regress)){
      SO_merge <- SO_merge
    }
    else{
      #SO_merge <- ScaleData(SO_merge, features=all.genes, assay = "RNA", vars.to.regress=vars_to_regress) 
    }
    DefaultAssay(SO_merge) <- "RNA"   
  }
  
  
  ############################
  ## Integrate data
  
  if (length(SO)>1) {
    all_features <- lapply(SO, row.names) %>% Reduce(intersect, .)
    if(integratedata==TRUE){
      integ_features <- SelectIntegrationFeatures(object.list = SO, nfeatures = 3000) 
      if(!is.null(SO[[1]]@assays$SCT)){
        SO <- PrepSCTIntegration(object.list = SO, anchor.features = integ_features)
        k.filter <- min(200, min(sapply(SO, ncol)))
        integ_anchors <- FindIntegrationAnchors(object.list = SO, normalization.method = "SCT", k.filter=k.filter, anchor.features = integ_features)
        SO_merge <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT",features.to.integrate = all_features)
        #SO_merge <- ScaleData(SO_merge,features=all_features)
      }
      else{
        k.filter <- min(200, min(sapply(SO, ncol)))
        integ_anchors <- FindIntegrationAnchors(object.list = SO, k.filter=k.filter, anchor.features = integ_features)
        SO_merge <- IntegrateData(anchorset = integ_anchors,features.to.integrate = all_features)
        #SO_merge <- ScaleData(SO_merge,features=all_features)  
      }}
  }
  
  ############################
  ## Dimension reduction
  
  SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = nfeatures, mean.cutoff = c(low_cut, high_cut), dispersion.cutoff = c(low_cut_disp, high_cut_disp), selection.method = selection_method, verbose = FALSE)
  SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = seed_for_PCA)
  SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=seed_for_UMAP)
  SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = seed_for_TSNE)
  SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)
  
  
  ############################
  ## Citeseq
  
  #check for CITE-seq data and if so, run reductions
  if(doCiteSeq) {
    # Moved below integration step. SO_merge is recreated and this information was lost
    #SO_merge <- ScaleData(SO_merge, assay = "Protein")
    
    print("finding protein variable features...")
    VariableFeatures(SO_merge,assay="Protein") <- rownames(SO_merge$Protein)
    #Support for partial
    if(all(sapply(seq_along(SO),function(i) "Protein" %in% names(SO[[i]]@assays)))){
      print("running protein pca...")
      SO_merge <- RunPCA(object = SO_merge, assay="Protein",npcs = npcs,verbose = FALSE,reduction.name="protein_pca",seed.use = seed_for_PCA)
      SO_merge <- RunUMAP(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein), reduction.name="protein_umap",seed.use=seed_for_UMAP)
      SO_merge <- RunTSNE(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein),seed.use = seed_for_TSNE,reduction.name="protein_tsne",check_duplicates=F)
      SO_merge <- FindNeighbors(SO_merge, assay="Protein",graph.name="Protein_snn",features=rownames(SO_merge$Protein))
    }else{
      doCiteSeq <- FALSE #set to false so we don't cluster protein
    }
    
  } else {
    doCiteSeq <- FALSE
  }
  
  ############################
  ## Cluster
  
  for (i in seq(clust_res_low,clust_res_high,clust_res_bin)){
    SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
    if(doCiteSeq){
      SO_merge <- FindClusters(SO_merge, graph.name="Protein_snn",resolution = i, algorithm = 1)
    }
  }
  print("Clustering successful!")
  
  ############################
  ## create plots
  
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
  cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  grobsList = list()
  if(draw_tsne){
    p1 <- DimPlot(SO_merge, reduction = "tsne", group.by = "orig.ident", repel = TRUE,          pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) +         ggtitle("RNA TSNE")
    grobsList[[length(grobsList)+1]] <- p1
    print("Added RNA TSNE")
    print(length(grobsList))
    
  }
  if(draw_umap){
    p2 <- DimPlot(SO_merge, reduction = "umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("RNA UMAP")
    grobsList[[length(grobsList)+1]] <- p2
    print("Added RNA UMAP")
    print(length(grobsList))
    
  }
  if(draw_tsne & doCiteSeq){ 
    p3 <- DimPlot(SO_merge, reduction = "protein_tsne", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody TSNE")
    grobsList[[length(grobsList)+1]] <- p3
    print("Added Antibody TSNE")
    print(length(grobsList))
    
    
  }
  if(draw_umap & doCiteSeq){ 
    p4 <- DimPlot(SO_merge, reduction = "protein_umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody UMAP")
    grobsList[[length(grobsList)+1]] <- p4
    print("Added Antibody UMAP")
    print(length(grobsList))
    
  }
  
  ############################
  ## create Figure output
  
  n = ceiling(length(grobsList)^0.5)
  m=ceiling(length(grobsList)/n)
  imageWidth = 1200*n
  imageHeight = 1200*m
  dpi = 300
  
  grobs=arrangeGrob(grobs=grobsList,ncol=n)
  
  if (imageType == 'png') {
    png(
      # filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=2,
      bg="white",
      res=dpi,
      type="cairo")
  } else {
    svglite::svglite(
      # file=graphicsFile,
      width=round(imageWidth/dpi,digits=2),
      height=round(imageHeight/dpi,digits=2),
      pointsize=1,
      bg="white")
  }

  
  # plot(grobs)
  
  grobs=grobs
  
  #slot(SO_merge,"commands") <- list()
  cat("\nPCA Object Checksum:\n")
  
  return(list(so=SO_merge,plot=grobs))
  
  # output <- new.output()
  # output_fs <- output$fileSystem()
  # saveRDS(SO_merge, output_fs$get_path("seurat_object.rds", 'w'))
  # 
  # return(output_fs)

}