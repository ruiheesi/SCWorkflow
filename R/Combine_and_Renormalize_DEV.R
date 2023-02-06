#' Combine_and_Renormalize template in single-cell-rna-seq-r4 NIDAP environment
#' from v178
#' 
#' 
#' @title Combine & Renormalize 
#' @description Combines samples, rescales and renormalizes, runs Dimensional Reduction, and returns a combined Seurat Object. This template will summarize the multi-dimensionality of your data into a set of "principal components" to allow for easier analysis. There is an option to use this template to perform Integration, as well. This is Step 4 in the canonical Single Cell pipeline after elbow plots and regression.   
#' @details
#' 
#' 
#' @param object description.
#' @param npcs Select the number of principal components for your analysis. Please see the elbow plot in the previous template to figure out what number of PCs explains your variance cut-off. For example, if the elbow plot has point at (15,0.02), it means that 15 PCs encapsulate 98% of the variance in your data.
#' @param vars_to_regress Subtract (‘regress out’) this source of heterogeneity from the data. For example, to regress out mitochondrial effects, input "percent.mt." Options: percent.mt, nCount_RNA, S.Score, G2M.Score, CC.Difference
#' @param integratedata Perform integration of cells across conditions using the most variant genes to identify cells most similar to each other.  N.B. Always look at cells before deciding whether to perform integration, and seek advice from bioinformatician.
#' @param clust.res.low Select minimum resolution for clustering plots. The lower you set this, the FEWER clusters will be generated.
#' @param clust.res.high Select the maximum resolution for clustering. The higher you set this number, the MORE clusters you will produced.
#' @param clust.res.bin Select the bins for your cluster plots. For example, if you input 0.2 as your bin, and have low/high resolution ranges of 0.2 and 0.6, then the template will produce cluster plots at resolutions of 0.2, 0.4 and 0.6.
#' @param only.var.genes If dataset is larger than ~40k filtered cells, toggle to TRUE. If TRUE, only variable genes will be available for downstream analysis. Default is FALSE.
#' @param draw.umap If TRUE, draw UMAP plot.
#' @param draw.tsne If TRUE, draw TSNE plot.
#' @param imageType Select output image type. Options: png, svg
#' @param nfeatures Number of variable features.
#' @param low.cut description.
#' @param high.cut description.
#' @param low.cut.disp description.
#' @param high.cut.disp description.
#' @param selection.method Method to choose top variable features. Options: vst, mean.var.plot, dispersion
#' @param cell.hashing.data Toggle "true" if you are using cell-hashed data.
#' @param project.name description.
# @param doMergeData Toggle FALSE to stop metadata from merging. Not recommended for standard pipeline. # doMergeData doesn't do anything last used in vesion 34 so removed as parameter
#' @param seed.for.PCA description.
#' @param seed.for.TSNE description.
#' @param seed.for.UMAP description.
#' @param SCTransform Set to TRUE to run SCTransform (recommended current v3 Seurat default). Set to FALSE to run ScaleData (previous v2 Seurat default) instead.
#' @param exclude.sample Exclude unwanted samples from the merge step. Include sample names to be removed. If you want to exclude several samples, separate each sample number by comma (e.g. sample1,sample2,sample3,sample4).
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


Combine_and_Renormalize_DEV <- function(Seurat.Object,
                                    npcs = 15,
                                    vars.to.regress = NULL,
                                    integratedata = FALSE,
                                    clust.res.low = 0.2,
                                    clust.res.high = 1.2,
                                    clust.res.bin = 0.2,
                                    only.var.genes = FALSE,
                                    draw.umap = TRUE,
                                    draw.tsne = TRUE,
                                    imageType = 'png',
                                    nfeatures = 2000,
                                    low.cut = 0.1,
                                    high.cut = 8,
                                    low.cut.disp = 1,
                                    high.cut.disp = 100000,
                                    selection.method = 'vst',
                                    cell.hashing.data = FALSE,
                                    project.name = 'scRNAProject',
                                    # doMergeData = TRUE,
                                    seed.for.PCA = 42,
                                    seed.for.TSNE = 1,
                                    seed.for.UMAP = 42,
                                    SCTransform = TRUE,
                                    exclude.sample = ""
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
  #   path <- fs$get.path("seurat.object.rds", 'r')
  #   SO <- readRDS(path)
  #   
  # }
  
  SO <- Seurat.Object
  
  #in case you want to redo this on a merged SO
  if (class(SO) =="Seurat") {
    x =list()
    x[[1]] <- SO
    SO <- x
  }
  
  # If exclude option is TRUE, filter out undesirable sample
  if (length(c(exclude.sample)) == 0){
    SO <- SO
  } else {
    # SO <- SO[-c(exclude.sample)]
    SO <- SO[!names(SO)%in%exclude.sample]
  }
  
  
  ###############################
  ## original settings v170
  # conserve.memory=only.var.genes
  
  ## Calculate total number of cells in input SO.
  cell.count <- sum(unlist((lapply(SO, function(x) dim(x)[2]))))
  
  ## Auto detect number of cells and turn on Conserve memory v173
  ## Setting a limit for cell numbers
  # if (dim(SO)[2] < 35000) {
  if (cell.count < 35000) {
    too.many.cells <- FALSE
  } else {
    too.many.cells <- TRUE
    # cat(“There are too many cells, only 2000 variable genes would be #reported\n”)
  }
  
  ## 
  if (too.many.cells || only.var.genes) {
    conserve.memory <- TRUE
  } else {
    conserve.memory <- FALSE
  }
  
  ################################
  if(SCT==1){
    
    if(is.null(vars.to.regress)){
      SOx <- lapply(SO,function(x){SCTransform(x,do.correct.umi = TRUE, conserve.memory = conserve.memory, return.only.var.genes = FALSE)})
    }else{       
      SO <- lapply(SO,function(x){SCTransform(x,do.correct.umi = TRUE,vars.to.regress=vars.to.regress, conserve.memory = conserve.memory, return.only.var.genes = FALSE)}) 
    }
    

  }
  ############################
  ## Detect Citeseq
  
  #initialize Citeseq functionality as false, 
  #later the template will check for a Protein assay and run if it finds it
  doCiteSeq <- FALSE
  
  dat = vector()
  
  if (length(SO) > 1) {
    for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
    SO.merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = project.name, merge.data = TRUE)
    allgenes <- rownames(SO.merge)
    #SO.merge <- ScaleData(SO.merge, assay = "RNA", features=allgenes)
  } else {
    SO.merge <- SO[[1]]
    allgenes <- rownames(SO.merge)
    #SO.merge <- ScaleData(SO.merge, assay = "RNA", features=allgenes)
  }
  
  if (!("orig.ident" %in% colnames(SO.merge@meta.data))) {
    SO.merge@meta.data$orig.ident <- SO.merge@meta.data$orig.ident
  }
  
  if ("Protein" %in% names(SO.merge@assays)){
    doCiteSeq <-TRUE
  }
  
  if(cell.hashing.data){
    SO.merge <- ScaleData(SO.merge, assay = "HTO")
  }
  
  ############################
  ## SCTransform
  
  
  if (SCTransform){
    if(is.null(vars.to.regress)){
      SO.merge <- SCTransform(SO.merge,do.correct.umi = TRUE, conserve.memory = conserve.memory, return.only.var.genes = FALSE)
      }else{       
      SO.merge <- SCTransform(SO.merge,do.correct.umi = TRUE,vars.to.regress=vars.to.regress, conserve.memory = conserve.memory, return.only.var.genes = FALSE) 
      }
  }else{
    all.genes <- rownames(SO.merge)
    if(is.null(vars.to.regress)){
      SO.merge <- SO.merge
    }else{
      #SO.merge <- ScaleData(SO.merge, features=all.genes, assay = "RNA", vars.to.regress=vars.to.regress) 
    }
    DefaultAssay(SO.merge) <- "RNA"   
  }
  
  
  ############################
  ## Integrate data
  
  if (length(SO)>1) {
    all.features <- lapply(SO, row.names) %>% Reduce(intersect, .)
    if(integratedata==TRUE){
      integ.features <- SelectIntegrationFeatures(object.list = SO, nfeatures = 3000) 
      if(!is.null(SO[[1]]@assays$SCT)){
        SO <- PrepSCTIntegration(object.list = SO, anchor.features = integ.features)
        k.filter <- min(200, min(sapply(SO, ncol)))
        integ.anchors <- FindIntegrationAnchors(object.list = SO, normalization.method = "SCT", k.filter=k.filter, anchor.features = integ.features)
        SO.merge <- IntegrateData(anchorset = integ.anchors, normalization.method = "SCT",features.to.integrate = all.features)
        #SO.merge <- ScaleData(SO.merge,features=all.features)
      }else{
        k.filter <- min(200, min(sapply(SO, ncol)))
        integ.anchors <- FindIntegrationAnchors(object.list = SO, k.filter=k.filter, anchor.features = integ.features)
        SO.merge <- IntegrateData(anchorset = integ.anchors,features.to.integrate = all.features)
        #SO.merge <- ScaleData(SO.merge,features=all.features)  
      }
    }
  }
  
  ############################
  ## Dimension reduction
  
  SO.merge <- FindVariableFeatures(object = SO.merge, nfeatures = nfeatures, mean.cutoff = c(low.cut, high.cut), dispersion.cutoff = c(low.cut.disp, high.cut.disp), selection.method = selection.method, verbose = FALSE)
  SO.merge <- RunPCA(object = SO.merge, npcs = npcs, verbose = FALSE,seed.use = seed.for.PCA)
  SO.merge <- RunUMAP(object = SO.merge, reduction = "pca", dims = 1:npcs, seed.use=seed.for.UMAP)
  SO.merge <- RunTSNE(object = SO.merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = seed.for.TSNE)
  SO.merge <- FindNeighbors(SO.merge, dims = 1:npcs)
  
  
  ############################
  ## Citeseq
  
  #check for CITE-seq data and if so, run reductions
  if(doCiteSeq) {
    # Moved below integration step. SO.merge is recreated and this information was lost
    #SO.merge <- ScaleData(SO.merge, assay = "Protein")
    
    print("finding protein variable features...")
    VariableFeatures(SO.merge,assay="Protein") <- rownames(SO.merge$Protein)
    #Support for partial
    if(all(sapply(seq.along(SO),function(i) "Protein" %in% names(SO[[i]]@assays)))){
      print("running protein pca...")
      SO.merge = ScaleData(SO.merge, assay = 'Protein', verbose = FALSE) ####Add this line here ################
      SO.merge <- RunPCA(object = SO.merge, assay="Protein",npcs = npcs,verbose = FALSE,reduction.name="protein.pca",seed.use = seed.for.PCA)
      SO.merge <- RunUMAP(object = SO.merge, assay="Protein", features=rownames(SO.merge$Protein), reduction.name="protein.umap",seed.use=seed.for.UMAP)
      SO.merge <- RunTSNE(object = SO.merge, assay="Protein", features=rownames(SO.merge$Protein),seed.use = seed.for.TSNE,reduction.name="protein.tsne",check.duplicates=F)
      SO.merge <- FindNeighbors(SO.merge, assay="Protein",graph.name="Protein.snn",features=rownames(SO.merge$Protein))
    }else{
      doCiteSeq <- FALSE #set to false so we don't cluster protein
    }
    
  } else {
    doCiteSeq <- FALSE
  }
  
  ############################
  ## Cluster
  
  for (i in seq(clust.res.low,clust.res.high,clust.res.bin)){
    SO.merge <- FindClusters(SO.merge, resolution = i, algorithm = 1)
    if(doCiteSeq){
      SO.merge <- FindClusters(SO.merge, graph.name="Protein.snn",resolution = i, algorithm = 1)
    }
  }
  print("Clustering successful!")
  
  ############################
  ## create plots
  
  n <- 60
  qual.col.pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual.col.pals = qual.col.pals[c(7,6,2,1,8,3,4,5),]
  cols = unlist(mapply(brewer.pal, qual.col.pals$maxcolors, rownames(qual.col.pals)))
  
  grobsList = list()
  if(draw.tsne){
    p1 <- DimPlot(SO.merge, reduction = "tsne", group.by = "orig.ident", repel = TRUE,          pt.size=0.02) + theme.classic() + scale.color.manual(values=cols) + theme(legend.position="top", legend.text=element.text(size=5)) +
      guides(colour = guide.legend(ncol=3, override.aes = list(size=1, alpha = 1))) +         ggtitle("RNA TSNE")
    grobsList[[length(grobsList)+1]] <- p1
    print("Added RNA TSNE")
    print(length(grobsList))
    
  }
  if(draw.umap){
    p2 <- DimPlot(SO.merge, reduction = "umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme.classic() + scale.color.manual(values=cols) + theme(legend.position="top", legend.text=element.text(size=5)) +
      guides(colour = guide.legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("RNA UMAP")
    grobsList[[length(grobsList)+1]] <- p2
    print("Added RNA UMAP")
    print(length(grobsList))
    
  }
  if(draw.tsne & doCiteSeq){ 
    p3 <- DimPlot(SO.merge, reduction = "protein.tsne", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme.classic() + scale.color.manual(values=cols) + theme(legend.position="top", legend.text=element.text(size=5)) +
      guides(colour = guide.legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody TSNE")
    grobsList[[length(grobsList)+1]] <- p3
    print("Added Antibody TSNE")
    print(length(grobsList))
    
    
  }
  if(draw.umap & doCiteSeq){ 
    p4 <- DimPlot(SO.merge, reduction = "protein.umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme.classic() + scale.color.manual(values=cols) + theme(legend.position="top", legend.text=element.text(size=5)) +
      guides(colour = guide.legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody UMAP")
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
  
  # if (imageType == 'png') {
  #   png(
  #     # filename=graphicsFile,
  #     width=imageWidth,
  #     height=imageHeight,
  #     units="px",
  #     pointsize=2,
  #     bg="white",
  #     res=dpi,
  #     type="cairo")
  # } else {
  #   svglite::svglite(
  #     # file=graphicsFile,
  #     width=round(imageWidth/dpi,digits=2),
  #     height=round(imageHeight/dpi,digits=2),
  #     pointsize=1,
  #     bg="white")
  # }
  
  # plot(grobs)
  
  grobs=grobs
  
  #slot(SO.merge,"commands") <- list()
  # cat("\nPCA Object Checksum:\n")
  
  return(list(so=SO.merge,plot=grobs))
  
  # output <- new.output()
  # output.fs <- output$fileSystem()
  # saveRDS(SO.merge, output.fs$get.path("seurat.object.rds", 'w'))
  # 
  # return(output.fs)
}