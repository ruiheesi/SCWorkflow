#' 
#' @title Combine & Renormalize 
#' @description Combines samples, rescales and renormalizes, 
#' runs Dimensional Reduction, and returns a combined Seurat Object. 
#' @details This is Step 4 in the basic Single-Cell RNA-seq workflow.
#'  This template will summarize the multi-dimensionality of your data into 
#'  a set of "principal components" to allow for easier analysis. 
#'  There is an option to use this template to perform Integration, as well.
#' 
#' 
#' @param object description.
#' @param npcs Select the number of principal components for your analysis.
#'  Please see the elbow plot in the previous template to figure out what
#'  number of PCs explains your variance cut-off. For example,
#'  if the elbow plot has point at (15,0.02), it means that 15 PCs
#'  encapsulate 98% of the variance in your data.(Default: 15)
#' @param vars.to.regress Subtract (‘regress out’) this source of heterogeneity
#'  from the data. For example, to Subtract mitochondrial effects, 
#'  input "percent.mt." Options: percent.mt, nCount.RNA, 
#'  S.Score, G2M.Score, CC.Difference
#' @param integrate.data Perform integration of cells across conditions using
#'  the most variant genes to identify cells most similar to each other.
#'  N.B. Always look at cells before deciding whether to perform integration,
#'  and seek advice from bioinformatician.
#' @param clust.res.low Select minimum resolution for clustering plots.
#'  The lower you set this, the FEWER clusters will be generated.
#' @param clust.res.high Select the maximum resolution for clustering.
#'  The higher you set this number, the MORE clusters you will produced.
#' @param clust.res.bin Select the bins for your cluster plots.
#'  For example, if you input 0.2 as your bin, and have low/high resolution 
#'  ranges of 0.2 and 0.6, then the template will produce cluster plots 
#'  at resolutions of 0.2, 0.4 and 0.6.
#' @param only.var.genes If dataset is larger than ~40k filtered cells,
#'  toggle to TRUE. If TRUE, only variable genes will be available for
#'  downstream analysis. Default is FALSE.
#' @param draw.umap If TRUE, draw UMAP plot.
#' @param draw.tsne If TRUE, draw TSNE plot.
#' @param nfeatures Number of variable features.
#' @param low.cut description.
#' @param high.cut description.
#' @param low.cut.disp description.
#' @param high.cut.disp description.
#' @param selection.method Method to choose top variable features.
#'  Options: vst, mean.var.plot, dispersion
#' @param cell.hashing.data Toggle "true" if you are using cell-hashed data
#' @param project.name description.
#' @param seed.for.pca description.
#' @param seed.for.tsne description.
#' @param seed.for.umap description.
#' @param sctransform Set to TRUE to run SCTransform
#'  (recommended current v3 Seurat default). Set to FALSE to run ScaleData
#'  (previous v2 Seurat default) instead.
#' @param exclude.sample Exclude unwanted samples from the merge step.
#'  Include sample names to be removed. If you want to exclude several samples,
#'  separate each sample number by comma (e.g. sample1,sample2,sample3,sample4).
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
#' @return Seurat Objects and QC plots


combineRenormalize <- function(object,
                               npcs = 15,
                               vars.to.regress = NULL,
                               integrate.data = FALSE,
                               clust.res.low = 0.2,
                               clust.res.high = 1.2,
                               clust.res.bin = 0.2,
                               only.var.genes = FALSE,
                               draw.umap = TRUE,
                               draw.tsne = TRUE,
                               nfeatures = 2000,
                               low.cut = 0.1,
                               high.cut = 8,
                               low.cut.disp = 1,
                               high.cut.disp = 100000,
                               selection.method = 'vst',
                               cell.hashing.data = FALSE,
                               project.name = 'scRNAProject',
                               seed.for.pca = 42,
                               seed.for.tsne = 1,
                               seed.for.umap = 42,
                               sctransform = TRUE,
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
  
  
  #in case you want to redo this on a merged SO
  if (class(object) =="Seurat") {
    x =list()
    x[[1]] <- object
    object <- x
  }
  
  # If exclude option is TRUE, filter out undesirable sample
  if (length(c(exclude.sample)) == 0){
    object <- object
  } else {
    # object <- object[-c(exclude.sample)]
    object <- object[!names(object)%in%exclude.sample]
  }
  
  
  ## Auto detect number of cells and turn on Conserve memory ####
  
  ## Calculate total number of cells in input SO.
  cell.count <- sum(unlist((lapply(object, function(x) dim(x)[2]))))
  
  ## Setting a limit for cell numbers
  # if (dim(object)[2] < 35000) {
  if (cell.count < 35000) {
    too.many.cells <- FALSE
  } else {
    too.many.cells <- TRUE
  }
  
  ## 
  if (too.many.cells || only.var.genes) {
    conserve.memory <- TRUE
  } else {
    conserve.memory <- FALSE
  }
  
  
  
  ## Merge samples into single SO ####
  
  dat = vector()
  
  if (length(object) > 1) {
    for(i in 2:length(object)){dat=c(dat,object[[i]]) }
    object.merge <- merge(object[[1]], y = dat, add.cell.ids = names(object), project = project.name, merge.data = TRUE)
    allgenes <- rownames(object.merge)
  } else {
    object.merge <- object[[1]]
    allgenes <- rownames(object.merge)
  }
  
  if (!("orig.ident" %in% colnames(object.merge@meta.data))) {
    object.merge@meta.data$orig.ident <- object.merge@meta.data$orig.ident
  }
  
  
  ## Detect Citeseq ####
  
  #initialize Citeseq functionality as false, 
  #later the template will check for a Protein assay and run if it finds it
  
  do.cite.seq <- FALSE
  
  if ("Protein" %in% names(object.merge@assays)){
    do.cite.seq <-TRUE
  }
  
  
  ## HTO ####
  if(cell.hashing.data){
    object.merge <- ScaleData(object.merge, assay = "HTO")
  }
  
  
  ## SCTransform on merge SO ####
  
  if (sctransform==T){
    if(is.null(vars.to.regress)){ 
      ## No Regression Variables
      object.merge <- SCTransform(object.merge,
                                  do.correct.umi = TRUE, 
                                  conserve.memory = conserve.memory, 
                                  return.only.var.genes = FALSE)
    }else{
      ## With Regression Variables  
      object.merge <- SCTransform(object.merge,
                                  do.correct.umi = TRUE,
                                  vars.to.regress=vars.to.regress, 
                                  conserve.memory = conserve.memory, 
                                  return.only.var.genes = FALSE) 
    }
  }else{
    all.genes <- rownames(object.merge)
    if(is.null(vars.to.regress)){
      object.merge <- object.merge
    }else{
      print("SCTransform not Preformed")
    }
    DefaultAssay(object.merge) <- "RNA"   
  }
  
  
  
  ## Integrate data ####
  
  if (length(object)>1) {
    all.features <- lapply(object, row.names) %>% Reduce(intersect, .)
    if(integrate.data==TRUE){
      integ.features <- SelectIntegrationFeatures(object.list = object, 
                                                  nfeatures = 3000) 
      
      if(!is.null(object[[1]]@assays$SCT)){
        ## Integrate SCTransform Data
        object <- PrepSCTIntegration(object.list = object, 
                                     anchor.features = integ.features)
        k.filter <- min(200, min(sapply(object, ncol)))
        integ.anchors <- FindIntegrationAnchors(object.list = object, 
                                                normalization.method = "SCT", 
                                                k.filter=k.filter, 
                                                anchor.features=integ.features)
        object.merge <- IntegrateData(anchorset = integ.anchors, 
                                      normalization.method = "SCT",
                                      features.to.integrate = all.features)
      }else{
        ## Integrate Counts Data
        k.filter <- min(200, min(sapply(object, ncol)))
        integ.anchors <- FindIntegrationAnchors(object.list = object, 
                                                k.filter=k.filter, 
                                                anchor.features = integ.features)
        object.merge <- IntegrateData(anchorset = integ.anchors,
                                      features.to.integrate = all.features)
      }
    }
  }
  
  
  ## Dimension reduction ####
  
  object.merge <- FindVariableFeatures(object = object.merge, 
                                       nfeatures = nfeatures, 
                                       mean.cutoff = c(low.cut, high.cut), 
                                       dispersion.cutoff=c(low.cut.disp,high.cut.disp), 
                                       selection.method = selection.method, 
                                       verbose = FALSE)
  object.merge <- RunPCA(object = object.merge, 
                         npcs = npcs, verbose = FALSE,
                         seed.use = seed.for.pca)
  object.merge <- RunUMAP(object = object.merge, 
                          reduction = "pca", 
                          dims = 1:npcs, 
                          seed.use=seed.for.umap)
  object.merge <- RunTSNE(object = object.merge, 
                          reduction = "pca", 
                          dim.embed = 2, 
                          dims = 1:npcs, 
                          seed.use = seed.for.tsne)
  object.merge <- FindNeighbors(object.merge, dims = 1:npcs)
  
  
  ## Citeseq ####
  
  #check for CITE-seq data and if so, run reductions
  if(do.cite.seq) {
    object.merge <- ScaleData(object.merge, assay = "Protein")
    
    print("finding protein variable features...")
    VariableFeatures(object.merge,assay="Protein") <- rownames(object.merge$Protein)
    
    #Support for partial
    if(all(sapply(seq.along(object),function(i) "Protein" %in% names(object[[i]]@assays)))){
      print("running protein pca...")
      object.merge = ScaleData(object.merge, 
                               assay = 'Protein',
                               verbose = FALSE) 
      object.merge <- RunPCA(object = object.merge, 
                             assay="Protein",
                             npcs = npcs,
                             reduction.name="protein.pca",
                             seed.use = seed.for.pca,
                             verbose = FALSE)
      object.merge <- RunUMAP(object = object.merge, 
                              assay="Protein", 
                              features=rownames(object.merge$Protein), 
                              reduction.name="protein.umap",
                              seed.use=seed.for.umap)
      object.merge <- RunTSNE(object = object.merge, 
                              assay="Protein", 
                              features=rownames(object.merge$Protein),
                              seed.use = seed.for.tsne,
                              reduction.name="protein.tsne",
                              check.duplicates=F)
      object.merge <- FindNeighbors(object.merge, assay="Protein",
                                    graph.name="Protein.snn",
                                    features=rownames(object.merge$Protein))
    }else{
      do.cite.seq <- FALSE #set to false so we don't cluster protein
    }
    
  } else {
    do.cite.seq <- FALSE
  }
  
  
  ## Cluster ####
  
  for (i in seq(clust.res.low,clust.res.high,clust.res.bin)){
    object.merge <- FindClusters(object.merge, resolution = i, algorithm = 1)
    if(do.cite.seq==TRUE){
      object.merge <- FindClusters(object.merge, 
                                   graph.name="Protein_snn",
                                   resolution = i, 
                                   algorithm = 1)
    }
  }
  print("Clustering successful!")
  
  
  ## create plots ####
  
  n <- 60
  qual.col.pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  qual.col.pals = qual.col.pals[c(7,6,2,1,8,3,4,5),]
  cols = unlist(mapply(brewer.pal, qual.col.pals$maxcolors, rownames(qual.col.pals)))
  
  
  grobsList = list()
  if(draw.tsne){
    p1 <- DimPlot(object.merge, 
                  reduction = "tsne",group.by = "orig.ident", 
                  repel = TRUE,pt.size=0.02) + 
      theme_classic() + 
      scale_color_manual(values=cols) + 
      theme(legend.position="top", legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, 
                                   override.aes = list(size=1, alpha = 1))) +         
      ggtitle("RNA TSNE")
    
    grobsList[[length(grobsList)+1]] <- p1
    print("Added RNA TSNE")
    print(length(grobsList))
    
  }
  if(draw.umap){
    p2 <- DimPlot(object.merge, 
                  reduction = "umap", 
                  group.by = "orig.ident", 
                  repel = TRUE, pt.size=0.02) + 
      theme_classic() + scale_color_manual(values=cols) + 
      theme(legend.position="top", legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, 
                                   override.aes = list(size=1, alpha = 1))) + 
      ggtitle("RNA UMAP")
    
    grobsList[[length(grobsList)+1]] <- p2
    print("Added RNA UMAP")
    print(length(grobsList))
    
  }
  
  ### CITEseq Figures
  if(draw.tsne & do.cite.seq){ 
    p3 <- DimPlot(object.merge, 
                  reduction = "protein_tsne", 
                  group.by = "orig.ident", 
                  repel = TRUE, pt.size=0.02) + 
      theme_classic() + scale_color_manual(values=cols) + 
      theme(legend.position="top", legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, 
                                   override.aes = list(size=1, alpha = 1))) + 
      ggtitle("Antibody TSNE")
    
    grobsList[[length(grobsList)+1]] <- p3
    print("Added Antibody TSNE")
    print(length(grobsList))
    
  }
  if(draw.umap & do.cite.seq==TRUE){ 
    p4 <- DimPlot(object.merge, 
                  reduction = "protein_umap", 
                  group.by = "orig.ident", 
                  repel = TRUE, pt.size=0.02) + 
      theme_classic() + 
      scale_color_manual(values=cols) + 
      theme(legend.position="top", 
            legend.text=element_text(size=5)) +
      guides(colour = guide_legend(ncol=3, 
                                   override.aes = list(size=1, alpha = 1))) + 
      ggtitle("Antibody UMAP")
    grobsList[[length(grobsList)+1]] <- p4
    print("Added Antibody UMAP")
    print(length(grobsList))
    
  }
  
  
  ## create Figure output ####
  
  grobs <- arrangeGrob(grobs=grobsList,ncol=n)
  
  
  return(list(object=object.merge,plot=grobs))
}

