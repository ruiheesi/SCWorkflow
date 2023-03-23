#' @title Annotating cell types using SingleR module
#' @description SingleR is an automatic annotation method for single-cell
#' RNA sequencing (scRNAseq) data (Aran et al. 2019). Given a reference dataset
#' of samples (single-cell or bulk) with known labels, it labels new cells
#' from a test dataset based on similarity to the reference.
#' @details This function is Step 5 of the basic Single-Cell RNA-seq workflow.
#' It is the starting point for downstream visualization, subsetting, and
#' analysis. It takes a combined seurat object as input, such as the one created
#' by the Combined&Renormalized function at the end of the Filter&QC Path
#'
#' @param object Object of class Seurat (your combined Seurat Object
#' after PC reduction has been performed)
#' @param species The species of your samples ("Human" or "Mouse").
#' Default is "Mouse"
#' @param reduction.type  Select the kind of clustering visualization you would
#' like to use to visualize the cell type results ("umap", "tsne", "pca")
#' @param legend.dot.size the size of the colored dots on your chart legend.
#' Default is 2
#' @param do.finetuning Performs the SingleR fine-tuning function.
#' Default is FALSE

#'
#' @import Seurat
#' @import gridExtra
#' @import tools
#' @import grid
#' @import gridBase
#' @import cowplot
#' @import ggplot2
#' @import RColorBrewer
#' @import magrittr
#' @import SingleR
#' @import celldex
#'
#'
#' @export
#'
#' @return a Seurat object with additional metadata



annotateCellTypes <- function(object,
                              species = "Mouse",
                              reduction.type = "umap",
                              legend.dot.size = 2,
                              do.finetuning = FALSE) {
  ## -------------------------------- ##
  ## Functions                        ##
  ## -------------------------------- ##
  
  
  .annotations <- function(so) {
    so.counts = GetAssayData(object = so)[, colnames(x = so)]
    if (species == "Human") {

      #HPCA block
      HPCA <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = HPCA,
        labels = HPCA$label.main
      )
      so[["HPCA_main"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = HPCA,
        labels = HPCA$label.fine
      )
      so[["HPCA"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
      #BP_encode block
      BP <- celldex::BlueprintEncodeData(ensembl = FALSE)
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = BP,
        labels = BP$label.main
      )
      so[["BP_encode_main"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = BP,
        labels = BP$label.fine
      )
      so[["BP_encode"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
    }
    if (species == "Mouse") {
      
      #mouseRNAseq block
      mousernaseq <- celldex::MouseRNAseqData(ensembl = FALSE)
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = mousernaseq,
        labels = mousernaseq$label.main
      )
      so[["mouseRNAseq_main"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = mousernaseq,
        labels = mousernaseq$label.fine
      )
      so[["mouseRNAseq"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
      #ImmGen block
      immgen <- celldex::ImmGenData(ensembl = FALSE)
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = immgen,
        labels = immgen$label.main
      )
      so[["immgen_main"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
      singler = SingleR(
        test = so.counts,
        genes = 'de',
        clusters = so@meta.data$seurat_clusters,
        fine.tune = do.finetuning,
        ref = immgen,
        labels = immgen$label.fine
      )
      so[["immgen"]] <-
        singler$labels[match(rownames(so[[]]), rownames(singler))]
      
    }
    return(so)
  }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  
  # Getting Cluster List from MetaData of Seurat Object:
  clusterList <- object@meta.data$seurat_clusters
  
  # Running Annotation:
  object <- .annotations(object)
  print("done")
  
  # Assigning Colors (choice of datasets depends on species):
  if (species == "Human") {
    numColors = max(length(unique(object@meta.data$BP_encode_main)), length(unique(object@meta.data$HPCA_main)))
  } else {
    numColors = max(length(unique(object@meta.data$mouseRNAseq_main)), length(unique(object@meta.data$immgen_main)))
  }
  colpaired = colorRampPalette(brewer.pal(12, "Paired"))
  cols = c(
    "#e6194B",
    "#3cb44b",
    "#4363d8",
    "#f58231",
    "#911eb4",
    "#42d4f4",
    "#f032e6",
    "#bfef45",
    "#fabebe",
    "#469990",
    "#e6beff",
    "#9A6324",
    "#800000",
    "#aaffc3",
    "#808000",
    "#000075",
    colpaired(numColors)
  )
  
  # Creating plots (choice of datasets depends on species):
  if (species == "Human") {
    p1 = DimPlot(object, reduction = reduction.type, group.by = "HPCA_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                     colour = guide_legend(ncol = 4)) + ggtitle("HPCA Main Cell Type Annotations")
    p2 = DimPlot(object, reduction = reduction.type, group.by="BP_encode_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                          colour = guide_legend(ncol = 4)) + ggtitle("BP Encode Main Cell Type Annotations")
  } else {
    p1 = DimPlot(object, reduction = reduction.type, group.by = "immgen_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                       colour = guide_legend(ncol = 4)) + ggtitle("Immgen Main Cell Type Annotations")
    p2 = DimPlot(object, reduction=reduction.type, group.by="mouseRNAseq_main") +
         scale_color_manual(values = cols) +
         theme(legend.position = "top") +
         guides(override.aes = list(size = legend.dot.size),
                                                                                                                                                            colour = guide_legend(ncol = 4)) + ggtitle("Mouse RNAseq Main Cell Type Annotations")
  }
  
  #  Adjusting Seurat Object MetaData:
  #  print(plot_grid(p1,p2,nrow=1))
  object@meta.data$Barcode <- rownames(object@meta.data)
  object@meta.data$sample_name <- object@meta.data$orig.ident
  object@meta.data$sample_name <-
    gsub("-", "_", object@meta.data$sample_name)
  object@meta.data$seurat_clusters <- NULL
  rownames(object@meta.data) <- object@meta.data$Barcode
  
  slot(object, "commands") <- list()
  #  cat("\nSingleR Object Checksum:\n")
  #  print(digest::digest(object))
  
  # Returning Seurat Object and 2 plots:
  return(list(
    "object" = object,
    "p1" = p1,
    "p2" = p2
  ))
}
