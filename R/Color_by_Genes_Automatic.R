#' @title Color by Gene List
#' @description Returns a panel of reduction plots colored by marker expression
#' @details Takes in a gene table inputted by the user, displays a panel of
#'          tsne, umap, or pca colored by marker expression. The panel will be 
#'          organized in a similar format as the gene table, but with the 
#'          omission of genes not found in the data
#'          
#' @param object Seurat-class object
#' @param samples.subset List of samples to subset the data by
#' @param samples.to.display List of samples to depict on dimension plot, 
#'                           samples not in the list would be colored gray in 
#'                           the background
#' @param marker.table Table of marker genes for each celltype 
#'                    (column names of the table), append "_prot" or "_neg" for 
#'                    proteins or negative markers
#' @param cells.of.interest Celltypes from geneset_dataframe to screen for
#' @param protein.presence Set to TRUE if protein markers are used
#' @param assay Assay to extract gene expression data from (Default: "SCT")
#' @param reduction.type Choose among tsne, umap, and pca (Default: "umap")
#' @param point.transparency Set to lower values for more see through points on 
#'                           dimension plot (Default: 0.5)
#' @param point.shape Change the shape of points for between visualization
#'                    (Default: 16)
#' @param cite.seq Set to TRUE to use CITE-seq embedding for dimension reduction
#'
#' @import Seurat
#' @import tidyverse
#' @import gridExtra
#' @import ggpubr
#' @import ggplot2
#'
#' @export
#' @example Do not run: colorByMarkerTable(object = seurat,
  #'                                       samples.subset = c("mouse1","mouse2),
  #'                                       samples.to.display = c("mouse1"),
  #'                                       marker.table = immuneCellMarkers,
  #'                                       cells.of.interest = c("CD4","Treg")
  #'                                       )

#' @return arranged grob of dimension reduction plots colored by individual 
#'         marker expression

colorByMarkerTable <- function(object,
                               samples.subset,
                               samples.to.display,
                               marker.table,
                               cells.of.interest,
                               protein.presence = FALSE,
                               assay = "SCT",
                               reduction.type = "umap",
                               point.transparency = 0.5,
                               point.shape = 16,
                               cite.seq = FALSE) {

  # Error Messages
  if (!assay %in% Assays(object)) {
    stop("assay type not found in seurat")
  } else if (!reduction.type %in% Reductions(object)) {
    stop("reduction type not found in seurat")
  }
  
  # Functions
  .plotMarkers <- function(markers) {
    if (is.na(markers) == TRUE) {
      g <- ggplot() + theme_void()
      return(g)
    } else {
      markers.mat = object.sub[[assay]]@scale.data[markers, ]
      markers.quant = quantile(markers.mat[markers.mat > 1], probs = c(.1, .5, .90))
      markers.mat[markers.mat > markers.quant[3]] = markers.quant[3]
      markers.mat[markers.mat < markers.quant[1]] = 0
      
      if (!(cite.seq)) {
        if (reduction.type == "tsne") {
          p1 <- DimPlot(object.sub, reduction = "tsne", group.by = "ident")
          clusmat = data.frame(
            umap1 = p1$data$tSNE_1,
            umap2 = p1$data$tSNE_2,
            markers = markers.mat,
            ident = as.factor(p1$data$ident)
          )
        }
        else if (reduction.type == "umap") {
          p1 <- DimPlot(object.sub, reduction = "umap", group.by = "ident")
          clusmat = data.frame(
            umap1 = p1$data$UMAP_1,
            umap2 = p1$data$UMAP_2,
            markers = markers.mat,
            ident = as.factor(p1$data$ident)
          )
        }
        else{
          p1 <- DimPlot(object.sub, reduction = "pca", group.by = "ident")
          clusmat = data.frame(
            umap1 = p1$data$PC_1,
            umap2 = p1$data$PC_2,
            markers = markers.mat,
            ident = as.factor(p1$data$ident)
          )
        } #if CITEseq is chosen then:
      } else {
        if (reduction.type == "tsne") {
          p1 <-
            DimPlot(object.sub, reduction = "protein_tsne", group.by = "ident")
          clusmat = data.frame(
            umap1 = p1$data$protein_tsne_1,
            umap2 = p1$data$protein_tsne_2,
            markers = markers.mat,
            ident = as.factor(p1$data$ident)
          )
        }
        else if (reduction.type == "umap") {
          p1 <-
            DimPlot(object.sub, reduction = "protein_umap", group.by = "ident")
          clusmat = data.frame(
            umap1 = p1$data$protein_umap_1,
            umap2 = p1$data$protein_umap_2,
            markers = markers.mat,
            ident = as.factor(p1$data$ident)
          )
        }
        else{
          p1 <- DimPlot(object.sub, reduction = "protein_pca", 
                        group.by = "ident")
          clusmat = data.frame(
            umap1 = p1$data$protein_pca_1,
            umap2 = p1$data$protein_pca_2,
            markers = markers.mat,
            ident = as.factor(p1$data$ident)
          )
        }
      }
      
      # Samples caption
      samples.caption <-
        paste(samples.to.display,
              sep = "",
              collapse = "\n")
      final_caption <-
        paste(
          "Samples Displayed: ",
          samples.caption,
          sep = "",
          collapse = "\n"
        )
      
      clusmat <-
        mutate(clusmat,
               sample.markers = clusmat$markers * grepl(paste(
                 samples.to.display, collapse = "|"), clusmat$ident))
      
      clusmat %>% dplyr::arrange(sample.markers) -> clusmat
      if (grepl("_neg", markers) == TRUE) {
        clusmat %>% dplyr::arrange(desc(sample.markers)) -> clusmat
        g <- ggplot(clusmat, aes(
          x = umap1,
          y = umap2,
          group = ident
        )) +
          theme_bw() +
          theme(legend.title = element_blank()) +
          ggtitle(markers) +
          geom_point(
            aes(color = sample.markers, shape = ident),
            alpha = point.transparency,
            shape = point.shape,
            size = 1
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text = element_text(size = rel(0.5))
          ) +
          scale_color_gradient(
            limits = c(0, markers.quant[3]),
            low = "lightgrey",
            high = "red"
          ) +
          xlab("umap-1") + ylab("umap-2")
        return(g)
      } else {
        clusmat %>% dplyr::arrange(sample.markers) -> clusmat
        g <- ggplot(clusmat, aes(
          x = umap1,
          y = umap2,
          group = ident
        )) +
          theme_bw() +
          theme(legend.title = element_blank()) +
          ggtitle(markers) +
          geom_point(
            aes(color = sample.markers, shape = ident),
            alpha = point.transparency,
            shape = point.shape,
            size = 1
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text = element_text(size = rel(0.5))
          ) +
          scale_color_gradient(
            limits = c(0, markers.quant[3]),
            low = "lightgrey",
            high = "red"
          ) +
          xlab("umap-1") + ylab("umap-2")
        return(g)
      }
    }
  }
  
  # Main Code Block
  if (length(samples.subset) == 0) {
    samples.subset = unique(object@meta.data$sample.name)
  }
  
  if ("active.ident" %in% slotNames(object)) {
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name) = names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples.subset)
  } else {
    sample.name = as.factor(object@meta.data$orig.ident)
    names(sample.name) = names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- sample.name
    object.sub = subset(object, ident = samples.subset)
  }
  
  marker.table <- marker.table[cells.of.interest]
  
  # Remove columns with all missing values
  present.marker.ls <- list()
  
  for (celltype in colnames(marker.table)) {
    print(names(marker.table[celltype]))
    present = lapply(marker.table[[celltype]], function(x)
      x %in% rownames(object.sub$SCT@scale.data))
    absent.genes = unlist(marker.table[[celltype]])[present == FALSE]
    present.genes = unlist(marker.table[[celltype]])[present == TRUE]
    print(paste0("Genes not present: ", paste0(absent.genes, collapse = ",")))
    print(paste0("Genes present: ", paste0(present.genes, collapse = ",")))
    
    if (length(present.genes) == 0) {
      print(paste0(
        names(marker.table[celltype]),
        " genes were not found in object and will not be analyzed"
      ))
    } else {
      present.marker.ls[[celltype]] <- present.genes
    }
  }
  
  # Padd processed list containing only the present genes
  padded.ls <- lapply(present.marker.ls, `length<-`, 
                      max(lengths(present.marker.ls)))
  markers.from.list <- do.call(cbind, padded.ls)
  
  # Recognize any markers designated as proteins and add protein expression 
  if (protein.presence) {
    protein.markers <- markers.from.list[grepl("_prot", markers.from.list)]
    
    protein.orig <- gsub("_prot.*", "", protein.markers)
    
    protein.name <- paste(protein.orig,
                                  "_prot", sep = "")
    
    protein.array <- list()
    for (i in seq_along(protein.orig)) {
      protein.array[[i]] <-
        object.sub@assays$Protein[protein.orig[i], ]
      rownames(protein.array[[i]]) <- protein.name[i]
    }
    protein.array.comp <- do.call(rbind, protein.array)
    object.sub@assays$SCT@scale.data <-
      rbind(object.sub@assays$SCT@scale.data, protein.array.comp)
  }
  
  # Add negative/low identifiers
  neg.markers.names <-
    markers.from.list[grepl("_neg", markers.from.list)]
  orig.markers <- gsub("_neg.*", "", neg.markers.names)
  
  # Append neg.markers.names to rownames of object.sub
  neg.markers.ls <- list()
  
  # Calculate adjusted expression for negative markers
  for (i in seq_along(orig.markers)) {
    if (orig.markers[i] %in% rownames(object.sub@assays$SCT@scale.data)) {
      # Format the data so that it can rbinded with object$SCT@scale.data
      neg.markers.ls[[i]] <-
        t(matrix(
          max(object.sub@assays$SCT@scale.data[orig.markers[i], ]) - 
            object.sub@assays$SCT@scale.data[orig.markers[i], ]
        ))
      row.names(neg.markers.ls[[i]]) <- neg.markers.names[i]
      colnames(neg.markers.ls[[i]]) <-
        colnames(object.sub@assays$SCT@scale.data)
      
      # Append new Negative/low marker (w Expression Count) to object slot
      object.sub@assays$SCT@scale.data <-
        rbind(object.sub@assays$SCT@scale.data, neg.markers.ls[[i]])
    } else {
      print(
        paste(
          orig.markers[i],
          " is not found in Seurat, cannot calculate negative expression",
          sep = ""
        )
      )
    }
  }
  
  markers.present = unlist(markers.from.list)
  
  if (!length(markers.present) > 0) {
    print("No markers found in dataset")
    return(NULL)
  }
  
  # Create list for storing color by gene plots of each celltype column
  gg.storage <- list()
  
  for (cell in colnames(markers.from.list)) {
    title <- cell
    
    markers.to.analyze <- as.character(markers.from.list[, cell])
    
    grob <- lapply(markers.to.analyze, function(x) .plotMarkers(x))
    
    gg.storage[[cell]] <-
      gridExtra::arrangeGrob(
        grobs = grob,
        ncol = 1,
        newpage = F,
        as.table = F,
        top = text_grob(title, size = 15, face = "bold")
      )
    
  }
  
  final.figures <-
    do.call(arrangeGrob, c(gg.storage, ncol = ncol(markers.from.list)))
  
  return(final.figures)
}
