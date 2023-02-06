# This code comes from NIDAP 'Color by Genes Automatic [scRNA-seq][CCBR]' code template
# Template Manager https://nidap.nih.gov/workspace/vector/templates/ri.vector.main.template.d71ed4e6-a25d-4f66-a186-27c00a50a703
# Documentation https://nidap.nih.gov/workspace/notepad/view/ri.notepad.main.notepad.8101be48-85a5-40ec-8fa8-9fe28e5d31cb

#' @title Create gene expression annotated dimension plots based on user inputted table
#' @description Returns a panel of dimension plots colored by individual gene expression
#' @details Takes in a list of genes inputted by the user, displays gene expression information on tsne, umap, or pca
#' 
#' @param SO Seurat-class object
#' @param samples_to_include List of samples to subset the data by
#' @param samples_to_display List of samples to depict on dimension plot, samples not in the list would be colored gray in the background
#' @param marker_list Table of marker genes for each celltype (column names of the table), append "_prot" or "_neg" for proteins or negative markers
#' @param cells_of_interest Vector of celltypes from geneset_dataframe to screen for
#' @param protein_presence Set to TRUE if there are CITE-seq markers in geneset_dataframe
#' @param assay Name of the assay to extract gene expression data from
#' @param reduction_type Choose among tsne, umap, and pca
#' @param point_transparency Set to lower value for more see through points on dimension plot
#' @param point_shape Change the shape of points for between visualization
#' @param number_of_rows Set the number of rows to be displayed on the compiled gene expression dimension plot
#' @param doCiteSeq Set to TRUE if using Cite-seq data
#' 
#' @import Seurat 
#' @import tidyverse
#' @import gridExtra
#' @import ggpubr
#' @import ggplot2
#'   
#' @export

#' @return compiled dimension plots of markers, in the same layout as the user-inputted marker table

color_by_genes <- function(SO, 
                             samples_to_include,
                             samples_to_display,
                             marker_list,
                             cells_of_interest,
                             protein_presence = FALSE,
                             assay = "SCT",
                             reduction_type = "umap",
                             point_transparency = 0.5,
                             point_shape = 16,
                             number_of_rows = 0,
                             doCiteSeq = FALSE){
  
  ## -------------- ##
  ## Error Messages ##
  ## -------------- ##
  
  if(!assay %in% Assays(SO)){
    stop("assay type not found in seurat")
  } else if (!reduction_type %in% Reductions(SO)){
    stop("reduction type not found in seurat")
  }
    
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  plotMarkers <- function(Markers){
    if (is.na(Markers) == TRUE) {
      g <- ggplot() + theme_void()
      return(g)
    } else {
      Markers.mat=SO.sub[[assay]]@scale.data[Markers,]
      Markers.quant=quantile(Markers.mat[Markers.mat>1],probs=c(.1,.5,.90))
      Markers.mat[Markers.mat>Markers.quant[3]]=Markers.quant[3]
      Markers.mat[Markers.mat<Markers.quant[1]]=0
      
      if (!(doCiteSeq)) {
        if(reduction_type == "tsne"){
          p1 <- DimPlot(SO.sub, reduction = "tsne", group.by = "ident")
          clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, Markers=Markers.mat, ident = as.factor(p1$data$ident))
        }
        else if(reduction_type == "umap"){
          p1 <- DimPlot(SO.sub, reduction = "umap", group.by = "ident")
          clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, Markers=Markers.mat,ident = as.factor(p1$data$ident))
        }
        else{
          p1 <- DimPlot(SO.sub, reduction = "pca", group.by = "ident")
          clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, Markers=Markers.mat,ident = as.factor(p1$data$ident))
        } #if CITEseq is chosen then:
      } else {
        if(reduction_type == "tsne"){
          p1 <- DimPlot(SO.sub, reduction = "protein_tsne", group.by = "ident")
          clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, Markers=Markers.mat,ident = as.factor(p1$data$ident))
        }
        else if(reduction_type == "umap"){
          p1 <- DimPlot(SO.sub, reduction = "protein_umap", group.by = "ident")
          clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, Markers=Markers.mat,ident = as.factor(p1$data$ident))
        }
        else{
          p1 <- DimPlot(SO.sub, reduction = "protein_pca", group.by = "ident")
          clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, Markers=Markers.mat,ident = as.factor(p1$data$ident))
        }
      }
      
      # Samples caption
      samples_displayed_message <- paste(samples_to_display, sep = "", collapse = "\n")
      final_caption <- paste("Samples Displayed: ", samples_displayed_message, sep = "", collapse = "\n")
      
      clusmat <- mutate(clusmat, sample_Markers = clusmat$Markers * grepl(paste(samples_to_display, collapse = "|"), clusmat$ident))
      
      clusmat %>% dplyr::arrange(sample_Markers) -> clusmat
      if (grepl("_neg",Markers) == TRUE){
        
        clusmat %>% dplyr::arrange(desc(sample_Markers)) -> clusmat
        g <- ggplot(clusmat, aes(x=umap1, y=umap2, group=ident)) +
          theme_bw() +
          theme(legend.title=element_blank()) +
          ggtitle(Markers) +
          geom_point(aes(color=sample_Markers, shape=ident),alpha=point_transparency,shape=point_shape, size=1) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(),legend.text=element_text(size=rel(0.5)) )+
          scale_color_gradient(limits = c(0, Markers.quant[3]),low = "lightgrey", high = "red") +
          xlab("umap-1") + ylab("umap-2")
        return(g)
      } else {
        clusmat %>% dplyr::arrange(sample_Markers) -> clusmat
        g <- ggplot(clusmat, aes(x=umap1, y=umap2, group=ident)) +
          theme_bw() +
          theme(legend.title=element_blank()) +
          ggtitle(Markers) +
          geom_point(aes(color=sample_Markers, shape=ident),alpha=point_transparency,shape=point_shape, size=1) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(),legend.text=element_text(size=rel(0.5)) )+
          scale_color_gradient(limits = c(0, Markers.quant[3]),low = "lightgrey", high = "red") +
          xlab("umap-1") + ylab("umap-2")
        return(g)}
    }
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  if (length(samples_to_include) == 0) {
    samples_to_include = unique(SO@meta.data$sample_name)
  }
  
  if("active.ident" %in% slotNames(SO)){
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples_to_include)
  } else {
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples_to_include)
  }
  
  marker_list <- marker_list[cells_of_interest]
  
  # Remove columns with all missing values
  Present_Markers_ls <- list()
  
  for (celltype in colnames(marker_list)) {
    print(names(marker_list[celltype]))
    present=lapply(marker_list[[celltype]], function(x) x %in% rownames(SO.sub$SCT@scale.data)) 
    absentgenes = unlist(marker_list[[celltype]])[present==FALSE]
    presentgenes = unlist(marker_list[[celltype]])[present==TRUE]
    print(paste0("Genes not present: ",paste0(absentgenes,collapse=",")))
    print(paste0("Genes present: ",paste0(presentgenes,collapse=",")))
    
    if(length(presentgenes) == 0){
      print(paste0(names(marker_list[celltype]), " genes were not found in SO and will not be analyzed"))
    } else {Present_Markers_ls[[celltype]] <- presentgenes}
  }
  
  # Padd processed list containing only the present genes
  padded_list <- lapply(Present_Markers_ls, `length<-`, max(lengths(Present_Markers_ls)))
  Markers_from_list <- do.call(cbind, padded_list)
  
  # Recognize any markers designated as proteins and insert protein expression data into slot where plots are created
  if (protein_presence){
    protein_markers <- Markers_from_list[grepl("_prot",Markers_from_list)]
    
    protein_orig_markers <- gsub("_prot.*","",protein_markers)
    
    protein_markers_name <- paste(protein_orig_markers,
                                  "_prot", sep = "")
    
    i = 0
    protein_array <- list()
    for (i in seq_along(protein_orig_markers)){
      protein_array[[i]] <- SO.sub@assays$Protein[protein_orig_markers[i],]
      rownames(protein_array[[i]]) <- protein_markers_name[i]
    }
    protein_array_comp <- do.call(rbind,protein_array)
    SO.sub@assays$SCT@scale.data <- rbind(SO.sub@assays$SCT@scale.data,protein_array_comp)
  }
  
  # Add negative/low identifiers
  neg_markers_names <- Markers_from_list[grepl("_neg",Markers_from_list)]
  orig_markers <- gsub("_neg.*","",neg_markers_names)
  
  # Append neg_markers_names to rownames of SO.sub
  neg_markers_list <- list()
  
  # Calculate adjusted expression for negative markers
  for (i in seq_along(orig_markers)){
    if (orig_markers[i] %in% rownames(SO.sub@assays$SCT@scale.data)){
      # Format the data so that it can rbinded with SO$SCT@scale.data
      neg_markers_list[[i]] <- t(matrix(max(SO.sub@assays$SCT@scale.data[orig_markers[i],]) - SO.sub@assays$SCT@scale.data[orig_markers[i],]))
      row.names(neg_markers_list[[i]]) <- neg_markers_names[i]
      colnames(neg_markers_list[[i]]) <- colnames(SO.sub@assays$SCT@scale.data)
      
      # Append new Negative/low marker (w Expression Count) to SO slot
      SO.sub@assays$SCT@scale.data <- rbind(SO.sub@assays$SCT@scale.data, neg_markers_list[[i]]) 
    } else {
      print(paste(orig_markers[i], " is not found in Seurat, cannot calculate negative expression", sep = ""))
    }} 
  
  Markers_present = unlist(Markers_from_list)

  if(!length(Markers_present)>0){
    print("No Markers found in dataset")
    return(NULL)
  }
  
  # Create list for storing color by gene plots of each celltype column
  gg_storage <- list()
  
  for (cell in colnames(Markers_from_list)){
    
    title <- cell
    
    markers_to_analyze <- as.character(Markers_from_list[,cell])
    
    grob <- lapply(markers_to_analyze,function(x) plotMarkers(x))
    gg_storage[[cell]] <- gridExtra::arrangeGrob(grobs=grob,ncol = 1,newpage=F, as.table = F, top = text_grob(title,size = 15, face = "bold"))
    
  }
  
  final_figures <- do.call(arrangeGrob, c(gg_storage, ncol = ncol(Markers_from_list)))

  return(final_figures)
}
  