  # This code comes from NIDAP 'ModuleScore [scRNA-seq][CCBR]' code template
  # Template Manager https://nidap.nih.gov/workspace/vector/templates/ri.vector.main.template.10cf059e-0bd0-4a1c-9a3a-65b8688dab23
  # Documentation https://nidap.nih.gov/workspace/notepad/view/ri.notepad.main.notepad.1d5957fa-ef90-4473-b627-989b9620278f

  #' @title Calculate average gene expression (subtracted by aggregated control feature set) and label cells according to ModuleScores
  #' @description Returns Seurat-class object with metadata containing ModuleScores and Likely_CellType calls
  #' @details Analyzed features are binned based on averaged expression; control features are randomly selected from each bin.
  #' 
  #' @param SO Seurat-class object
  #' @param sample_names List of samples to subset the data by
  #' @param sample_to_display List of samples to depict on dimension plot, samples not in the list would be colored gray in the background
  #' @param geneset_dataframe Table of marker genes for each celltype (column names of the table), append "_prot" or "_neg" for proteins or negative markers
  #' @param proteins_presence Set to TRUE if there are CITE-seq markers in geneset_dataframe
  #' @param celltypes_to_analyze Vector of celltypes from geneset_dataframe to screen for
  #' @param manual_threshold Specify bimodal thresholds for cell classification, should be of the same length as celltypes_to_analyze vector
  #' @param general_class Base population of cells to classify
  #' @param multi_level_class Toggle to TRUE if there are subpopulations of cells you want to screen for
  #' @param levels_dataframe Table of subpopulation levels (column names of the table) and parent-child information (e.g. Tcells-CD4, Tcells-CD8)
  #' @param reduction Choose among tsne, umap, and pca
  #' @param nbins Number of bins for storing control features and analyzing average expression
  #' @param gradient_density_font_size Set size of axis labels on gradient density plot of ModuleScore distribution
  #' @param violinplot_font_size Set size of axis labels on violin plot of ModuleScore distribution
  #' @param step_size Set step size of distribution plots
  #' @param image_width Set Image width (pixel units)
  #' @param image_height Set Image height (pixel units)
  #' @param image_type Choose between png and svg
  #' 
  #' @import Seurat
  #' @import tidyverse
  #' @import gridExtra
  #' @import quantmod
  #' @import grid
  #' @import data.table
  #' @import utils
  #' @importFrom dplyr select
  #'   
  #' @export

  #' @return Annotated dimension plot with ModuleScore distribution of cell marker gene, Seurat Object with cell classification metadata

ModuleScore <- function(SO, 
                         sample_names, 
                         sample_to_display, 
                         geneset_dataframe, 
                         proteins_presence = FALSE, 
                         celltypes_to_analyze,
                         manual_threshold = c(0), 
                         general_class, 
                         multi_level_class = FALSE, 
                         levels_dataframe, 
                         reduction = "tsne",
                         nbins = 24,
                         gradient_density_font_size = 6,
                         violinplot_font_size = 6,
                         step_size = 0.1
                         ){

  ## ---------------- ##
  ## Helper Functions ##
  ## ---------------- ##
  
  Predict_Cell_from_ModScore <- function(ModScore_Metadata,manual_threshold,rejection){
    
    thres_ls <- list()
    for (i in 1:ncol(ModScore_Metadata)){
      thres_ls[[i]]<- rep(manual_threshold[i],nrow(ModScore_Metadata))
    }
    thres_df <- data.frame(matrix(unlist(thres_ls),nrow = nrow(ModScore_Metadata)))
    
    thres_filter <- ModScore_Metadata > thres_df
    ModScore_Metadata_post_thres_filter <- ModScore_Metadata * thres_filter
    
    # Find column number with highest modscore
    max_col_vector <- max.col(ModScore_Metadata_post_thres_filter)
    
    # If a row contains all zeroes, they will be labeled with unknown
    all_zero_filter <- as.integer(!apply(ModScore_Metadata_post_thres_filter, 1, function(find_zero_rows) all(find_zero_rows == 0)))
    
    # Final filtering: 
    final_filter <- (max_col_vector * all_zero_filter) + 1
    
    # Original names appended to "unknown" classification for cells with ModScores below threshold
    appended_names <- c(rejection, names(ModScore_Metadata))
    
    # Added the names into a Likely_CellType Column
    dupl_data <- ModScore_Metadata
    dupl_data[,"Likely_CellType"] <- appended_names[final_filter]
    return(dupl_data)
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  geneset_for_module_scores <- unlist(geneset_dataframe)
  
  # Create a Barcode column if none is detected - required for matching cell calls in downstream steps
  if (!"Barcode" %in% colnames(SO@meta.data)){
    SO@meta.data$Barcode <- rownames(SO@meta.data)
  }
  
  if (length(sample_names) == 0) {
    sample_names = unique(SO@meta.data$sample_name)
  }
  
  colnames(SO@meta.data) <- gsub("orig_ident","orig.ident",colnames(SO@meta.data))
  
  if("active.ident" %in% slotNames(SO)){
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = sample_names)
  } else {
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = sample_names)
  } 
  
  # Remove original unprocessed SO
  rm(SO)
  
  # Adding protein marker expression
  if (proteins_presence){
    protein_markers <- geneset_for_module_scores[grepl("_prot",geneset_for_module_scores)]
    
    protein_orig_markers <- gsub("_prot.*","",protein_markers)
    
    protein_markers_name <- paste(protein_orig_markers,
                                  "_prot", sep = "")
    
    prot_indx = 0
    protein_array <- list()
    for (prot_indx in seq_along(protein_orig_markers)){
      protein_array[[prot_indx]] <- SO.sub@assays$Protein[protein_orig_markers[prot_indx],]
      rownames(protein_array[[prot_indx]]) <- protein_markers_name[prot_indx]
    }
    protein_array_comp <- do.call(rbind,protein_array)
    SO.sub@assays$SCT@data <- rbind(SO.sub@assays$SCT@data,protein_array_comp)
  }
  
  # Recognize any negative markers in marker list
  neg_markers_names <- geneset_for_module_scores[grepl("_neg",geneset_for_module_scores)]
  orig_markers <- gsub("_neg.*","",neg_markers_names)
  
  # Retrieve markers found in counts data
  orig_markers <- orig_markers[orig_markers %in% rownames(SO.sub@assays$SCT@data)]
  
  neg_markers_list <- list()
  
  # Calculate adjusted expression for negative markers
  for (neg_indx in seq_along(orig_markers)){
    
    # Format the data so that it can rbinded with SO$SCT@scale.data
    neg_markers_list[[neg_indx]] <- t(matrix(max(SO.sub@assays$SCT@data[orig_markers[neg_indx],]) - SO.sub@assays$SCT@data[orig_markers[neg_indx],]))
    row.names(neg_markers_list[[neg_indx]]) <- neg_markers_names[neg_indx]
    colnames(neg_markers_list[[neg_indx]]) <- colnames(SO.sub@assays$SCT@data)
    
    # Append new Negative/low marker (w Expression Count) to SO slot
    SO.sub@assays$SCT@data <- rbind(SO.sub@assays$SCT@data, neg_markers_list[[neg_indx]]) 
  }
  
  # Retrive markers from list
  marker = select(geneset_dataframe, celltypes_to_analyze)
  marker.list = as.list(marker)
  
  # Error checking and messages for when threshold vector is zero or of non-matching length to number of celltypes to analyze
  if (sum(unlist(marker.list) %in% rownames(SO.sub@assays$SCT@data)) == 0){
    stop("No genes from list was found in data")
  }
  
  if (length(manual_threshold) != length(celltypes_to_analyze)){
    if (sum(manual_threshold) == 0){
      manual_threshold <- rep(0, length(celltypes_to_analyze))
      print("Manual threshold set to zero - outputing preliminary data")
    } else {
      stop("Manual threshold length does not match number of celltypes to analyze - please check manual thresholds")
    }}
  
  names(manual_threshold) <- celltypes_to_analyze
  
  # Check if all markers for particular celltype is present, if none are present - corresponding celltype will not be analyzed
  figures <- list()
  exclude_cells <- c()
  
  h = 0
  j = 1
  
  for (h in seq_along(marker.list)) {
    print(names(marker.list[h]))
    present=lapply(marker.list[[h]], function(x) x %in% rownames(SO.sub)) # apply function(x) x %in% rownames(SO.sub) to each element of marker.list
    absentgenes = unlist(marker.list[[h]])[present==FALSE];absentgenes=absentgenes[is.na(absentgenes)==F]
    presentgenes = unlist(marker.list[[h]])[present==TRUE];presentgenes=presentgenes[is.na(presentgenes)==F]
    print(paste0("Genes not present: ",paste0(absentgenes,collapse=",")))
    print(paste0("Genes present: ",paste0(presentgenes,collapse=",")))
    
    if(length(presentgenes) == 0){
      print(paste0(names(marker.list[h]), " genes were not found in SO and will not be analyzed"))
      exclude_cells[j] <- h
      j = j + 1
    }}  
  
  if (length(exclude_cells) > 0){
    marker.list <- marker.list[-exclude_cells]} else {
      marker.list <- marker.list
    }  
  
  # Calculate aggregate expression of markers with Seurat's Module Score function
  for (i in seq_along(marker.list)) { 
    SO.sub=AddModuleScore(SO.sub,marker.list[i],name = names(marker.list[i]),
                          nbin = nbins)
    
    m = 0
    
    m = paste0(names(marker.list[i]),"1")
    SO.sub@meta.data[[m]] <- scales::rescale(SO.sub@meta.data[[m]], to=c(0,1))
    
    clusid = SO.sub@meta.data[[m]] 
    
    # Calculate density data for ModScore vs Number of cells
    d <- density(clusid)
    
    # Calculate dimension reduction
    if(reduction=="tsne"){
      p1 <- DimPlot(SO.sub, reduction = "tsne", group.by = "ident") 
    } else if(reduction=="umap"){
      p1 <- DimPlot(SO.sub, reduction = "umap", group.by = "ident")
    } else { 
      p1 <- DimPlot(SO.sub, reduction = "pca", group.by = "ident")
    }
    
    if(reduction=="tsne"){
      clusmat=data.frame(ident=p1$data$ident,umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
    } else if(reduction=="umap"){
      clusmat=data.frame(ident=p1$data$ident,umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
    } else { 
      clusmat=data.frame(ident=p1$data$ident,umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.numeric(SO.sub@meta.data[[m]]))
    }
    clusmat <- mutate(clusmat, sample_clusid = clusmat$clusid * grepl(paste(sample_to_display, collapse = "|"), clusmat$ident))
    
    clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
    title=as.character(m)
    clusmat %>% dplyr::arrange(clusid) -> clusmat
    
    # Plot data onto colored dimension reduction plot
    clusid.df <- data.frame(id=SO.sub@meta.data$orig.ident,ModuleScore=SO.sub@meta.data[[m]])
    
    g <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
      theme_bw() +
      theme(legend.title=element_blank()) +
      geom_point(aes(colour=sample_clusid),alpha=0.5,shape = 20,size=1) +
      scale_color_gradientn(colours = c("blue4","lightgrey", "red"), values = scales::rescale(c(0,manual_threshold[i]/2,manual_threshold[i],(manual_threshold[i]+1)/2,1), limits = c(0, 1))) + guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
      xlab("tsne-1") + ylab("tsne-2")
    
    g1 = RidgePlot(SO.sub,features=m,group.by="orig.ident") + theme(legend.position = "none", title = element_blank(), axis.text.x = element_text(size = gradient_density_font_size)) + geom_vline(xintercept = manual_threshold[i], linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0,1,step_size))
    
    # Violin Plot
    g2 = ggplot(clusid.df,aes(x=id,y=ModuleScore)) + geom_violin(aes(fill=id)) +  theme_classic() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.title = element_blank(), panel.background = element_blank(), axis.text.x=element_blank(),legend.text=element_text(size=rel(0.8)),legend.position="top", axis.text.y = element_text(size = violinplot_font_size)) + guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) + geom_hline(yintercept = manual_threshold[i], linetype = "dashed", color = "red3") + scale_y_continuous(breaks = seq(0,1,step_size))
    
    # Color gradient density plot
    g3 = ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + xlab("ModuleScore") + ylab("Density") + geom_line() + 
      geom_segment(aes(xend = d$x, yend = 0, colour = x)) + scale_y_log10() +
      scale_color_gradientn(colours = c("blue4","lightgrey", "red"), values = scales::rescale(c(0,manual_threshold[i]/2,manual_threshold[i],(manual_threshold[i]+1)/2,1), limits = c(0, 1))) + geom_vline(xintercept = manual_threshold[i], linetype = "dashed", color = "red3") + geom_vline(xintercept = manual_threshold[i], linetype = "dashed", color = "red3") + scale_x_continuous(breaks = seq(0,1,step_size)) + theme(legend.title = element_blank(), axis.text.x = element_text(size = 6))
    
    # Set title for arrangeGrob final figure
    figures[[i]] = arrangeGrob(g,g1,g2,g3, ncol=2, top=textGrob(names(marker.list[i]), gp = gpar(fontsize = 14, fontface = "bold")))
    }
  
  # Get rid of "1" at the end of MS columns
  colnames(SO.sub@meta.data)[colnames(SO.sub@meta.data) %in% paste0(names(marker.list),1)] <- names(marker.list)
  
  # Classification of general cell class
  
  # Analyze and plot only the class of cells found in the metadata
  general_class <- general_class[general_class %in% colnames(SO.sub@meta.data)]
  
  # Subset the columns of the metadata containing module scores only
  SO_Trunc_Metadata_General <- SO.sub@meta.data[general_class]
  
  General_thres_vec <- manual_threshold[general_class]
  
  # See if elements in each ModScore column exceeds CellType threshold, set elements below threshold to zero. Keep values of elements above threshold
  storage_list_MS_calls <- list()
  
  Calls_output <- Predict_Cell_from_ModScore(SO_Trunc_Metadata_General,General_thres_vec,rejection = "unknown")
  Calls_output$Barcode <- rownames(Calls_output)
  
  if (multi_level_class){   
    
    for (k in 1:ncol(levels_dataframe)){ 
      
      # Initialize list for temporarily keeping results from subpopulation calls
      Sub_class_calls <- list()
      
      ## Subclass Identification
      # Remove any NAs from comparisons
      Sub_class_storage <- levels_dataframe[[k]][!is.na(levels_dataframe[[k]])]
      
      parent_class <- unique(gsub("(.*)-(.*)","\\1",Sub_class_storage))
      
      for (parent in parent_class){
        Sub_class <- Sub_class_storage[grepl(parent,Sub_class_storage)]
        children_class <- gsub("(.*)-(.*)","\\2",Sub_class)
        
        # Subset out cells predicted to be parents. Child cells will be screened from this population
        parents <- Calls_output$Barcode[Calls_output$Likely_CellType == parent]
        SO_Trunc_Metadata_parents <- SO.sub@meta.data[parents,] %>% select(children_class)
        
        # Stores a new density plot containing MS information of children cells within parent population
        for (children in children_class){
          
          plot_title <- paste("Density plot for",children,"Module Scores within", parent,"population", sep = " ")
          adjusted_density_plot <- ggplot(SO_Trunc_Metadata_parents, aes_string(x = children)) + geom_density() + ggtitle(plot_title) +                  geom_vline(xintercept = manual_threshold[children], linetype = "dashed", color = "red3")
          
          figures[[length(figures) + 1]] <- adjusted_density_plot
        }
        
        # Create output table without information of parent population, will be the table with updated cell calls will be reinserted
        SO_Trunc_Metadata_no_parents <- Calls_output[!Calls_output$Likely_CellType == parent,]
        non_parents <- rownames(SO_Trunc_Metadata_no_parents)
        
        gen_annot_table <- data.frame(cells = non_parents, identity = SO_Trunc_Metadata_no_parents$Likely_CellType)
        
        # Repeat Module Score Comparison and Cell Prediction with Child Subset:
        children_thres_vec <- manual_threshold[children_class]
        
        Sub_class_calls[[match(parent,parent_class)]] <- Predict_Cell_from_ModScore(SO_Trunc_Metadata_parents,children_thres_vec,rejection = parent) %>%  select(Likely_CellType)
      }
      
      # Reappend subclassification results back to general output for higher level classification
      Sub_class_calls <- do.call(rbind,Sub_class_calls)
      Sub_class_calls$Barcode <- rownames(Sub_class_calls)
      
      # Continuously update Likely_CellType column in Calls Output with new results from each level subclassification(s)
      Calls_output$Likely_CellType_updated <- Sub_class_calls$Likely_CellType[match(Calls_output$Barcode,Sub_class_calls$Barcode)]
      
      Calls_output <- Calls_output %>% mutate(Likely_CellType = case_when(
        is.na(Likely_CellType_updated) ~ Likely_CellType,
        TRUE ~ Likely_CellType_updated
      ))
      
      # Remove Likely_CellType_updated
      Calls_output$Likely_CellType_updated <- NULL
    }}
  
  ## Updating CellType(s) in metadata with subclass calls
  SO.sub@meta.data$Likely_CellType <- Calls_output$Likely_CellType[match(SO.sub@meta.data$Barcode,Calls_output$Barcode)]
  
  modscore_res <- list(Arranged_figures = do.call(arrangeGrob, c(figures)),
                       Seurat_Object = SO.sub)
  
  return(modscore_res)
 }
 