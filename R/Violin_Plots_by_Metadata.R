# This code comes from NIDAP 'Violin Plots by Metadata Column [scRNA-seq][CCBR]' code template
# Template Manager https://nidap.nih.gov/workspace/vector/templates/ri.vector.main.template.1011399a-426d-4aa1-bb3f-4430e97ea9ec
# Documentation https://nidap.nih.gov/workspace/notepad/view/ri.notepad.main.notepad.fbef08da-93ba-4301-aac7-2545631aad4f
# splitFacet helper function comes from https://stackoverflow.com/questions/30510898/split-facet-plot-into-list-of-plots/52225543

#' @title Violin Plot by Metadata
#' @description Returns a violin plot of gene expression data that is partitioned by user defined metadata column
#' @details Takes in a list of genes inputted by the user, displays gene expression information in particular slot-assay with (optional) outliers removed
#' 
#' @param seurat_object Seurat-class object
#' @param assay Name of the assay to extract gene expression data from
#' @param slot Name of the slot within the assay to extract gene expression data from
#' @param ident_of_interest Split violin plot based on a column from the seurat object metadata
#' @param groups_of_interest Include only a specific subset from ident_of_interest
#' @param genes_of_interest Genes to visualize on the violin plot
#' @param filter_outliers Filter outliers from the data (TRUE/FALSE)
#' @param scale_data Scale data from 0 to 1 (TRUE/FALSE)
#' @param log_scale_data Transform data onto a log10 scale (TRUE/FALSE)
#' @param reorder_ident Numeric data will be ordered naturally by default. Toggling this option will order the groups to match the group list if non-numeric, and will have no effect if otherwise.
#' @param rename_ident Give alternative names to ident_of_interest displayed on the graph
#' @param ylimit Y-axis limit
#' @param plot_style Choose between grid, labeled, and row
#' @param outlier_low_lim Filter lower bound outliers (default = bottom 10 percentile)
#' @param outlier_up_lim Filter upper bound outliers (default = bottom 90 percentile)
#' @param jitter_points Scatter points on the plot (TRUE/FALSE)
#' @param jitter_width Set spread of jittered points 
#' @param jitter_dot_size Set size of individual points
#' @param print_outliers Print outliers as points in your graph that may be redundant to jitter 

#' @import Seurat 
#' @import reshape2
#' @import tidyverse
#' @import cowplot
#' @import rlang
#' @import ggplot2
#'   
#' @export

#' @return violin ggplot2 object

ViolinPlot <- function(so, 
                       assay = "SCT", 
                       slot = "scale.data", 
                       ident_of_interest, 
                       groups_of_interest, 
                       genes_of_interest,
                       filter_outliers = FALSE, 
                       scale_data = TRUE, 
                       log_scale_data = FALSE, 
                       reorder_ident = TRUE,
                       rename_ident = "", 
                       ylimit = 0, 
                       plot_style = "grid", 
                       outlier_low_lim = 0.1, 
                       outlier_up_lim = 0.9, 
                       jitter_points = FALSE,
                       jitter_width = 0.05, 
                       jitter_dot_size = 0.5, 
                       print_outliers = TRUE){

  ## ---------------##
  ## Error Messages ##
  ## ---------------##
  
  gene_filter <- genes_of_interest %in% rownames(GetAssayData(object = so, slot = slot, assay = assay))
  missing_genes <- genes_of_interest[!gene_filter]
  genes_of_interest <- genes_of_interest[gene_filter]
  
  if(length(missing_genes) > 0){
    print(paste("The following genes are missing from the dataset:", missing_genes, sep = " "))
  }
  
  if(length(genes_of_interest) == 0){
    stop("No query genes were found in the dataset.")
  }
  
  if(!ident_of_interest %in% colnames(so@meta.data)){
    colnames(so@meta.data) <- gsub("\\.","_",colnames(so@meta.data))
    if(!ident_of_interest %in% colnames(so@meta.data)){
      stop("Unable to find ident of interest in metadata.")
    }
  }
  
  group_filter <- groups_of_interest %in% so@meta.data[[ident_of_interest]]
  groups_of_interest <- groups_of_interest[group_filter]
  missing_groups <- groups_of_interest[!group_filter]
  if(length(missing_groups) > 0){
    cat("The following groups are missing from the selected ident:\n")
    print(missing_groups)
  }
  
  if(length(groups_of_interest) == 0){
    stop("No groups were found in the selected ident.")
  }
  
  if(rename_ident %in% c("Gene","Expression","scaled")){
    stop("New ident name cannot be one of Gene, Expression, or scaled.")
  }
  
  ## ---------------- ##
  ## Helper Functions ##
  ## ---------------- ##
  
  splitFacet <- function(x){
    facet_vars <- names(x$facet$params$facets)         # 1
    x$facet    <- ggplot2::ggplot()$facet              # 2
    datasets   <- split(x$data, x$data[facet_vars])    # 3
    new_plots  <- lapply(datasets,function(new_data) { # 4
      x$data <- new_data
      x})
  }
  
  # from rapportools
  is.empty <- function(x, trim = TRUE, ...) {
    if (length(x) <= 1) {
      if (is.null(x))
        return (TRUE)
      if (length(x) == 0)
        return (TRUE)
      if (is.na(x) || is.nan(x))
        return (TRUE)
      if (is.character(x) && nchar(ifelse(trim, trim.space(x), x)) == 0)
        return (TRUE)
      if (is.logical(x) && !isTRUE(x))
        return (TRUE)
      if (is.numeric(x) && x == 0)
        return (TRUE)
      return (FALSE)
    } else
      sapply(x, is.empty, trim = trim, ...)
  }
  
  # from rapportools
  trim.space <- function(x, what = c('both', 'leading', 'trailing', 'none'), space.regex = '[:space:]', ...){
    if (missing(x))
      stop('nothing to trim spaces to =(')
    re <- switch(match.arg(what),
                 both     = sprintf('^[%s]+|[%s]+$', space.regex, space.regex),
                 leading  = sprintf('^[%s]+', space.regex),
                 trailing = sprintf('[%s]+$', space.regex),
                 none     = {
                   return (x)
                 })
    vgsub(re, '', x, ...)
  }
  
  # from rapportools
  vgsub <- function(pattern, replacement, x, ...){
    for(i in 1:length(pattern))
      x <- gsub(pattern[i], replacement[i], x, ...)
    x
  }
  
  remove_outliers_func <- function(x, na.rm = TRUE){
    qnt <- quantile(x, probs=c(outlier_low_lim,outlier_up_lim), na.rm = na.rm)
    H <- 1.5*IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
  }
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  # deal with limits
  if(ylimit == 0){
    ylimit <- NULL
  }
  
  #set idents
  if("active.ident" %in% slotNames(so)){
    new_idents = as.factor(so@meta.data[[ident_of_interest]])
    names(new_idents) = names(so@active.ident)
    so@active.ident <- as.factor(vector())
    so@active.ident <- new_idents
    so.sub = subset(so, ident = groups_of_interest)
  } else {
    new_idents = as.factor(so@meta.data[[ident_of_interest]])
    names(sample_name)=so@meta.data[["Barcode"]]
    so@active.ident <- as.factor(vector())
    so@active.ident <- new_idents
    so.sub = subset(so, ident = groups_of_interest)
  }
  
  DefaultAssay(object = so) <- assay
  data <- FetchData(object = so.sub, vars = genes_of_interest, slot = slot)
  
  data[[ident_of_interest]] <- so.sub@meta.data[row.names(data),ident_of_interest]
  
  df.melt <- reshape2::melt(data)
  
  if(!is.empty(rename_ident)){
    ident_of_interest <- rename_ident
  }
  colnames(df.melt) <- c(ident_of_interest,"Gene","Expression")
  
  #check to see if ident of interest looks numeric
  if(suppressWarnings(all(!is.na(as.numeric(as.character(df.melt[[ident_of_interest]])))))){
    ident.values <- strtoi(df.melt[[ident_of_interest]])
    ident.levels <- unique(ident.values)[order(unique(ident.values))]
    df.melt[[ident_of_interest]] <- factor(ident.values, levels = ident.levels)
  }else if(reorder_ident){
    # if non-numeric, place in order of groups of interests
    df.melt[[ident_of_interest]] <- factor(df.melt[[ident_of_interest]], levels = groups_of_interest)        
  }
  
  # Filter outliers
  if(filter_outliers){
    for(gene in genes_of_interest){
      for(group in groups_of_interest){
        current.ind <- which(df.melt[["Gene"]] == gene & df.melt[[ident_of_interest]] == group)
        df.melt[current.ind,"Expression"] <- remove_outliers_func(df.melt[current.ind,"Expression", drop = TRUE])
      }
    }
  }
  
  # Scale data to y limit
  if(scale_data){
    expression_data = "scaled"
    axis.title.y = "Expression (scaled)"
    ylimit <- ylimit %||% 1
    df.melt <- df.melt %>% group_by(Gene) %>% mutate(scaled = scales::rescale(Expression, to=c(0,ylimit)))
  }else{
    expression_data <- axis.title.y <- "Expression"
  }
  
  g <- ggplot(df.melt, aes_string(x=ident_of_interest, y=expression_data)) +
    geom_violin(aes_string(fill = ident_of_interest), scale="width", trim = FALSE, show.legend = FALSE) + 
    #geom_jitter(height = 0, width = 0.05, size=0.1) +
    theme_classic() + 
    #scale_fill_manual(values=cols) + 
    labs(y=axis.title.y) +
    theme(strip.text.y = element_text( 
      color="blue", face="bold.italic", angle = -90))
  
  
  if(!is.null(ylimit)){
    g <- g + ylim(0,ylimit)
  }
  
  if(jitter_points){
    g <- g + geom_jitter(height = 0, width = jitter_width, size=jitter_dot_size)
  }
  
  if(log_scale_data){
    g <- g + scale_y_log10()
  }
  
  # Plot after jitter if wanted
  g <- g + geom_boxplot(width=0.1, fill="white", outlier.shape = ifelse(print_outliers,19,NA))
  
  # Plot styles
  ncol = ceiling(length(unique(df.melt$Gene))^0.5)
  nrow = ceiling(length(unique(df.melt$Gene)) / ncol)
  if(plot_style == "rows"){
    g <- g + facet_grid(rows = vars(Gene))
  }else{
    g <- g + facet_wrap(~Gene, nrow = nrow, ncol = ncol)
    if(plot_style == "labeled"){
      plots <- splitFacet(g)
      plots <- lapply(seq_along(plots), function(i) plots[[i]] + ggtitle(genes_of_interest[i]) + theme(plot.title = element_text(hjust = 0.5)) )
      g <- plot_grid(plotlist = plots, nrow = nrow, ncol = ncol, labels = LETTERS[seq( from = 1, to = length(plots) )])
    }
  }
  
  return(g)
}
