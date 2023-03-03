#' Updating cluster annotation in the meta.data slot of Seurat-class object
#' 
#' @description Maps cluster IDs to cluster names and adds a new metadata column 
#' called *clusternames*. 
#' 
#' @details Seurat object is returned if all provided cluster IDs 
#' and cluster number from the Seurat Object meta.data match within each cluster 
#' is returned from the Seurat object.
#' 
#' @param object Seurat-class object with cluster a cluster IDs column and "Likely_CellTypes" column present
#' @param cluster.numbers Vector containing cluster ID's (numeric)
#' @param cluster.names vector containing cluster labels
#' @param cluster.column column name containing cluster ID in the meta.data slot
#'  in the Seurat Object
#' @param labels.column column name containing labels (usually cell type) in the
#'  meta.data slot in the Seurat Object
#' @param order.clusters.by vector containing order of clusters in graph 
#' (default is NULL)
#' @param order.celltypes.by vector containing order of celltypes in graph 
#' (default is NULL)
#' @param interactive if TRUE, draw plotly plot (default is FALSE)
#' 
#' @importFrom dplyr pull
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_point aes theme_classic ylim scale_y_reverse 
#' theme ggtitle
#' @importFrom plotly ggplotly
#' @importFrom Seurat AddMetaData
#' @importFrom tibble deframe tibble
#' 
#' @export
#' @return Returns Seurat-class object with updated meta.data slot containing 
#' custom cluster annotation 

nameClusters <-function(object,
                        cluster.numbers,
                        cluster.names,
                        cluster.column,
                        labels.column,
                        order.clusters.by = NULL,
                        order.celltypes.by = NULL,
                        interactive = FALSE)
  {
    

    cluster.numbers <- factor(cluster.numbers,
                              str_sort(cluster.numbers,numeric = TRUE))
    levels(cluster.numbers) <- str_sort(cluster.numbers, numeric=TRUE)    

    # get metadata from SO and get cluster column
    metadata.df <- object@meta.data
    colval <- metadata.df[[cluster.column]]
    names(cluster.names) <- as.character(cluster.numbers)
    
    # If cluster numbers on input table match the cluster numbers selected:
    
    if (levels(metadata.df[[cluster.column]]) %in% levels(colval)) {
      object <- AddMetaData(object, 
                  metadata = deframe(tibble(metadata.df[[cluster.column]], 
                  cluster.names))[as.character(colval)], col.name = "Clusternames")
    } else{
      stop("Cluster ID's have to match metadata column. Please check entry in input table.")
    }
    
    cluster_num <- as.data.frame.matrix(table(SO@meta.data$Clusternames, SO@meta.data[[Celltype_Column_from_Seurat_Object_Metadata]]))
    clusnum.df <- melt(as.matrix(cluster_num))
    sums <- rowSums(cluster_num)
    cluster_perc <- (cluster_num / sums) * 100
    
    # draw plot with Likely cell type numbers per cluster
    clus.df <- melt(as.matrix(cluster_perc))
    clus.df$num <- clusnum.df$value
    colnames(clus.df) <- c("cluster", "celltype", "percent", "number")
    head(clus.df)
    clus.df$cluster <- factor(clus.df$cluster)
    
    if(length(Order_Clusters_by) == 0){
      levels(clus.df$cluster) <- str_sort(levels(clus.df$cluster),numeric=TRUE)
    } else {
      levels(clus.df$cluster) <- Order_Clusters_by
    }
    
    if(length(Order_Celltypes_by) == 0){
      levels(clus.df$celltype) <- str_sort(levels(clus.df$celltype),numeric=TRUE)
    } else {
      levels(clus.df$celltype) <- Order_Celltypes_by
    }
    
    clus.df %>% arrange(cluster,celltype) -> clus.df
    
    ## print clusters in log:
    clus.df %>% pull(cluster) %>% unique() %>% as.character() -> clus
    cat("clusters:\n")
    cat(clus, sep = "\n")
    cat("\n\n")
    
    ## print cell types in log
    clus.df %>% pull(celltype) %>% unique() %>% as.character() -> celltypes
    cat("celltypes:\n")
    cat(celltypes, sep = "\n")
    
    # if (length(cluster.numbers) == length(unique(colval))) {
    #   if (all(as.numeric(cluster.numbers) == sort(unique(colval)))) {
    #     cat('\nAdding "Clusternames" column to metadata in the Seurat Object
    #         \n\n')
    #     clusternames <- cluster.names[object@meta.data[[cluster.column]]]
    #     object <- AddMetaData(object, 
    #                           metadata = clusternames,
    #                           col.name="clusternames")
    #     metatable <- object@meta.data
    #     contable <- as.data.frame.matrix(table(metatable$clusternames, 
    #                                            metatable[[labels.column]]))
    #     table <- tableGrob(contable)
    #   } else {
    #     table <- cluster.num
    #     cat(
    #       sprintf(
    #         "\n\nNo update in the Seurat Object - Cluster numbers in the 
    #         metadata table are not identical with the Seurat Object's %s.\n\n",
    #         cluster.column
    #       )
    #     )
    #   }
    # } else {
    #   table <- cluster.num
    #   cat(
    #     sprintf(
    #       "\n\nNo update in the Seurat Object - The number of clusters in the 
    #       metadata table is not the same as in the Seurat Object's %s column.
    #       \n\n",
    #       cluster.column
    #     )
    #   )
    # }
    # 
    # cluster.num <-
    #   as.data.frame.matrix(table(colval, object@meta.data[[labels.column]]))
    # clusnum.df <- melt(as.matrix(cluster.num))
    # sums <- rowSums(cluster.num)
    # zero <- sums == 0
    # if(sum(zero) > 0){
    #   empty <- paste(names(sums[zero]),collapse = ",")
    #   sums <- sums[!zero]
    #   cluster.num <- cluster.num[!zero,]
    #   clusnum.df %>% filter(Var1 %in% names(zero[zero==FALSE])) -> clusnum.df
    #   warning(sprintf("Some clusters had no detected cell types: %s", empty))
    # }
    # cluster.perc <- (cluster.num / sums) * 100
    # 
    # # draw plot with Likely cell type numbers per cluster
    # clus.df <- melt(as.matrix(cluster.perc))
    # clus.df$num <- clusnum.df$value
    # colnames(clus.df) <- c("cluster", "celltype", "percent", "number")
    # clus.df$cluster <- cluster.names[as.character(clus.df$cluster)]
    # 
    # if(length(order.celltypes.by) > 0){
    #   #Check for factors not present in dataset
    #   add <- order.celltypes.by[!order.celltypes.by %in% levels(clus.df$celltype)]
    #   if(length(add) > 0){
    #     add <- paste(add,collapse = ", ")
    #     warning(sprintf("Some factors are not in data: %s",add))
    #   }
    #   clus.dat <- unique(clus.df$celltype)
    #   minus <- clus.dat[!clus.dat %in% order.celltypes.by]
    #   if(length(minus)> 0){
    #     minus <- paste(minus,collapse = ", ")
    #     warning(sprintf("Some factors were not included in the list: %s",minus))
    #   }
    #   clus.df <- clus.df %>% filter(celltype %in% order.celltypes.by) 
    #   clus.df$celltype <- factor(clus.df$celltype,levels = order.celltypes.by)
    # } 
    # if(length(order.clusters.by) > 0){
    #   add2 <- order.clusters.by[!order.clusters.by %in% levels(clus.df$celltype)]
    #   if(length(add2) > 0){
    #     add2 <- paste(add2,collapse = ", ")
    #     warning(sprintf("Some factors are not in data: %s",add2))
    #   }
    #   clus.dat2 <- unique(clus.df$cluster)
    #   minus2 <- clus.dat2[!clus.dat2 %in% order.clusters.by]
    #   if(length(minus2) > 0){
    #     minus2 <- paste(minus2,collapse = ", ")
    #     warning(sprintf("Some factors were not included in the list: %s",minus2))
    #   }
    #   clus.df <- clus.df %>% filter(cluster %in% order.clusters.by) 
    #   clus.df$cluster <- factor(clus.df$cluster,levels = order.clusters.by)
    # } else {
    #   clus.df$cluster <- factor(clus.df$cluster,
    #                             levels = str_sort(cluster.names,numeric = TRUE))
    # }
    # if (length(cluster.numbers) == length(unique(colval))) {
    #   if (all(as.numeric(cluster.numbers) == sort(unique(colval)))) {
    #     cat('\nAdding "Clusternames" column to metadata in the Seurat Object
    #         \n\n')
    #     clusternames <- cluster.names[object@meta.data[[cluster.column]]]
    #     object <- AddMetaData(object, 
    #                           metadata = clusternames,
    #                           col.name="clusternames")
    #     metatable <- object@meta.data
    #     contable <- as.data.frame.matrix(table(metatable$clusternames, 
    #                                            metatable[[labels.column]]))
    #     table <- tableGrob(contable)
    #   } else {
    #     table <- cluster.num
    #     cat(
    #       sprintf(
    #         "\n\nNo update in the Seurat Object - Cluster numbers in the 
    #         metadata table are not identical with the Seurat Object's %s.\n\n",
    #         cluster.column
    #       )
    #     )
    #   }
    # } else {
    #   table <- cluster.num
    #   cat(
    #     sprintf(
    #       "\n\nNo update in the Seurat Object - The number of clusters in the 
    #       metadata table is not the same as in the Seurat Object's %s column.
    #       \n\n",
    #       cluster.column
    #     )
    #   )
    # }
    # 
    # cluster.num <-
    #   as.data.frame.matrix(table(colval, object@meta.data[[labels.column]]))
    # clusnum.df <- melt(as.matrix(cluster.num))
    # sums <- rowSums(cluster.num)
    # zero <- sums == 0
    # if(sum(zero) > 0){
    #   empty <- paste(names(sums[zero]),collapse = ",")
    #   sums <- sums[!zero]
    #   cluster.num <- cluster.num[!zero,]
    #   clusnum.df %>% filter(Var1 %in% names(zero[zero==FALSE])) -> clusnum.df
    #   warning(sprintf("Some clusters had no detected cell types: %s", empty))
    # }
    # cluster.perc <- (cluster.num / sums) * 100
    # 
    # # draw plot with Likely cell type numbers per cluster
    # clus.df <- melt(as.matrix(cluster.perc))
    # clus.df$num <- clusnum.df$value
    # colnames(clus.df) <- c("cluster", "celltype", "percent", "number")
    # clus.df$cluster <- cluster.names[as.character(clus.df$cluster)]
    # 
    # if(length(order.celltypes.by) > 0){
    #   #Check for factors not present in dataset
    #   add <- order.celltypes.by[!order.celltypes.by %in% levels(clus.df$celltype)]
    #   if(length(add) > 0){
    #     add <- paste(add,collapse = ", ")
    #     warning(sprintf("Some factors are not in data: %s",add))
    #   }
    #   clus.dat <- unique(clus.df$celltype)
    #   minus <- clus.dat[!clus.dat %in% order.celltypes.by]
    #   if(length(minus)> 0){
    #     minus <- paste(minus,collapse = ", ")
    #     warning(sprintf("Some factors were not included in the list: %s",minus))
    #   }
    #   clus.df <- clus.df %>% filter(celltype %in% order.celltypes.by) 
    #   clus.df$celltype <- factor(clus.df$celltype,levels = order.celltypes.by)
    # } 
    # if(length(order.clusters.by) > 0){
    #   add2 <- order.clusters.by[!order.clusters.by %in% levels(clus.df$celltype)]
    #   if(length(add2) > 0){
    #     add2 <- paste(add2,collapse = ", ")
    #     warning(sprintf("Some factors are not in data: %s",add2))
    #   }
    #   clus.dat2 <- unique(clus.df$cluster)
    #   minus2 <- clus.dat2[!clus.dat2 %in% order.clusters.by]
    #   if(length(minus2) > 0){
    #     minus2 <- paste(minus2,collapse = ", ")
    #     warning(sprintf("Some factors were not included in the list: %s",minus2))
    #   }
    #   clus.df <- clus.df %>% filter(cluster %in% order.clusters.by) 
    #   clus.df$cluster <- factor(clus.df$cluster,levels = order.clusters.by)
    # } else {
    #   clus.df$cluster <- factor(clus.df$cluster,
    #                             levels = str_sort(cluster.names,numeric = TRUE))
    # }
    
    # do plot (suppressMessages for ggplot2 scale replacemnt)
      g <- ggplot(clus.df, aes(
              x = celltype,
              y = cluster,
              size = percent,
              color = celltype,
              label = number)) +
              theme_classic() +
              geom_point(alpha = 0.5) +
              theme(axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 1)) +
              ggtitle("Percentage of Cell Type within Clusters")
    
    if(interactive == TRUE){
      g <- ggplotly(g)
    }
    
    invisible(list(object = object, table = table, plot = g))
    
  }
