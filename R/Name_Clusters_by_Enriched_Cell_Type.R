#' Updating cluster annotation in the meta.data slot of Seurat-class object
#' 
#' Maps cluster IDs to cluster names and adds a new metadata column called *Clusternames*. Seurat object is returned if all provided cluster IDs and cluster number from the Seurat Object meta.data match, otherwise a data.frame with counts of *Likely_CellTypes* within each cluster is returned from the Seurat object.
#' 
#' 
#' @param SO Seurat-class object with cluster a cluster IDs column and "Likely_CellTypes" column present
#' @param cluster.identities.table a data.frame with 2 columns - one with Cluster IDs (numeric) and the other with Cluster names
#' @param cluster.column.from.SO name of the cluster ID column present in the meta.data slot in the Seurat Object
#' @param cluster.names column containing cluster labels
#' @param cluster.numbers column containing cluster ID's (numeric)
#' 
#' @importFrom dplyr pull
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_point aes theme_classic ylim scale_y_reverse theme ggtitle
#' @importFrom plotly ggplotly
#' @importFrom Seurat AddMetaData
#' @importFrom tibble deframe tibble
#' 
#' @export
#' @return Returns Seurat-class object with updated meta.data slot containing custom cluster annotation or extracts a data.frame from the Seurat meta.data slot with already annotated cell types

NameClusters <-
  function(SO = seurat.object,
           cluster.identities.table,
           cluster.column.from.SO,
           cluster.names,
           cluster.numbers)
  {
    # set metadata
    df <- cluster.identities.table
    metadata.df <- SO@meta.data
    colnames(metadata.df) <- gsub("\\.", "_", colnames(metadata.df))
    colval <- metadata.df[, cluster.column.from.SO]
    seurat.cluster <-
      colnames(SO@meta.data)[which(colnames(metadata.df) == cluster.column.from.SO)]
    cluster.num <-
      as.data.frame.matrix(table(colval, SO@meta.data$Likely_CellType))
    clusnum.df <- melt(as.matrix(cluster.num))
    sums <- rowSums(cluster.num)
    cluster.perc <- (cluster.num / sums) * 100
    
    # draw plot with Likely cell type numbers per cluster
    clus.df <- melt(as.matrix(cluster.perc))
    clus.df$num <- clusnum.df$value
    colnames(clus.df) <- c("cluster", "celltype", "percent", "number")
    minclus <- min(clus.df$cluster)
    maxclus <- max(clus.df$cluster)
    numclus <- length(unique(clus.df$cluster))
    
    # log
    ## clusters
    clus.df %>% pull(cluster) %>% unique() %>% as.character() -> clus
    cat("clusters in Seurat Object:\n")
    cat(clus, sep = "\n")
    ## cell types
    clus.df %>% pull(celltype) %>% unique() %>% as.character() -> celltypes
    cat("\ncelltypes in Seurat Object:\n")
    cat(celltypes, sep = "\n")
    ## cluster identities
    cat("\nclusters in cluster identities table:\n")
    cat(df[, cluster.numbers], sep = "\n")
    
    
    # do plot
    suppressMessages(g <-
      ggplot(clus.df,
             aes(
               x = celltype,
               y = cluster,
               size = percent,
               color = celltype,
               label = number
             )) +
      theme_classic() +
      geom_point(alpha = 0.5) +
      ylim(c(minclus, maxclus)) +
      scale_y_reverse(n.breaks = numclus) +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) +
      ggtitle("Percentage of Cell Type within Clusters")
    )
    
    g <- ggplotly(g)
    print(g)
    
    # If cluster numbers on input table match the cluster numbers selected:
    
    if (nrow(df) == length(unique(colval))) {
      if (all(as.numeric(df[, cluster.numbers]) == sort(unique(colval)))) {
        cat('\nAdding "Clusternames" column to metadata in the Seurat Object\n\n')
        
        output <-
          AddMetaData(SO, metadata = deframe(tibble(df[, cluster.numbers], df[, cluster.names]))[as.character(colval)], col.name = "Clusternames")
        print(table(output$Clusternames, output@meta.data[, seurat.cluster]))
        print(table(output$Clusternames, output$Likely_CellType))
        
      } else {
        output <- cluster.num
        cat(
          sprintf(
            "\n\nNo update in the Seurat Object - Cluster numbers in the metadata table are not identical with the Seurat Object's %s.\n\n",
            seurat.cluster
          )
        )
      }
      
    } else {
      output <- cluster.num
      cat(
        sprintf(
          "\n\nNo update in the Seurat Object - The number of clusters in the metadata table is not the same as in the Seurat Object's %s column.\n\n",
          seurat.cluster
        )
      )
    }
    
    cat("\nReturning dataset:\n\n")
    print(output)
    invisible(list(output = output, plot = g))
    
  }
