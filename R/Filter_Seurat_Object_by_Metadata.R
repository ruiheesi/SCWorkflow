# [scRNA-Seq][CCBR] Filter Seurat Object by Metadata (ec3f23f9-bcba-4f3a-8a08-8ba611fbb6c7): v124

#' @title Filter Seurat Object by Metadata
#' @description Filter and subset your Seurat object based on some metadata column.
#' @details This is a downstream template that should be loaded after Step 5 of the pipeline (SingleR Annotations on Seurat Object).
#' 
#' @param object A dataset containing your SingleR annotated/merged seurat object.
#' @param samples.to.include select which samples to include.
#' @param keep.or.remove TRUE to filter down to selected values, FALSE to filter out selected values.
#' @param greater.less.than Decide if you want to keep cells above or below the threshold.
#' @param seed Set seed for colors.
#' @param sample.name Sample Name Column.
#' @param cut.off The cut-off you want to use for your greater than/less than filter.
#' @param category.to.filter What kind of metadata you want to subset by. This should be one column in your Metadata table. 
#' @param legend.position Select "none" if legend takes up too much space in plot.
#' @param values.to.filter One or more values where you want to filter.
#' @param reduction What kind of clustering visualization you would like to use for your summary plot (umap, tsne, pca).
#' @param plot.as.interactive.plot TRUE for interactive, FALSE for static.
#' @param legend.symbol.size Default is 2.
#' @param colors User-selected colors from palette of 62 unique colors from ColorBrewer.
#' @param dot.size Size of dots on TSNE/UMAP projection plot.
#' @param number.of.legend.columns Default is 1. If legend is too long, provide more legend columns
#' @param dot.size.highlighted.cells Dot size for cells in the after-filter plot which have been highlighted.
#' @param use.cite.seq.data TRUE if you would like to plot Antibody clusters from CITEseq instead of scRNA.

#' 
#' @import Seurat 
#' @import ggplot2 
#' @import gridExtra 
#' @import grid 
#' @import gridBase 
#' @import cowplot 
#' @import RColorBrewer 
#' @import colorspace 
#' @import tidyverse 
#' 
#' 
#'   
#' @export
#' 
#' @return a subset Seurat object




filterSeuratObjectByMetadata <- function(
  object,
  samples.to.include,
  keep.or.remove = TRUE,
  greater.less.than = "greater than",
  seed = 10,
  sample.name,
  cut.off = 0.5,
  category.to.filter,
  legend.position = "top",
  values.to.filter,
  reduction = "umap",
  plot.as.interactive.plot = FALSE,
  legend.symbol.size = 2,
  colors = c("aquamarine3","salmon1","lightskyblue3","plum3","darkolivegreen3","goldenrod1","burlywood2","gray70","firebrick2","steelblue","palegreen4","orchid4","darkorange1","yellow","sienna","palevioletred1","gray60","cyan4","darkorange3","mediumpurple3","violetred2","olivedrab","darkgoldenrod2","darkgoldenrod","gray40","palegreen3","thistle3","khaki1","deeppink2","chocolate3","paleturquoise3","wheat1","lightsteelblue","salmon","sandybrown","darkolivegreen2","thistle2","gray85","orchid3","darkseagreen1","lightgoldenrod1","lightskyblue2","dodgerblue3","darkseagreen3","forestgreen","lightpink2","mediumpurple4","lightpink1","thistle","navajowhite","lemonchiffon","bisque2","mistyrose","gray95","lightcyan3","peachpuff2","lightsteelblue2","lightyellow2","moccasin","gray80","antiquewhite2","lightgrey"),
  dot.size = 0.1,
  number.of.legend.columns = 1,
  dot.size.highlighted.cells = 0.5,
  use.cite.seq.data = FALSE
  ){

  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  #image: image.type
  

  
  samples = eval(parse(text=gsub('\\[\\]','c()',samples.to.include)))
  
  if (length(samples) == 0) {
##    samples = unique(object@meta.data$sample.name)
    samples = unique(object@meta.data[[sample.name[1]]])
  }
  
  ## Replace dots in metadata column names with underscores.
  colnames(object@meta.data) = gsub("\\.", "_", colnames(object@meta.data))
  new_sample_name <- gsub("\\.", "_", sample.name[1])
  
  ## If you have protien data, then ...
  if (use.cite.seq.data) {
    reduction =paste("protein_", reduction,sep='')
  }
  

  ## Original color-picking code.
  n <- 2e3
  set.seed(seed)
  ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
  ourColorSpace <- as(ourColorSpace, "LAB")
  distinctColorPalette <-function(k=1,seed) {
    currentColorSpace <- ourColorSpace@coords
    # Set iter.max to 20 to avoid convergence warnings.
    set.seed(seed)
    km <- kmeans(currentColorSpace, k, iter.max=20)
    colors <- unname(hex(LAB(km$centers)))
    return(colors)
  }
  
  ## User-selected metadata column is used to set idents.
  Filter.orig =object@meta.data[category.to.filter[1]]
  colname <- category.to.filter[1]
  
  ident_of_interest = as.factor(object@meta.data[[colname]])
  names(ident_of_interest)=names(object@active.ident)
  object@active.ident <- as.factor(vector())
  object@active.ident <- ident_of_interest
  
  ## Get colors from user parameter and add more if the default list is too short.
  if(class(object@meta.data[[category.to.filter[1]]]) != "numeric"){
    q = length(levels(as.factor(Filter.orig[[colname]])))
    if(length(colors) < q) {
      r = q - length(colors)
      more_cols = distinctColorPalette(r,10)
      colors <- c(colors, more_cols)
    }
    names(colors) <- levels(as.factor(Filter.orig[[colname]]))
    
    ## Keep or remove cells based on user input values.
    if (keep.or.remove) {
      subsetValue <- values.to.filter
      metaCol <- unique(object@meta.data[[category.to.filter[1]]])
      print("Missing values:")
      print(setdiff(subsetValue,metaCol))
      subsetValue <- intersect(metaCol,subsetValue)
    } else {
      metaCol <- unique(object@meta.data[[colname]])
      valsToRemove <- values.to.filter
      subsetValue <- setdiff(metaCol, valsToRemove)
    }
    
    ## Subset Seurat object.
    #SO.sub <-SubsetData(object, ident.use=subsetValue)
    SO.sub <- subset(object, idents = subsetValue)
    
    ## Log output of tables of cell types by samples before and after filtes.
    print("Breakdown of filtered data:")
##    print(table(object@meta.data[[category.to.filter[1]]],object@meta.data$sample_name))
    print(table(object@meta.data[[category.to.filter[1]]],object@meta.data[[new_sample_name]]))
    
    cat("\n")
    print("After Filtering:")
##    print(table(SO.sub@meta.data[[category.to.filter[1]]],SO.sub@meta.data$sample_name))
    print(table(SO.sub@meta.data[[category.to.filter[1]]],SO.sub@meta.data[[new_sample_name]]))
    
    ## Set filter for the subsetted SO.
    SO.sub@meta.data[[colname]] <- as.factor(as.character(SO.sub@meta.data[[colname]])) #Relevel Factors
    
    Filter.sub = SO.sub@meta.data[[colname]]
    
    ## More color stuff.
    #Set colors for unfiltered and filtered data by sample name:
    n = length(levels(as.factor(Filter.sub)))
    idx = vector("list", n)
    names(idx) <- levels(as.factor(Filter.sub))
    for (i in 1:n) {
      id = Filter.orig %in% levels(as.factor(Filter.sub))[i]
      idx[[i]] <- rownames(object@meta.data)[id]
    }
    cols2 <- colors[levels(as.factor(Filter.sub))]
    
    ## Make before and after plots.
    title <- paste0("filtered by ", category.to.filter[1], " and split by ", category.to.filter[2])
    p1 = DimPlot(object, reduction=reduction, 
                 group.by=colname,
                 pt.size=dot.size) + 
      theme_classic() + 
      scale_color_manual(values=colors) + 
      theme(legend.position=legend.position) +
      guides(colour=guide_legend(ncol=number.of.legend.columns,override.aes = list(size = legend.symbol.size))) +
      ggtitle(colname)
    p2 = DimPlot(object, reduction=reduction, 
                 cells.highlight = idx, 
                 cols.highlight= rev(cols2[1:n]), 
                 sizes.highlight = dot.size.highlighted.cells) + 
      theme_classic() + 
      theme(legend.position=legend.position)+
      guides(colour=guide_legend(ncol=number.of.legend.columns,reverse=TRUE,override.aes = list(size = legend.symbol.size))) +
      ggtitle(title)
    
    ## Else, filter on numeric data with a user defined threshold and direction.
  } else {
    filterDirection <-greater.less.than
    metaCol <- unique(object@meta.data[[category.to.filter]])
    value <- cut.off
    if (filterDirection =="greater than") {
      SO.sub <- subset(object, subset = category.to.filter > cut.off)
    } else {
      SO.sub <- subset(object, subset = category.to.filter < cut.off)
    }
    
    
    drawtsne <- function(SO,reduction,scalecol,colgrad){
      
      SO.clus <- SO@meta.data[[m]]
      
      p1 <- DimPlot(SO, reduction = reduction, group.by = "ident")
      class(p1$data$ident) <- "numeric"
      
      if(reduction=="tsne"){
        clusmat=data.frame(umap1=p1$data$tSNE_1,umap2=p1$data$tSNE_2, clusid=as.numeric(SO@meta.data[[m]]))
      } else if(reduction=="umap"){
        clusmat=data.frame(umap1=p1$data$UMAP_1,umap2=p1$data$UMAP_2, clusid=as.numeric(SO@meta.data[[m]]))
      } else if (reduction=="pca"){
        clusmat=data.frame(umap1=p1$data$PC_1,umap2=p1$data$PC_2, clusid=as.numeric(SO@meta.data[[m]]))
      } else if (reduction=="protein_tsne"){
        clusmat=data.frame(umap1=p1$data$protein_tsne_1,umap2=p1$data$protein_tsne_2, clusid=as.numeric(SO@meta.data[[m]]))
      } else if (reduction=="protein_umap"){
        clusmat=data.frame(umap1=p1$data$protein_umap_1,umap2=p1$data$protein_umap_2, clusid=as.numeric(SO@meta.data[[m]]))
      } else {
        clusmat=data.frame(umap1=p1$data$protein_pca_1,umap2=p1$data$protein_pca_2, clusid=as.numeric(SO@meta.data[[m]]))
      }
      
      clusmat %>% group_by(clusid) %>% summarise(umap1.mean=mean(umap1), umap2.mean=mean(umap2)) -> umap.pos
      title=as.character(m)
      clusmat %>% dplyr::arrange(clusid) -> clusmat
      
      p2 <- ggplot(clusmat, aes(x=umap1, y=umap2)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        geom_point(aes(colour=clusid),alpha=0.5,shape = 20,size=dot.size) +
        #scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = scales::rescale(c(min, midpt2,midpt,midpt3, max), limits = c(0, 1))) +
        #scale_color_gradientn(colours = c("blue4", "lightgrey", "red"), values = c(min, midpt2,midpt,midpt3, max)) + 
        #scale_color_gradientn(colors=brewer.pal(n = 5, name = "RdBu"), values = scales::rescale(c(min, midpt2,midpt,midpt3, max))) +
        
        scale_color_gradientn(colors=brewer.pal(n = 5, name = colgrad), values = scalecol) +
        
        guides(colour = guide_legend(override.aes = list(size=5, alpha = 1))) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
        ggtitle(title) +
        xlab("umap-1") + ylab("umap-2")
      return(p2)
    }
    
    m = category.to.filter
    clusid = object@meta.data[[m]]
    maxclus = max(clusid)
    clusmid = 0.01/maxclus        
    min = min(clusid)
    midpt1 = 0.99*value
    midpt = value
    midpt2 = 1.01*value
    max = max(clusid)
    colpoints <- c(min,midpt1,midpt,midpt2,max)
    colpoints <- scales::rescale(colpoints,c(0,1))
    
    p1 <- drawtsne(object,reduction,colpoints,"RdBu")
    
    clusid = scales::rescale(SO.sub@meta.data[[m]], to=c(0,1))
    clus.quant=quantile(clusid[clusid>0],probs=c(0,.25,.5,.75,1))
    min = clus.quant[1]
    midpt = clus.quant[3]
    midpt3 = clus.quant[2]
    midpt4 = clus.quant[4]
    max = clus.quant[5]  
    colpoints2 <- c(min,midpt3,midpt,midpt4,max)
    
    p2 <- drawtsne(SO.sub,reduction,colpoints2,"Blues")
    
  }
  
#  ## If interactive plot requested, then ...
#  if (plot.as.interactive.plot == TRUE) {
#    gp1 <- ggplotly(p1)
#    gp2 <- ggplotly(p2)
#    p <- subplot(gp1, gp2, nrows=2)
#    print(p)
#  } else {
#    ## Else, print non-interactive plot.
#    print(plot_grid(p1,p2,nrow=1))
#  }
  
#  ## Return the subsetted Seurat object.
#  output <- new.output()
#  output_fs <- output$fileSystem()
#  saveRDS(SO.sub, output_fs$get_path("seurat_object.rds", 'w'))
  
  #return(list(SO.sub, p1, p2))
  
  result.list <- list("object" = SO.sub,"plot1" = p1,"plot2" = p2)
  return(result.list)
  }
