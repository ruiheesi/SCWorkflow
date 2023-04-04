#' @title Filter & QC Samples 
#' @description Creates Seurat Objects, Filters and Plots QC (before and after
#'  filtering) on samples. Reads h5 files (each representing one sample),
#'   and performs basic QC across several metrics. 
#' @details This is Step 1 in the basic Single-Cell RNA-seq workflow.
#'  Returns data as a Seurat Object, the basic data structure for 
#'  Seurat Single Cell analysis. 
#' @param input Vector of scRNA-Seq H5 files, including file path.
#' @param organism Please select species. Choices are Human or Mouse.
#'  (Default is human)
#' @param rename If FALSE, keep original sample names. If TRUE, rename samples
#'  (see below to input your new sample names).
#' @param new.sample.names Vector of sample names to replace originals in the
#'  order of H5 files. eg c("Sample_1", "Sample_2", "Sample_3")
#' @param min.cells Filter out genes found in less than this number of cells.
#'  E.g. Setting to 3 will remove genes found in fewer than 3 cells of a sample.
#'  (Default is 3)
#' @param min.genes Filter out cells with less than this number of genes found
#'  in them. E.g. setting to 200 will remove cells that have fewer than 
#'  200 genes from those analyzed for each sample (Default is 200).
#' @param complexity Number of genes detected per UMI. The more genes detected
#'  per UMI, the more complex the data. Cells that have a high number of UMIs 
#'  but only a low number of genes could be dying cells, but also could 
#'  represent a population of a low complexity cell type (i.e red blood cells). 
#'  We suggest that you increase to 0.8 if samples have suspected 
#'  RBC contamination. (Default is 0.5)
#' @param mad.genes.value How many Median Absolute Deviations do you want to
#'  use to filter out cells with too many genes? For example, entering "3" will 
#'  remove all cells with 3 absolute deviations greater a number of genes than 
#'  the median cell is calculated to have. (Default is 3)
#' @param mad.mitoch.value How many Median Absolute Deviations do you want to
#'  use to filter out cells with too high a percentage of mitochondrial RNA?
#'  For example, entering "3" will remove all cells with 3 absolute deviations
#'  greater a percent mitonchondrial content than the median cell is calculated 
#'  to have. (Default is NA)
#' @param min.umi Filter out cells with low number of reads  (Default is 500)
#' @param mad.gene Filter by number of genes: If TRUE, uses median absolute
#'  deviation to detect outliers. If FALSE, uses an absolute threshold set 
#'  below (see Filter maximum number of genes) (Default TRUE)
#' @param mad.mitoch Filter by mitochondrial percentage: If TRUE, uses median
#'  absolute deviation to remove outliers. If FALSE, uses set value 
#'  (below, Filter maximum percentage mitochondrial content) (Default TRUE)
#' @param max.genes To remove potential doublets, set maximum number of genes
#'  per cell. E.g. Setting to 2,500 will remove cells with 
#'  more than 2,500 genes. (Default is 2500)
#' @param max.mitoch Filter out cells whose proportion of mitochondrial genes
#'  exceed this threshold. E.g. setting to 10 removes cells with more 
#'  than 10 mitochondrial RNA. (Default is 10)
#' @param filter.vdj.genes If FALSE to remove VDJ genes from the scRNA
#'  transcriptome assay. This is to prevent clustering bias in T-cells of the 
#'  same clonotype. 
#'  Only recommended if you are also doing TCR-seq. (Default is FALSE)
#' @param keep If TRUE, keep files when pattern is found in sample name.
#'  If FALSE, remove files when pattern is found in sample name.
#'  The pattern is set in the file.filter.regex parameter (below).
#' @param file.filter.regex Pattern or regular expression in file 
#' name. Use the keep parameter (above) to keep or remove files
#' that contain pattern.
#' @param split.h5 If TRUE, split H5 into individual files. (Default is FALSE)
#' @param protein If TRUE, dataset contains CITE-Seq data (default is FALSE)
#' @param cell.hash If TRUE, dataset contains cell hashtags. (Default is FALSE)
#' @param plot.histogram Set to TRUE to plot QC graphs as histograms.
#'  Keep as FALSE to plot as violin plots. (Default is FALSE)
#' 
#' @import Seurat 
#' @import reshape2
#' @import tidyverse 
#' @import RColorBrewer
#' @import stringr
#' @import svglite 
#' @import ggplot2
#' @import png
#' @import grid
#' @import svglite
#' @importFrom Seurat CreateAssayObject
#' @importFrom gridExtra arrangeGrob
#' @importFrom Seurat Idents
#' @importFrom svglite svglite
#' @importFrom digest digest

#' 
#' @export
#' 
#' @return Seurat Object and QC plots

filterQC <- function(input,
                     organism = "Human",
                     rename = FALSE,
                     new.sample.names = NULL,
                     min.cells = 3,
                     min.genes = 200,
                     complexity = 0.5,
                     mad.genes.value = 3,
                     mad.mitoch.value = NA,
                     min.umi = 500,
                     mad.gene = TRUE,
                     mad.mitoch = TRUE,
                     max.genes = 2500,
                     max.mitoch = 10,
                     filter.vdj.genes = FALSE,
                     keep = TRUE,
                     file.filter.regex = "",
                     split.h5 = FALSE,
                     protein = FALSE,
                     cell.hash = FALSE,
                     plot.histogram = FALSE
){
  
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  ###################################
  ## Process SO object
  
  .seuratObject <- function(i) {
    ## Normalize Data
    so.nf <- so.orig.nf[[i]]
    so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
    so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, 
                                                  pattern = mitoch)
    so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA)/log10(so.nf$nCount_RNA)
    
    ## Detect CITEseq
    if ("protein" %in% names(so.nf)){
      so.nf <- NormalizeData(so.nf, assay = "protein",
                             normalization.method = "CLR")
    }
    
    ## Detect HTO data
    if ("HTO" %in% names(so.nf)){
      so.nf <- NormalizeData(so.nf, assay = "HTO", normalization.method = "CLR")
    }
    
    so <- so.nf
    
    ## Filter out VDJ genes
    if (filter.vdj.genes==TRUE) {
      allGenes = rownames(so)
      VDJgenes = c("TRBV","TRAV","TRBD","TRAJ","TRBJ")
      print("Removing VDJ genes. Genes removed...")
      for (j in VDJgenes) {
        print(allGenes[grepl(j, allGenes)])
        allGenes = allGenes[!grepl(j, allGenes)]  
      }
      so <- subset(so,features = allGenes)
    }
    
    cat("\n\n")
    cat(names(obj.list)[i],":\n")
    so.origcount = dim(so.nf)[2]
    cat(paste0("Original Cell Count=", so.origcount),"\n")
    
    ###################################
    ## Start with filtering here:
    
    ngenestdev <- mad(so@meta.data$nFeature_RNA)
    ngenemed <- median(so@meta.data$nFeature_RNA)
    ngenemaxlim <- ngenemed+(mad.genes.value*ngenestdev)
    gl <- format(round(ngenemaxlim,0),nsmall=0)
    
    
    mitostdev <- mad(so@meta.data$percent.mt)
    mitomed <- median(so@meta.data$percent.mt)
    mitomaxlim <- mitomed+(mad.mitoch.value*mitostdev)
    ml <- format(round(mitomaxlim,2),nsmall=2)
    ## Use Median Absolute Deviations for standard genes and 
    ## mitochondrial Genes to filter data
    if (mad.gene == TRUE & mad.mitoch == TRUE) {
      cat(paste0("Gene Count Filter = low:",min.genes," high:",gl),"\n")
      cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
      cat(paste0("Complexity Filter =",complexity,"\n"))
      so <- subset(so, cells = 
                     rownames(so@meta.data[which(
                       so@meta.data$nFeature_RNA < ngenemaxlim & 
                         so@meta.data$percent.mt <= mitomaxlim & 
                         so@meta.data$log10GenesPerUMI > complexity), ]
                     ))
      perc.remain = (dim(so)[2]/so.origcount)*100
      perc.remain=formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
      
      ##Use only Median Absolute Deviations for mitochondrial Genes to filter data  
    }else if (mad.gene == FALSE & mad.mitoch == TRUE) {
      cat(paste0("Gene Count Filter = low:", min.genes," high:", max.genes),"\n")
      cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
      so <- subset(so, cells = 
                     rownames(so@meta.data[which(
                       so@meta.data$nFeature_RNA < max.genes & 
                         so@meta.data$percent.mt <= mitomaxlim & 
                         so@meta.data$log10GenesPerUMI > complexity & 
                         so@meta.data$nCount_RNA > min.umi), ]
                     ))
      perc.remain <- (dim(so)[2]/so.origcount)*100
      perc.remain=formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
      
      ##Use only Median Absolute Deviations for standard Genes to filter data  
    }else if (mad.gene == TRUE & mad.mitoch == FALSE){
      cat(paste0("Gene Count Filter = low:",min.genes," high:",gl),"\n")
      cat(paste0("Mitochondrial Percentage Filter =", max.mitoch,"\n"))
      so <- subset(so, cells = 
                     rownames(so@meta.data[which(
                       so@meta.data$nFeature_RNA < ngenemaxlim & 
                         so@meta.data$nFeature_RNA > min.genes & 
                         so@meta.data$percent.mt < max.mitoch & 
                         so@meta.data$log10GenesPerUMI > complexity & 
                         so@meta.data$nCount_RNA > min.umi), ]
                     ))
      perc.remain <- (dim(so)[2]/so.origcount)*100
      perc.remain <- formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
      
      ##Dont use Median Absolute Deviations to filter data  
    }else{
      cat(paste0("Gene Count Filter = low:", min.genes," high:", max.genes),"\n")
      cat(paste0("Mitochondrial Percentage Filter =", max.mitoch,"\n"))
      so <- subset(so, cells = 
                     rownames(so@meta.data[which(
                       so@meta.data$nFeature_RNA < max.genes & 
                         so@meta.data$nFeature_RNA > min.genes & 
                         so@meta.data$percent.mt < max.mitoch & 
                         so@meta.data$log10GenesPerUMI > 
                         complexity & so@meta.data$nCount_RNA > min.umi), ]
                     ))
      perc.remain <- (dim(so)[2]/so.origcount)*100
      perc.remain <- formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
    }
    
    
    ## Set up data for Violin or Histogram plot  
    .runPlots <- function(x,name){
      count.df <- df.m %>% dplyr::filter(variable == x) 
      count2.df <- df2.m %>%  plyr::filter(variable == x) 
      
      qc.df <- array(0,dim=c(0,4))
      qc.df <- rbind(qc.df,count2.df,count.df)
      
      if(plot.histogram==TRUE){
        gg <- .plotHist(qc.df,name)
      }else{
        gg <- .plotViolin(qc.df,name)
      }
    }
    
    ## Set up data for Scatter plot  
    .runScatter <- function(x,name){
      x <- as.character(x)
      scplot.m <- so@meta.data %>% 
        dplyr::select("nCount_RNA",x) %>% 
        dplyr::mutate(filt = "filt")
      
      scplot2.m <- so.nf@meta.data %>% 
        dplyr::select("nCount_RNA",x) %>% 
        dplyr::mutate(filt = "raw") 
      
      sc.plot.all <- rbind(scplot2.m,scplot.m)
      
      g <- ggplot(sc.plot.all,aes_string(x="nCount_RNA",y=x,color="filt")) + 
        geom_point(size = 0.5) + 
        theme_classic() +
        ggtitle(paste(name)) 
      
      return(g)
    }
    
    
    ## input data for plots
    df.m <- melt(so@meta.data)
    df.m$filt <- "filt"
    df.m$filt <- as.factor(df.m$filt)
    df2.m <- melt(so.nf@meta.data)
    df2.m$filt <- "raw"
    df2.m$filt <- as.factor(df2.m$filt)
    
    ## create plots
    v <- unique(df.m$variable)
    grob.list <- lapply(v,function(x){.runPlots(x,so@project.name)})
    grob2.list <- lapply(v,function(x){.runScatter(x, so@project.name)})
    grob.all <- arrangeGrob(grobs = grob.list, ncol = length(grob.list))
    grob2.all <- arrangeGrob(grobs = grob2.list, ncol = length(grob2.list))
    so2.list <- list(so,so.nf,grob.all,grob2.all)
    
    return(so2.list)
  }
  
  
  ###################################
  ## Plotting Functions
  
  ## Create Historgram
  .plotHist <- function(count.df,name){ 
    g <- ggplot(count.df,aes(x=value,fill=filt)) + 
      theme_bw() +
      geom_histogram(binwidth=.05, alpha = 0.7, position="identity") +
      scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
      scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
      labs(x = NULL) +
      theme(plot.title = element_text(size=6),
            legend.position='right',
            legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      ggtitle(paste(name,count.df$variable[1])) +
      scale_x_continuous(trans='log10') + 
      scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),6))
    return(g)
  }
  
  ## Create Violin Plot
  .plotViolin <- function(count.df,name){ 
    axislab <- unique(count.df$filt)
    col1 <- brewer.pal(8, "Set3")[-2] 
    col2 <- c(col1,brewer.pal(8,"Set2")[3:6])
    v <- ggplot(count.df, aes(x=filt, y=value)) +
      ggtitle(paste(name,count.df$variable[1])) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text=element_text(size=rel(1.5)),
            legend.title=element_blank(), axis.text=element_text(size=10),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            plot.title = element_text(size = 12, face = "bold")) +
      geom_violin(aes(fill=as.factor(filt))) +  
      scale_fill_manual(values = c("#00AFBB", "#FC4E07")) + 
      geom_boxplot(width=.1) +
      scale_x_discrete(limits = as.vector(axislab)) 
    return(v)
    
  }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  ## Get Cell Cycle information
  if (organism == "Human"){
    mitoch = "^MT-"
  }else{
    mitoch = "^mt-"
    cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
    cc.genes$s.genes = str_to_title(cc.genes$s.genes)
  }
  
  ###################################
  ## Process H5 files
  
  obj.list <- lapply(input, 
                     function(x) { return(Read10X_h5(x, use.names=TRUE)) })
  
  
  ## Rename Samples   
  if (rename == FALSE){
    names(obj.list) <- lapply(input, basename)
    names(obj.list) <- sapply(names(obj.list), 
                              function(x) gsub("_filtered(\\w+)?.h5","", x))
    names(obj.list) <- sapply(names(obj.list), 
                              function(x) gsub("\\.(\\w+)?","", x))
  }else{
    names(obj.list) <- new.sample.names
    obj.list <- obj.list[sort(names(obj.list))]
  }
  
  ## Remove Sample files
  subsetRegex <- file.filter.regex
  if (length(subsetRegex) > 0) {
    if (keep == TRUE){
      for (i in length(subsetRegex)) {
        obj.list <- obj.list[grepl(subsetRegex[[i]],names(obj.list))]
      }
    }else{
      for (i in length(subsetRegex)) {
        obj.list <- obj.list[!grepl(subsetRegex[[i]],names(obj.list))]
      }
    }
  }
  
  
  ## Create Seurat Object from original H5, splitting to multiple ones as necessary.  
  ## For splitting H5's only RNA slot is expected and supported.
  
  if(split.h5 == TRUE){
    i <- 1
    if (class(obj.list[[i]]) == "dgCMatrix"){
      so.orig.nf <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", 
                                       project=names(obj.list)[[i]], 
                                       min.cells = min.cells)
    }else{
      so.orig.nf <- CreateSeuratObject(counts=obj.list[[i]][1]$`Gene Expression`, 
                                       assay="RNA", project=names(obj.list)[[i]], 
                                       min.cells = min.cells)
    }
    so.orig.nf <- SplitObject(so.orig.nf, split.by = "ident")
  }else{
    ## Create Seurat Object
    so.orig.nf <- list()
    for(i in seq_along(names(obj.list))){
      ## From dgCMatrix
      if (class(obj.list[[i]]) == "dgCMatrix"){
        so.orig.nf[[i]] <- CreateSeuratObject(counts = obj.list[[i]], 
                                              assay = "RNA", 
                                              project=names(obj.list)[[i]], 
                                              min.cells = min.cells)
      }else{
        ## From gene Expression Matrix
        k <- names(obj.list[[i]])
        for(j in 1:length(k)){
          if(names(obj.list[[i]][j]) == "Gene Expression"){
            so.orig.nf[[i]] <- CreateSeuratObject(counts = obj.list[[i]][k][[j]], 
                                                  assay = "RNA", 
                                                  project=names(obj.list)[[i]], 
                                                  min.cells = min.cells)
            
          }else if(names(obj.list[[i]][j]) == "Antibody Capture"){
            ## CITEseq data and HTO data
            protein <- rownames(
              obj.list[[i]][k][[j]])[
                !grepl("HTO*",rownames(obj.list[[i]][k][[j]]))]
            HTO <- rownames(
              obj.list[[i]][k][[j]])[
                grepl("HTO*",rownames(obj.list[[i]][k][[j]]))]
            
            so.orig.nf[[i]][["protein"]] <- 
              CreateAssayObject(obj.list[[i]][k][[j]][protein, 
                                                      colnames(so.orig.nf[[i]]
                                                      )])
            
            if(length(HTO)>0){
              so.orig.nf[[i]][['HTO']] <- 
                CreateAssayObject(counts=obj.list[[i]][k][[j]][HTO, 
                                                               colnames(
                                                                 so.orig.nf[[i]]
                                                               )])}
          }else{
            ## Error Report if h5 not correctly formated
            print(paste(names(obj.list[[i]][j]),"found, not stored"))
          }
        }
      }
      names(so.orig.nf)[[i]] <- names(obj.list)[[i]]
    }
  }
  
  #############################
  ## Run Filtering Function
  so.list <- lapply(seq_along(so.orig.nf), .seuratObject)
  
  ### Filter Data
  so.f.list <- lapply(so.list,function(x) x[[1]])
  names(so.f.list) <- unlist(lapply(so.list, 
                                    function(x)as.character(
                                      Seurat::Idents(x[[1]])[1])))
  
  ### No filter Filter Data
  so.nf.list <- lapply(so.list,function(x) x[[2]])
  names(so.nf.list) <- unlist(lapply(so.list, 
                                     function(x)as.character(
                                       Seurat::Idents(x[[1]])[1])))
  
  so.final.list <- so.f.list
  
  so.grobs.list <- lapply(so.list,function(x) x[[3]])
  so.grobs2.list <- lapply(so.list,function(x) x[[4]])
  
  cat("Final filtered samples:\n")
  print(so.f.list)
  cat("Final filtered sample names:\n")
  print(names(so.f.list))   
  
  #############################
  ## Plot Image 
  
  grobdat <- list()
  for(i in 1:length(so.grobs.list)){
    grobdat <- append(grobdat,list(so.grobs.list[[i]]))}
  for(i in 1:length(so.grobs2.list)){
    grobdat <- append(grobdat,list(so.grobs2.list[[i]]))}
  
  grob <- arrangeGrob(grobs=grobdat,nrow=length(grobdat))
  
  
  #############################
  ## Output
  
  cellcount.nf <- lapply(so.nf.list, function(x) dim(x)[2])
  cellcount.f <- lapply(so.f.list, function(x) dim(x)[2])
  sum.before <- sum(unlist(cellcount.nf))
  sum.after <- sum(unlist(cellcount.f))
  cat("\n\nTotal number of cells before filtering:", sum.before, "\n")
  cat("Total number of cells after filtering: ", sum.after,"\n\n")
  cat("Percentage cells remaining after filtering:", 
      (sum.after/sum.before)*100,"\n")
  
  rm(so.list)
  rm(so.nf.list)
  rm(so.f.list)
  
  gc(full = TRUE) 
  
  return(list(object=so.final.list,plot=grob))
  
  
}










