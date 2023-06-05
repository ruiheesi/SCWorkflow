#' @title Process Raw Data
#' @description Creates Seurat Objects from h5 files for individual or
#' merged samples. Will log normalize and produce QC figures for 
#' individual samples
#' @details This is Step 1 in the basic Single-Cell RNA-seq workflow.
#'  Returns data as a Seurat Object, the basic data structure for 
#'  Seurat Single Cell analysis.
#'   
#' @param input Input can be a vector of scRNA-Seq .h5 files, or a list of 
#'  seurat objects for each sample. vector should include file path.
#' @param sample.metadata.table A table of sample metadata that you want to 
#' append to the already-existing metadata within the input Seurat Object(s). 
#' (optional)
#' @param sample.name.column The column of the input metadata.to.append table
#' that contains sample names matching the orig.idents in the input object(s).
#' (optional)
#' @param organism Please select species. Choices are Human or Mouse.
#'  (Default: Human).
#' @param rename.col Select column name from metadata table that contains new 
#' samples name (optional).
#' @param keep If TRUE, keep files when pattern is found in sample name.
#'  If FALSE, remove files when pattern is found in sample name.
#'  The pattern is set in the file.filter.regex parameter (below).
#' @param file.filter.regex Pattern or regular expression in file 
#' name. Use the keep parameter (above) to keep or remove files
#' that contain pattern. If samples have been renamed set
#' regular expression based on new names
#' @param split.h5 If TRUE, split H5 into individual files. (Default: FALSE)
#' @param cell.hash If TRUE, dataset contains cell hashtags. (Default: FALSE)
#' @param do.normalize.data If TRUE counts table will be log2 normalized. If 
#' input contains counts that are already normalzed set to FALSE. 
#' (Default: TRUE)
#' 
#' 
#' @import Seurat 
#' @import reshape2
#' @import scDblFinder
#' @import tidyverse 
#' @import RColorBrewer
#' @import stringr
#' @import svglite 
#' @import ggplot2
#' @import ggpubr
#' @import grid
#' @import png
#' @import svglite
#' @importFrom Seurat CreateAssayObject
#' @importFrom svglite svglite
#' @importFrom digest digest
#' 
#' 

processRawData <- function(input,
                           sample.metadata.table=NULL,
                           sample.name.column=NULL,
                           organism,
                           rename.col=NULL,
                           keep=T,
                           file.filter.regex=c(),
                           split.h5=F,
                           cell.hash=F,
                           do.normalize.data=T                
){          
  ## --------- ##
  ## Functions ####
  ## --------- ##
  
  # Cell Cycle Scoring and Find Variable Features
  CC_FVF_so <- function(so){
    so <- CellCycleScoring(object = so, 
                           g2m.features = cc.genes$g2m.genes,
                           s.features = cc.genes$s.genes)
    so$CC.Difference <- so$S.Score - so$G2M.Score
    return(so)
  }
  
  #### Figures ####
  
  .plotScatterPost2=function(count.df,xaxis,yaxis){	
    ylab = as.character(xaxis)	
    xlab = as.character(yaxis)	
    name = paste(ylab,"vs.",xlab)          
    g = ggplot(count.df, aes(x=.data[[xaxis]], y=.data[[yaxis]],color=Sample)) +
      geom_point(size = 0.5) + 
      theme_classic() +
      theme(legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      guides(colour = guide_legend(override.aes = list(size=2))) +
      scale_color_manual(values = col2) +
      labs( x = xlab, y = ylab)
    
    return(g)
  }
  
  .plotHistPost2 <- function(count.df,xaxis){	
    g=ggplot(count.df) + 
      theme_bw() +
      geom_density(aes(x = .data[[xaxis]], colour = Sample)) +
      # labs(x = NULL) +
      theme(legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      # ggtitle(xaxis) +
      scale_x_continuous(trans='log10') + 
      scale_color_manual(values = col2) %>% 
      suppressMessages()%>%suppressWarnings()
    return(g)
  }
  
  .plotViolinPost2=function(count.df,yaxis){
    axis.lab = unique(count.df$Sample)
    
    g=ggplot(count.df, aes_string(x='Sample', y=(yaxis))) +
      ggtitle(yaxis) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            legend.text=element_text(size=rel(1)),
            legend.title=element_blank(), 
            axis.text=element_text(size=10),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            #axis.text.x=element_text(angle=45,hjust=1),
            plot.title = element_text(size = 20, face = "bold")) +
      geom_violin(aes(fill=as.factor(Sample))) +  
      scale_fill_manual(values = col2) +
      geom_boxplot(width=.1) +
      scale_x_discrete(limits = as.vector(axis.lab)) 
    return(g)
    
  }
  
  
  ### log Normalized data
  .logNormSeuratObject <- function(i) {
    so.nf <- so.orig.nf[[i]]
    
    ## Normalize Data
    so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
    
    
    ## Detect & Normalize CITEseq 
    if ("protein" %in% names(so.nf)){
      so.nf <- NormalizeData(so.nf, assay = "protein",
                             normalization.method = "CLR")
    }
    
    ## Detect & Normalize HTO data
    if ("HTO" %in% names(so.nf)){
      so.nf <- NormalizeData(so.nf, assay = "HTO", normalization.method = "CLR")
    }
    
    return(so.nf)
  }
  
  
  ### Add Metrics to Metadata
  .calcMetrics <- function(i) {
    so.nf <- so.orig.nf[[i]]
    
    ## calculate Percent Mito
    so.nf[["percent.mt"]] = PercentageFeatureSet(object = so.nf, 
                                                 pattern = mitoch)
    ## calculate Genes per UMI
    so.nf$log10GenesPerUMI = log10(so.nf$nFeature_RNA)/log10(so.nf$nCount_RNA)
    
    return(so.nf)
  }
  
  
  
  ## --------------- ##
  ## Main Code Block ####
  ## --------------- ##
  
  
  ## Get Cell Cycle information
  if (organism == "Human"){
    mitoch = "^MT-"
  }else{
    mitoch = "^mt-"
    cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
    cc.genes$s.genes = str_to_title(cc.genes$s.genes)
  }
  
  
  ### Process files h5, rds ####
  
  ### Create SO object depending on class of input SOs. 
  if(class(input)=='RFilePaths'){
    print(paste0('File Type: ',class(input)))
    input <- input$value[grepl(".h5",input$value)]
    
    obj.list <- lapply(input, 
                       function(x) { return(Read10X_h5(x, use.names=TRUE)) })
  }else if(class(input)=='FoundryTransformInput'){
    print(paste0('File Type: ',class(input)))
    
    input=nidapGetFiles(input,'h5')
    
    obj.list <- lapply(input, 
                       function(x) { return(Read10X_h5(x, use.names=TRUE)) })
    
  } else if(class(input)=='character'){
    if (sum(grepl('.rds',input))==1) {
      ## Log output.
      cat("1. Reading Seurat Object from dataset: seurat_object.rds\n\n")
      
      obj.list <- readRDS(input)
      
    } else if (sum(grepl('.h5',input))>0){
      ## Log output.
      cat("1. Processing .h5 files from dataset \n\n")
      
      obj.list <- lapply(input, 
                         function(x) { return(Read10X_h5(x, use.names=TRUE)) })
      
    }else {
      stop("Incorrect input format")
    }
  } else {
    stop("Incorrect input format")
  }
  
  
  ## Clean up sample names
  names(obj.list) <- lapply(input, basename)
  names(obj.list) <- sapply(names(obj.list), 
                            function(x) gsub("_filtered(\\w+)?.h5","", x))
  names(obj.list) <- sapply(names(obj.list), 
                            function(x) gsub("\\.(\\w+)?","", x))
  
  ### Process Metadata table ####
  if(is.null(sample.metadata.table)==F){
    meta_class <- getClass(class(sample.metadata.table)) 
    if (meta_class@className=='character'){
      
      meta.table=read.delim(file = sample.metadata.table,
                            header = T,sep = '\t')%>%
                            suppressWarnings()
      
    } else {
      meta.table=sample.metadata.table
    }
  } else { print("No Metadata provided")}
  
  ### Create Seurat Object ####
  so.orig.nf <- list()
  for(i in seq_along(names(obj.list))){
    ## From dgCMatrix
    if (class(obj.list[[i]]) == "dgCMatrix"){
      so.orig.nf[[i]] <- CreateSeuratObject(counts = obj.list[[i]], 
                                            assay = "RNA", 
                                            project=names(obj.list)[[i]], 
                                            min.cells = 0)
    }else{
      ## From gene Expression Matrix
      k <- names(obj.list[[i]])
      for(j in 1:length(k)){
        if(names(obj.list[[i]][j]) == "Gene Expression"){
          so.orig.nf[[i]] = CreateSeuratObject(counts = obj.list[[i]][k][[j]], 
                                               assay = "RNA", 
                                               project=names(obj.list)[[i]], 
                                               min.cells = 0)
          
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
            so.orig.nf[[i]][['HTO']] = 
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
  
  
  ### Split H5 ####
  
  if(split.h5 == TRUE){
    if (length(so.orig.nf)==1) {
      cat('Splitting seurat object by idents')
      so.orig.nf <- SplitObject(so.orig.nf[[1]], split.by = "ident")
      
      cat('Split Sample Names:\n',paste(names(so.orig.nf),collapse = '\n'))
      
    } else {
      stop('imported data is list of muiltiple objects and cannot be split: 
    set split.h5 to FALSE'
      )
    }
  }
  sample.names=names(so.orig.nf)
  
  
  ### Rename Samples ####
  if(is.null(sample.metadata.table)==F&is.null(rename.col)==F){
    if(sample.name.column!=rename.col){
      if(identical(sort(names(so.orig.nf)), 
                   sort(meta.table[[sample.name.column]]))){
        for (i in names(so.orig.nf)) {
          
        nname=meta.table[meta.table[,sample.name.column]%in%i,]
        ## add original name to metadata table
        so.orig.nf[[i]]@meta.data[,sample.name.column]=
          nname[,sample.name.column]
        ## change orig.ident col to new name
        so.orig.nf[[i]]@meta.data$orig.ident=nname[,rename.col]
        
        names(so.orig.nf)[names(so.orig.nf)%in%i]=nname[,rename.col]
        }
        ## remove original sample names from meta.data table
        meta.table=meta.table[,!colnames(meta.table)%in%
                                sample.name.column,drop=F]
        sample.name.column=rename.col
        sample.names=names(so.orig.nf)
        cat("Renamed Samples:\n",paste(names(so.orig.nf),collapse = '\n'))
        
        
      }else {
        stop("ERROR: Your metadata sample names do not mach the sample 
                      names of your samples")
      }
    }else{
      cat("Sample Names:\n",paste(names(so.orig.nf),collapse = '\n'))
    }
  }else{
    cat("Sample Names:\n",paste(names(so.orig.nf),collapse = '\n'))
  }
  
  
  
  
  ### log Normalize Data ####
  if (do.normalize.data) {
    so.orig.nf <- lapply(seq_along(so.orig.nf), .logNormSeuratObject)
    names(so.orig.nf)=sample.names
  }else{
    print('Did not run Normalization, Input data is already Log Normalized')
  }   
  
  ### Calculate metrics ####
  so.orig.nf <- lapply(seq_along(so.orig.nf), .calcMetrics)
  names(so.orig.nf)=sample.names
  
  ### add cell cycle information
  so.orig.nf <- lapply(so.orig.nf, CC_FVF_so) %>%suppressWarnings()
  
  ### Add Metadata ####
  if(is.null(sample.metadata.table)==F ){
    metacols <- colnames(meta.table)
    metacols <- metacols[!metacols %in% 
                           unique(c(rename.col,sample.name.column))]
    if (length(metacols)>0) {
      
      so.orig.nf=appendMetadataToSeuratObject(
        so.orig.nf,
        meta.table[,c(sample.name.column,metacols)],
        sample.name.column)[['object']]
      
    } else { print("No Metadata Columns were included in Metadata table")}
  } else { print("No Metadata table was supplied")}
  
  
  
  
  
  ### Remove Sample files ####
  subsetRegex <- file.filter.regex
  if (length(subsetRegex) > 0) {
    if (keep == TRUE){
      for (i in length(subsetRegex)) {
        so.orig.nf <- so.orig.nf[grepl(subsetRegex[[i]],names(so.orig.nf))]
      }
    }else{
      for (i in length(subsetRegex)) {
        so.orig.nf <- so.orig.nf[!grepl(subsetRegex[[i]],names(so.orig.nf))]
      }
    }
  }
  
  
  ### Create QC figure ####
  col1 <- brewer.pal(8, "Set3")[-2] 
  col2 <- c(col1,brewer.pal(8,"Set2")[3:6])
  
  features=c("orig.ident",
             "nCount_RNA",
             "nFeature_RNA",
             "percent.mt",
             "log10GenesPerUMI")
  v=features[features%in%c('nCount_RNA',
                           'nFeature_RNA',
                           'percent.mt',
                           'log10GenesPerUMI')]
  
  
  #### Combine SO meta.data tables ####
  table.meta=data.frame()
  for(s in names(so.orig.nf)) {
    ftable=so.orig.nf[[s]]@meta.data
    ftable=ftable[,colnames(ftable)%in%features,drop=F]
    ftable$Sample=s
    
    table.meta=rbind(table.meta,ftable)
  }
  
  table.meta$nFeature_RNA=as.numeric(table.meta$nFeature_RNA)
  
  
  ### Post Filter Summary - Scatter
  
  scatter.allsamples=lapply(v,
                      function(y){.plotScatterPost2(table.meta,'nCount_RNA',y)})
  names(scatter.allsamples)=v
  
  scatter.allsamples.grob=ggarrange(plotlist=scatter.allsamples,
                                    ncol=1,
                                    common.legend = T,
                                    legend = 'right')
  
  
  ### Post Filter Summary - histogram
  
  hist.allsamples=lapply(v,function(x){.plotHistPost2(table.meta,x)})
  names(hist.allsamples)=v
  
  hist.allsamples.grob=ggarrange(plotlist=hist.allsamples,
                                 ncol=1,
                                 common.legend = T,
                                 legend = 'right') %>% 
    suppressMessages()%>%suppressWarnings()
  
  
  ### Post Filter Summary - Violin
  
  violin.allsamples=lapply(v,function(x){.plotViolinPost2(table.meta,x)})
  names(violin.allsamples)=v
  violin.allsamples.grob=ggarrange(plotlist=violin.allsamples,
                                   ncol=1,
                                   common.legend = T,
                                   legend = 'right')
  
  violin.allsamples.grob=annotate_figure(violin.allsamples.grob, 
                                     top = text_grob("", 
                                                     face = "bold", size = 14))
  
  ### Post Filter Summary - combined Scatter + Histogram
  
  raw.grobs=ggarrange(
    ggarrange(plotlist=scatter.allsamples,
              ncol=1,legend = 'none'),
    ggarrange(plotlist=hist.allsamples,
              ncol=1,legend = 'none'),
    # ggarrange(plotlist=violin.allsamples,ncol=1,legend = 'none'),
    legend.grob=get_legend(scatter.allsamples[[1]]),
    ncol=2,
    legend='right') %>% 
    suppressMessages()%>%suppressWarnings()
  
  raw.grobs=annotate_figure(raw.grobs, 
                            top = text_grob(" Summary ", 
                                            face = "bold", size = 14))
  
  
  
  ### Output
  return(
    list(
      object=so.orig.nf,
      plots=list(
        CombinedQC=raw.grobs,
        ViolinQC=violin.allsamples.grob,
        ScatterQC=scatter.allsamples.grob,
        HistogramQC=hist.allsamples.grob
      )
    )
  )
  
}




