#' Filter and QC template in single-cell-rna-seq-r4 NIDAP environment
#' function from v 167
#' 
#' @title Filter & QC Samples 
#' @description Creates Seurat Objects, Filters and Plots QC (before and after filtering) on samples. Takes a dataset of filtered h5 files (each representing one sample), and performs basic QC across several metrics. Also puts the data into a Seurat Object, the basic data structure for Seurat Single Cell analysis. This template is Step 1 of the basic Single-Cell RNA-seq workflow.
#' @details 
#' 
#' 
#' @param localFilePaths Please input your h5 files that have been uploaded to a single dataset. If you uploaded your h5 files through the Uploader page, this dataset should be the only node in the graph of your code workbook when you first launch a new analysis.
#' @param organism Please select what species your data comes from.
#' @param rename If FALSE, keep original sample names.  If TRUE, rename samples (see below to input your new sample names).
#' @param New_Sample_Names Input a list of new sample names to name your samples. Note: number of new names MUST match the number of samples. For reference, consult log output for Filter QC template for a list of old sample names in order.
#' @param mincells Filter out genes found in less than this number of cells. E.g. Setting to 3 will remove genes found in fewer than three cells of a sample.
#' @param mingenes Filter out cells with less than this number of genes found in them. E.g. setting to 200 will remove cells that have fewer than 200 genes from those analyzed for each sample.
#' @param complexity Number of genes detected per UMI. The more genes detected per UMI, the more complex the data. Cells that have a high number of UMIs but only a low number of genes could be dying cells, but also could represent a population of a low complexity cell type (i.e red blood cells). We suggest that you increase to 0.8 if samples have suspected RBC contamination.
#' @param MAD_genes_value How many Median Absolute Deviations do you want to use to filter out cells with too many genes? For example, entering "3" will remove all cells with 3 absolute deviations greater a number of genes than the median cell is calculated to have.
#' @param MAD_mitoch_value How many Median Absolute Deviations do you want to use to filter out cells with too high a percentage of mitochondrial RNA? For example, entering "3" will remove all cells with 3 absolute deviations greater a percent mitonchondrial content than the median cell is calculated to have.
#' @param minUMI Filter out cells with low number of reads
#' @param MAD_gene Filter by number of genes: If TRUE, uses median absolute deviation to detect outliers. If FALSE, uses an absolute threshold set below (see Filter maximum number of genes)
#' @param MAD_mitoch Filter by mitochondrial percentage: If TRUE, uses median absolute deviation to remove outliers. If FALSE, uses set value (below, Filter maximum percentage mitochondrial content)
#' @param maxgenes To remove potential doublets, set maximum number of genes per cell. E.g. Setting to 2,500 will remove cells with more than 2,500 genes.
#' @param maxmitoch Filter out cells whose proportion of mitochondrial genes exceed this threshold. E.g. setting to 10 removes cells with more than 10% mitochondrial RNA.
#' @param Filter_VDJ_Genes Toggle off to remove VDJ genes from the scRNA transcriptome assay. This is to prevent clustering bias in T-cells of the same clonotype. Only recommended if you are also doing TCR-seq.
#' @param Keep TRUE is Keep, FALSE is Remove when pattern is found in the sample name. The pattern is set in the File Filter Regex parameter (above).
#' @param File_Filter_Regex Filter for pattern or regular expression in file name. Using this, you can filter on the name of the input files. Use the Keep or Remove parameter (below) to choose to keep or remove sample files whose names match this pattern.
#' @param Split_H5 If TRUE, split H5 into individual files.  By default FALSE
#' @param Protein Does this dataset contain CITE-Seq antibodies? Set to TRUE if your dataset does contain CITE-Seq antibodies. Set to FALSE if your dataset does not. The default is FALSE.
#' @param Cell_hash Does this dataset contain Cell Hashing Tags? set to TRUE if your dataset does contain Cell Hashing Tags. Set to FALSE if your dataset does not. The default is FALSE.
#' @param imageType Remember that svgs are much larger than pngs, so we recommend doing everything first in png, then rerunning to output specific svgs as needed.
#' @param plot_histogram Set to TRUE to plot QC graphs as histograms. Keep as FALSE to plot as violin plots.
#' 
#' 
#' @import Seurat 
#' @import reshape2
#' @import tidyverse 
#' @import gridExtra 
#' @import RColorBrewer
#' @import stringr
#' @import svglite 
#' @import ggplot2
#' @import png
#' @import grid
#' @import svglite
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom Seurat Idents
#' @importFrom svglite svglite
#' @importFrom digest digest

#' 
#' @export
#' 
#' @return Seurat Objects, Filters and Plots QC (before and after filtering) on samples. Takes a dataset of filtered h5 files (each representing one sample), and performs basic QC across several metrics. Also puts the data into a Seurat Object, the basic data structure for Seurat Single Cell analysis. This template is Step 1 of the basic Single-Cell RNA-seq workflow.


Filter_and_QC <- function(localFilePaths,
                               organism = "Human",
                               rename = FALSE,
                               New_Sample_Names = c("Sample_1", "Sample_2", "Sample_3"),
                               mincells = 3,
                               mingenes = 200,
                               complexity = 0.5,
                               MAD_genes_value = 3,
                               MAD_mitoch_value = NA,
                               minUMI = 500,
                               MAD_gene = TRUE,
                               MAD_mitoch = TRUE,
                               maxgenes = 2500,
                               maxmitoch = 10,
                               Filter_VDJ_Genes = FALSE,
                               Keep = TRUE,
                               File_Filter_Regex = "",
                               Split_H5 = FALSE,
                               Protein = FALSE,
                               Cell_hash = FALSE,
                               imageType = "png",
                               plot_histogram = FALSE
){
  
  
  
  ## --------- ##
  ## Functions ##
  ## --------- ##
  
  ###################################
  ## Process SO object
  
  seurat_object <- function(i) {
    
    so.nf <- so.orig.nf[[i]]
    so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
    so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
    so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
    
    if ("Protein" %in% names(so.nf)){
      so.nf <- NormalizeData(so.nf, assay = "Protein", normalization.method = "CLR")
    }
    
    if ("HTO" %in% names(so.nf)){
      so.nf <- NormalizeData(so.nf, assay = "HTO", normalization.method = "CLR")
    }
    
    so <- so.nf
    
    if (Filter_VDJ_Genes==TRUE) {
      allGenes = rownames(so)
      VDJgenes = c("TRBV","TRAV","TRBD","TRAJ","TRBJ")
      print("Removing VDJ genes. Genes removed...")
      for (j in VDJgenes) {
        print(allGenes[grepl(j, allGenes)])
        allGenes = allGenes[!grepl(j, allGenes)]  
      }
<<<<<<< HEAD
<<<<<<< HEAD
      so <- subset(so,features = allGenes)
=======
      so <- SubsetData(so,features = allGenes,assay="RNA")
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
=======
      so <- SubsetData(so,features = allGenes,assay="RNA")
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
    }
    
    cat("\n\n")
    cat(names(obj.list)[i],":\n")
    so.origcount = dim(so.nf)[2]
    cat(paste0("Original Cell Count=", so.origcount),"\n")
    
    ###################################
    ## Start with filtering here:
    
    ngenestdev <- mad(so@meta.data$nFeature_RNA)
    ngenemed <- median(so@meta.data$nFeature_RNA)
    ngenemaxlim <- ngenemed+(MAD_genes_value*ngenestdev)
    gl = format(round(ngenemaxlim,0),nsmall=0)
    
    
    mitostdev <- mad(so@meta.data$percent.mt)
    mitomed <- median(so@meta.data$percent.mt)
    mitomaxlim <- mitomed+(MAD_mitoch_value*mitostdev)
    ml = format(round(mitomaxlim,2),nsmall=2)
    
    if (MAD_gene == TRUE & MAD_mitoch == TRUE) {
      cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
      cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
      cat(paste0("Complexity Filter =",complexity,"\n"))
      so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
      perc.remain = (dim(so)[2]/so.origcount)*100
      perc.remain=formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
    } else if (MAD_gene == FALSE & MAD_mitoch == TRUE) {
      cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
      cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
      so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
      perc.remain = (dim(so)[2]/so.origcount)*100
      perc.remain=formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
    } else if (MAD_gene == TRUE & MAD_mitoch == FALSE){
      cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
      cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
      so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
      perc.remain = (dim(so)[2]/so.origcount)*100
      perc.remain=formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
    } else {
      cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
      cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
      so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
      perc.remain = (dim(so)[2]/so.origcount)*100
      perc.remain=formatC(perc.remain,format = "g",digits=3)
      cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
      cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
    }
    
    
    
    Runplots <- function(x,name){
      df.m %>% dplyr::filter(variable == x) -> count.df
      df2.m %>% dplyr::filter(variable == x) -> count2.df
      qc.df <- array(0,dim=c(0,4))
      qc.df <- rbind(qc.df,count2.df,count.df)
      if(plot_histogram==TRUE){
        gg <- plothist(qc.df,name)
        } else{
        gg <- plotviolin(qc.df,name)
      }
    }
    
    RunScatter <- function(x,name){
      x <- as.character(x)
      scplot.m = so@meta.data %>% dplyr::select("nCount_RNA",x) %>% dplyr::mutate(filt = "filt")
      scplot2.m = so.nf@meta.data %>% dplyr::select("nCount_RNA",x) %>% dplyr::mutate(filt = "raw") 
      sc.plot.all = rbind(scplot2.m,scplot.m)
      g=ggplot(sc.plot.all,aes_string(x="nCount_RNA",y=x,color="filt")) + 
        geom_point(size = 0.5) + 
        theme_classic() +
        ggtitle(paste(name)) 
      #print(g)
      #return(NULL)
      return(g)
    }
    
    df.m <- melt(so@meta.data)
    df.m$filt <- "filt"
    df.m$filt <- as.factor(df.m$filt)
    df2.m <- melt(so.nf@meta.data)
    df2.m$filt <- "raw"
    df2.m$filt <- as.factor(df2.m$filt)
    
    v <- unique(df.m$variable)
    grob.list <- lapply(v,function(x){Runplots(x,so@project.name)})
    grob2.list <- lapply(v,function(x){RunScatter(x, so@project.name)})
    grob.all <- arrangeGrob(grobs = grob.list, ncol = length(grob.list))
    grob2.all <- arrangeGrob(grobs = grob2.list, ncol = length(grob2.list))
    so2.list <- list(so,so.nf,grob.all,grob2.all)
    
    return(so2.list)
  }
  
  
  ###################################
  ## Plotting Functions
  
  plothist <- function(count.df,name){
    g=ggplot(count.df,aes(x=value,fill=filt)) + 
      theme_bw() +
      geom_histogram(binwidth=.05, alpha = 0.7, position="identity") +
      #geom_density(aes(x = value, colour = filt, linetype = filt)) +
      scale_fill_manual(values = c("#00AFBB", "#FC4E07")) +
      scale_color_manual(values = c("#00AFBB", "#FC4E07")) +
      labs(x = NULL) +
      theme(plot.title = element_text(size=6),legend.position='right',legend.text=element_text(size=10),
            legend.title=element_blank()) + 
      ggtitle(paste(name,count.df$variable[1])) +
      scale_x_continuous(trans='log10') + 
      scale_linetype_manual(values=rep(c('solid', 'dashed','dotted'),6))
    return(g)
    #print(g)
    #return(NULL)
  }
  
  plotviolin <- function(count.df,name){
    axislab = unique(count.df$filt)
    col1=brewer.pal(8, "Set3")[-2] 
    col2=c(col1,brewer.pal(8,"Set2")[3:6])
    v = ggplot(count.df, aes(x=filt, y=value)) +
      ggtitle(paste(name,count.df$variable[1])) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),legend.text=element_text(size=rel(1.5)),
            legend.title=element_blank(), axis.text=element_text(size=10),
            axis.title.x=element_blank(),axis.title.y=element_blank(),
            axis.text.x=element_text(angle=45,hjust=1),
            plot.title = element_text(size = 12, face = "bold")) +
      geom_violin(aes(fill=as.factor(filt))) +  
      scale_fill_manual(values = c("#00AFBB", "#FC4E07")) + 
      geom_boxplot(width=.1) +
      #labs(colour = n, y=m) +
      #geom_jitter(height = 0, width = 0.1, size = 0.1) +
      #scale_y_continuous(trans='log2') + 
      scale_x_discrete(limits = as.vector(axislab)) 
    return(v)
    #print(v)
    #return(NULL)
  }
  
  
  ## --------------- ##
  ## Main Code Block ##
  ## --------------- ##
  
  
  if (organism == "Human"){
    mitoch = "^MT-"
  } else{
    mitoch = "^mt-"
    # cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
    # cc.genes$s.genes = str_to_title(cc.genes$s.genes)
  }
  
  ###################################
  ## Process H5 files
  
  # if(class(localFilePaths)=='RFilePaths'){
  #   print(paste0('File Type: ',class(localFilePaths)))
  #   localFilePaths <- h5_files$value[grepl(".h5",h5_files$value)]
  #   
  # }else if(class(localFilePaths)=='FoundryTransformInput'){
  #   print(paste0('File Type: ',class(localFilePaths)))
  #   fs <- localFilePaths$fileSystem()
  #   localFilePaths <- fs$ls(glob="*.h5")
  #   localFilePaths <- lapply(localFilePaths, function(file) { return(fs$get_path(file$path, 'r')) })
  #   
  # }else{
  #   stop("\n \n h5 files not correctly imported into Workbook")
  # }
  
  obj.list <- lapply(localFilePaths, function(x) { return(Read10X_h5(x, use.names=TRUE)) })
  
  
  ## Rename    
  if (rename == FALSE){
    names(obj.list) <- lapply(localFilePaths, basename)
    names(obj.list) <- sapply(names(obj.list), function(x) gsub("_filtered(\\w+)?.h5","", x))
    names(obj.list) <- sapply(names(obj.list), function(x) gsub("\\.(\\w+)?","", x))
  }else{
    names(obj.list) <- New_Sample_Names
    obj.list <- obj.list[sort(names(obj.list))]
  }
  
  ## Remove Sample files
<<<<<<< HEAD
<<<<<<< HEAD
  # subsetRegex = eval(parse(text=gsub('\\[\\]','c()',File_Filter_Regex)))
  subsetRegex=File_Filter_Regex
=======
  subsetRegex = eval(parse(text=gsub('\\[\\]','c()',File_Filter_Regex)))
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
=======
  subsetRegex = eval(parse(text=gsub('\\[\\]','c()',File_Filter_Regex)))
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  if (length(subsetRegex) > 0) {
    if (Keep == TRUE){
      for (i in length(subsetRegex)) {
        obj.list <- obj.list[grepl(subsetRegex[[i]],names(obj.list))]
      }
    }else{
      for (i in length(subsetRegex)) {
        obj.list <- obj.list[!grepl(subsetRegex[[i]],names(obj.list))]
      }
    }
  }
  
  
  ## Create Seurat Object from original H5, splitting to multiple ones as necessary.  For splitting H5's only RNA slot is expected and supported.
  
  if(Split_H5 == TRUE){
    i = 1
    if (class(obj.list[[i]]) == "dgCMatrix"){
      so.orig.nf <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells)
    } else {
      so.orig.nf <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells)
    }
    so.orig.nf <- SplitObject(so.orig.nf, split.by = "ident")
  } else {
    so.orig.nf <- list()
    for(i in seq_along(names(obj.list))){
      if (class(obj.list[[i]]) == "dgCMatrix"){
        so.orig.nf[[i]] <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells)
      } else {
        k = names(obj.list[[i]])
        for(j in 1:length(k)){
          if(names(obj.list[[i]][j]) == "Gene Expression"){
            so.orig.nf[[i]] <- CreateSeuratObject(counts = obj.list[[i]][k][[j]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells)
          } else if(names(obj.list[[i]][j]) == "Antibody Capture"){
            Protein = rownames(obj.list[[i]][k][[j]])[!grepl("HTO*",rownames(obj.list[[i]][k][[j]]))]
            HTO = rownames(obj.list[[i]][k][[j]])[grepl("HTO*",rownames(obj.list[[i]][k][[j]]))]
            so.orig.nf[[i]][["Protein"]] <- CreateAssayObject(obj.list[[i]][k][[j]][Protein, colnames(so.orig.nf[[i]])])
            if(length(HTO)>0){
              so.orig.nf[[i]][['HTO']] <- CreateAssayObject(counts = obj.list[[i]][k][[j]][HTO, colnames(so.orig.nf[[i]])])}
          } else{
            print(paste(names(obj.list[[i]][j]),"found, not stored"))
          }
        }
      }
      names(so.orig.nf)[[i]] <- names(obj.list)[[i]]
    }
  }
  
  #############################
  ## Run Filtering Function
  so.list <- lapply(seq_along(so.orig.nf), seurat_object)
  
  so.f.list <- lapply(so.list,function(x) x[[1]])
  names(so.f.list) <- unlist(lapply(so.list, function(x) as.character(Seurat::Idents(x[[1]])[1])))
  
  so.nf.list <- lapply(so.list,function(x) x[[2]])
  names(so.nf.list) <- unlist(lapply(so.list, function(x) as.character(Seurat::Idents(x[[1]])[1])))
  
  so.final.list <- so.f.list
  
  so.grobs.list <- lapply(so.list,function(x) x[[3]])
  so.grobs2.list <- lapply(so.list,function(x) x[[4]])
  
  cat("Final filtered samples:\n")
  print(so.f.list)
  cat("Final filtered sample names:\n")
  print(names(so.f.list))   
  
  #############################
  ## Plot Image 
  
<<<<<<< HEAD
<<<<<<< HEAD
  # ## Set Image Size
  # imageWidth = min(1000*length(so.list[[1]][[3]]),15000)
  # imageHeight = min(1000*length(so.grobs.list)*2,24000)
  # dpi = 300
  # if (imageType == 'png') {
  #   png(
  #     # filename=graphicsFile,
  #     width=imageWidth,
  #     height=imageHeight,
  #     units="px",
  #     pointsize=4,
  #     bg="white",
  #     res=dpi,
  #     type="cairo")
  # } else {
  #   svglite::svglite(
  #     # file=graphicsFile,
  #     width=round(imageWidth/dpi,digits=2),
  #     height=round(imageHeight/dpi,digits=2),
  #     pointsize=1,
  #     bg="white")
  # }
=======
=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  ## Set Image Size
  imageWidth = min(1000*length(so.list[[1]][[3]]),15000)
  imageHeight = min(1000*length(so.grobs.list)*2,24000)
  dpi = 300
  if (imageType == 'png') {
    png(
      # filename=graphicsFile,
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
  } else {
    svglite::svglite(
      # file=graphicsFile,
      width=round(imageWidth/dpi,digits=2),
      height=round(imageHeight/dpi,digits=2),
      pointsize=1,
      bg="white")
  }
<<<<<<< HEAD
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
=======
>>>>>>> 30b7146f0a89f5ab7f8ae790ea33038fe5ca58de
  
  #grid.arrange(grobs = so.grobs.list, nrow = length(so.grobs.list))
  grobdat = list()
  for(i in 1:length(so.grobs.list)){grobdat=append(grobdat,list(so.grobs.list[[i]]))}
  for(i in 1:length(so.grobs2.list)){grobdat=append(grobdat,list(so.grobs2.list[[i]]))}
  
  # grid.arrange(
    grob = arrangeGrob(grobs=grobdat,nrow=length(grobdat))
    # ,nrow=1)
  
  #############################
  ## Output
  
  cellcount.nf <- lapply(so.nf.list, function(x) dim(x)[2])
  cellcount.f <- lapply(so.f.list, function(x) dim(x)[2])
  sum.before <- sum(unlist(cellcount.nf))
  sum.after <- sum(unlist(cellcount.f))
  cat("\n\nTotal number of cells before filtering:", sum.before, "\n")
  cat("Total number of cells after filtering: ", sum.after,"\n\n")
  cat("Percentage cells remaining after filtering:", (sum.after/sum.before)*100,"\n")
  
  rm(so.list)
  rm(so.nf.list)
  rm(so.f.list)
  
  gc(full = TRUE)
  
  return(list(so=so.final.list,plot=grob))
  
  #output <- new.output()
  #output_fs <- output$fileSystem()
  #saveRDS(so.final.list, output_fs$get_path("seurat_object.rds", 'w')) 
  #return(output_fs) 
}










