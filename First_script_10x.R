#### Clear Environnement !!!
############################
rm(list=ls(all=TRUE))  # broom variables
gc()        # garbage collector
cat("\014")   # clear consol

C3M <- T  # whether I'm using the Lab's OS or not

################################################
########                                ########
####       10X Genomics Pre-Analysis        ####
########                                ########
################################################

##  Directories :     ----------------------------------------------------------
if(C3M == T){
  setwd("~/Bureau/Juba/Projects/SingleCell/")
} else {
  setwd("/Projects/T-ALL scRNAseq/")
}

##  Packages :        ----------------------------------------------------------
library(Seurat)
library(dplyr)
library(tibble)
library(ggplot2)

##  Functions :       ----------------------------------------------------------
Implement_VDJ <- function(seurat,
                          VDJ_data){
  
  ## initiate the data
  meta <- seurat@meta.data
  meta$clone_id <- NA
  meta$v_call <- NA
  meta$j_call <- NA
  meta$c_call <- NA
  
  ## implement the values
  for (i in 1:nrow(VDJ_data)) {
    cell <- VDJ_data$cell_id[i]
    if (cell %in% rownames(meta)) {
      meta[cell, "clone_id"] <- VDJ_data$clone_id[i]
      meta[cell, "v_call"] <- VDJ_data$v_call[i]
      meta[cell, "j_call"] <- VDJ_data$j_call[i]
      meta[cell, "c_call"] <- VDJ_data$c_call[i]
    }
  }
  
  ## update the object
  seurat@meta.data <- meta
  
  return(seurat)
}

## Palettes :     --------------------------------------------------------------
MyPalette <- c("#9933aa","#ffdd22","#aa4400","#ff0000","#337722","#00ff66","#005566","#002277",
               "#441144","#aa0077","#00bbff","#003333","#4422cc","#116611","#330077","#111111",
               "#667700","#ddaa00","#33ffff","#ff22ff","#ffff33","#00ff00","#0000ff","#444444")


################################################
########                                ########
####              Get started               ####
########                                ########
################################################

Libraries <- grep("Library", list.files(), value = T)

for(Lib in Libraries){
  
  ##  Setting directories
  Raw_dir <- paste0("data/",Lib,"/")
  Seurat_outdir <- paste0(Lib,"/SeuratObjects/")
  Out_dir_QC <- paste0(Lib,"/Out_figures/QC/")
  Out_dir_fig <- paste0(Lib,"/Out_figures/")
  
  
  ##  Importing 10x data      --------------------------------------
  Rawdata <- Read10X(paste0(Raw_dir,"filtered_feature_bc_matrix/"))
  VDJdata <- read.delim(paste0(Raw_dir,"airr_rearrangement.tsv"), header = T, sep = "\t") %>% 
    dplyr::select(-c(sequence, sequence_aa, sequence_alignment, germline_alignment))
  
  barcodes <- intersect(colnames(Rawdata$`Gene Expression`),
                        colnames(Rawdata$`Antibody Capture`))
  
  umi <- Rawdata$`Gene Expression`[,barcodes]
  hto <- Rawdata$`Antibody Capture`[,barcodes]
  
  ##  Creating seurat object & normalizing RNA data (log normalizing)   -------------
  SeuratObj <- CreateSeuratObject(counts = umi)
  SeuratObj <- NormalizeData(SeuratObj)
  
  ##  Find and scale variable features        ---------------------------------
  SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "mean.var.plot", nfeatures = 2000)
  SeuratObj <- ScaleData(SeuratObj, features = VariableFeatures(SeuratObj))
  
  ##  Add HTO data as an independent assay & normalize      --------------------------
  SeuratObj[["HTO"]] <- CreateAssayObject(counts = hto)
  SeuratObj <- NormalizeData(SeuratObj, assay = "HTO", normalization.method = "CLR")
  
  ##  Demultiplexing the HTO data     --------------------------------------------
  SeuratObj <- HTODemux(SeuratObj, 
                        assay = "HTO", positive.quantile = 0.99, nstarts = 100, nsamples = 100) # default params
  table(SeuratObj$hash.ID)
  table(SeuratObj$HTO_classification.global)
  
  ##  Group cell by their max HTO signal --> identification   -----------------------
  Idents(SeuratObj) <- "HTO_maxID"
  RidgePlot(SeuratObj, assay = "HTO", features = rownames(SeuratObj[["HTO"]]), ncol = 2)
  
  ##  Compare number of UMIs for singlets, doublets and negative cells    ---------------
  Idents(SeuratObj) <- "HTO_classification.global"
  VlnPlot(SeuratObj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
  
  ##  Subset only singlets      -------------------------------------------------------
  SeuratObj <- subset(SeuratObj, idents = "Singlet", invert = FALSE)
  HTOHeatmap(SeuratObj, assay = "HTO", ncells = 3000)
  
  ##  Mitochondiral gene expression : QC      ---------------------------------
  SeuratObj <- PercentageFeatureSet(SeuratObj,
                                    pattern = "^MT-",
                                    col.name = "percent.mt",
                                    assay = "RNA")
  
  ##  Visualization --> Importing as PDF    ---------------------------
  pdf(file = paste0(Out_dir_QC,Lib,"_preprocessing.pdf"), 
      width = 6, 
      height = 6,
      title = "Quality Control")
  
  print(VlnPlot(object = SeuratObj,
                features=c("nCount_RNA","nFeature_RNA","percent.mt"), pt.size = 0.01))
  
  print(ggplot(SeuratObj@meta.data, 
               aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
          geom_point(size = 0.8) +
          scale_y_log10(breaks = c(90,100,300,1000,3000)) +
          scale_color_gradientn(colours = c("green","gold","orange","darkred"),
                                values = c(0,0.2,0.5,1),
                                limits = c(0,50)) +
          ggtitle("QC plot", "Number of detected genes in function of number of UMI") +
          labs(x = "Number of UMI by cell", y = "Number of detected genes by cell")+
          theme_bw())
  dev.off()
  
  ##  Filtering cells and visualization as PDF    -----------------------------
  SeuratObj <- subset(SeuratObj, subset = percent.mt < 5 & nFeature_RNA < 3500)
  
  pdf(file = paste0(Out_dir_QC,Lib,"_preprocessing_filtered.pdf"), 
      width = 6, 
      height = 6,
      title = "Quality Control - filtered")
  
  print(VlnPlot(object = SeuratObj,
                features=c("nCount_RNA","nFeature_RNA","percent.mt"), pt.size = 0.01))
  
  print(ggplot(SeuratObj@meta.data, 
               aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
          geom_point(size = 0.8) +
          scale_y_log10(breaks = c(90,100,300,1000,3000)) +
          scale_color_gradientn(colours = c("green","gold","orange","darkred"),
                                values = c(0,0.2,0.5,1),
                                limits = c(0,50)) +
          ggtitle("QC plot", "Number of detected genes in function of number of UMI") +
          labs(x = "Number of UMI by cell", y = "Number of detected genes by cell")+
          theme_bw())
  dev.off()
  
  ##  SCtransforme the data     -------------------------------------------------
  SeuratObj <- SCTransform(SeuratObj, 
                           vars.to.regress = "percent.mt", 
                           verbose = F,
                           ncells = 5000,
                           variable.features.n = 3000)
  
  ##  Normalization         -----------------------------------------------
  SeuratObj <- NormalizeData(SeuratObj, 
                             normalization.method = "LogNormalize", 
                             scale.factor = median(SeuratObj$nCount_SCT))
  
  ##  Find variable features        -------------------------------------------
  SeuratObj <- FindVariableFeatures(SeuratObj, selection.method = "vst", nfeatures = 2000)
  SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj))
  
  ##  Add vdj data          -----------------------------------------------------
  SeuratObj <- Implement_VDJ(SeuratObj, VDJdata)
  
  ##  Save seuratobject       ------------------------------------------------
  saveRDS(SeuratObj, file = paste0(Seurat_outdir,"Seurat1.rds"))
  
  #   Clustering        ----------------------------------------------------------
  ##  Linear dim reduction      --------------------------------------------------
  SeuratObj <- RunPCA(SeuratObj, 
                      features = VariableFeatures(object = SeuratObj),
                      npcs = 50,
                      ndims.print = 1:2,
                      nfeatures.print = 10)
  
  ##  non-linear dim reduction      -----------------------------------------------
  SeuratObj <- RunUMAP(SeuratObj, dims = 1:30)
  
  ##  Clustering      -------------------------------------------------------------
  SeuratObj <- FindNeighbors(SeuratObj, dims = 1:30, verbose = FALSE)
  SeuratObj <- FindClusters(SeuratObj, verbose = FALSE, resolution = seq(0.1, 1, by= 0.1))
  
  ##  Add umap coordinates to metadata      ---------------------------------------
  SeuratObj@meta.data <- merge(SeuratObj@meta.data, 
                               as.data.frame(SeuratObj@reductions$umap@cell.embeddings),
                               by = 0)
  SeuratObj@meta.data <- SeuratObj@meta.data %>% 
    column_to_rownames(var = "Row.names")
  
  ##  Cell Cycle      --------------------------------------------------------
  SeuratObj <- CellCycleScoring(SeuratObj,
                                s.features = cc.genes.updated.2019$s.genes,
                                g2m.features = cc.genes.updated.2019$g2m.genes)
  
  
  ##  Plotting the clustering as pdf    ------------------------------------------
  pdf(file = paste0(Out_dir_fig,Lib,"_Clusters_no_filt.pdf"), 
      width = 6, 
      height = 5,
      title = "Clusters - no.filt")
  
  print(DimPlot(SeuratObj, 
                group.by = "hash.ID",
                cols = MyPalette))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.1",
                cols = MyPalette))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.3",
                cols = MyPalette))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.5",
                cols = MyPalette))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.7",
                cols = MyPalette))
  print(DimPlot(SeuratObj, 
                group.by = "SCT_snn_res.0.9",
                cols = MyPalette))
  dev.off()
  
  ##  Plotting some metadata (clonotypes) as pdf    ------------------------------------------
  pdf(file = paste0(Out_dir_fig,Lib,"_Clonotypes_no_filt.pdf"), 
      width = 6, 
      height = 5,
      title = "Clonotypes - no.filt")
  
  print(ggplot(SeuratObj@meta.data)+
          geom_point(aes(x= umap_1,
                         y= umap_2,
                         colour= hash.ID),
                     size = 0.6)+
          labs(title = "General conditions",
               colour = "Conditions")+
          theme_bw()+
          scale_colour_manual(values = MyPalette)+
          theme(legend.position = "right",
                panel.grid = element_blank()))
  
  print(ggplot(SeuratObj@meta.data)+
          geom_point(aes(x= umap_1,
                         y= umap_2,
                         colour= Phase),
                     size = 0.6)+
          labs(title = "Cell Cycle",
               colour = "CC Phases")+
          theme_bw()+
          scale_colour_manual(values = MyPalette)+
          theme(legend.position = "right",
                panel.grid = element_blank()))
  
  print(ggplot(SeuratObj@meta.data)+
          geom_point(aes(x= umap_1,
                         y= umap_2,
                         colour= clone_id),
                     size = 0.6)+
          labs(title = "Clonotypes")+
          theme_bw()+
          theme(legend.position = "none",
                panel.grid = element_blank()))
  
  print(ggplot(SeuratObj@meta.data)+
          geom_point(aes(x= umap_1,
                         y= umap_2,
                         colour= clone_id),
                     size = 0.6)+
          labs(title = "TCR")+
          theme_bw()+
          theme(legend.position = "none",
                panel.grid = element_blank()))
  
  dev.off()
  
  
  ##  Save seuratobject       ------------------------------------------------
  saveRDS(SeuratObj, file = paste0(Seurat_outdir,"Seurat2.rds"))
}



rm(hto, Rawdata, umi, VDJdata, barcodes, Implement_VDJ, Lib)

# cowplot::plot_grid(p1,p2,p3,p4, ncol = 2, labels = letters[1:4], label_size = 16)




