library(Seurat)
library(scDist)
library(dplyr)                
library(miloR)                
library(SingleCellExperiment) 

rm(list=ls())
setwd('/data01/common/001_LNP/001_processed/')

dat1 <- Read10X(data.dir = './TAC_1/output/filter_matrix/',gene.column = 1)
dat2 <- Read10X(data.dir = './TAC_2/output/filter_matrix/',gene.column = 1)
dat3 <- Read10X(data.dir = './TAC_3/output/filter_matrix/',gene.column = 1)

dat4 <- Read10X(data.dir = './Sham_1/output/filter_matrix/',gene.column = 1)
dat5 <- Read10X(data.dir = './Sham_2/output/filter_matrix/',gene.column = 1)
dat6 <- Read10X(data.dir = './Sham_3/output/filter_matrix/',gene.column = 1)

dat7 <- Read10X(data.dir = './TAC+SLNP@Fgf20_1/output/filter_matrix/',gene.column = 1)
dat8 <- Read10X(data.dir = './TAC+SLNP@Fgf20_1/output/filter_matrix/',gene.column = 1)
dat9 <- Read10X(data.dir = './TAC+SLNP@Fgf20_3/output/filter_matrix/',gene.column = 1)

colnames(dat1) <- paste0('TAC_1_',colnames(dat1))
colnames(dat2) <- paste0('TAC_2_',colnames(dat2))
colnames(dat3) <- paste0('TAC_3_',colnames(dat3))

colnames(dat4) <- paste0('Sham_1_',colnames(dat4))
colnames(dat5) <- paste0('Sham_2_',colnames(dat5))
colnames(dat6) <- paste0('Sham_3_',colnames(dat6))

colnames(dat7) <- paste0('TAC_1_',colnames(dat7))
colnames(dat8) <- paste0('TAC_2_',colnames(dat8))
colnames(dat9) <- paste0('TAC_3_',colnames(dat9))

#--get intersect gene_name---
gene_list <- list(rownames(dat1),rownames(dat2),rownames(dat3),
                  rownames(dat4),rownames(dat5),rownames(dat6),
                  rownames(dat7),rownames(dat8),rownames(dat9))
interG <- Reduce(intersect, gene_list)
length(interG)

combined_dat <- cbind(dat1[interG,],dat2[interG,],dat3[interG,],
                      dat4[interG,],dat5[interG,],dat6[interG,],
                      dat7[interG,],dat8[interG,],dat9[interG,])


identical(rownames(meta_dat),colnames(combined_dat))
seurat_obj <- CreateSeuratObject(counts = combined_dat,meta.data = meta_dat)
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")  
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 500 & nCount_RNA > 500 & percent.mt < 15)
#---processing---
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2500)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
ElbowPlot(seurat_obj,ndims = 50)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:40)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:40)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.6)
save(seurat_obj, file = 'seurat_combined.RData')