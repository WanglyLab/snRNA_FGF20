library(Seurat)
library(scDist)
library(dplyr)
library(miloR)
library(SingleCellExperiment)

rm(list=ls())
load('./seurat_combined.RData')

seurat_sub <- subset(seurat_obj, subset = Group %in% c('TAC', 'TAC+SLNP@Fgf20'))
sce <- as.SingleCellExperiment(seurat_sub)

milo <- Milo(sce)

milo <- buildGraph(milo, k = 40, d = 40, reduced.dim = "PCA")
milo <- makeNhoods(milo, prop = 0.05, k = 40, d = 40, refined = TRUE)

milo <- countCells(
  milo,
  meta.data = as.data.frame(colData(milo)),
  sample = "Sample"
)

milo_design <- distinct(as.data.frame(colData(milo))[, c("Sample", "Group")])
rownames(milo_design) <- milo_design$Sample

design_mat <- model.matrix(~ Group, data = milo_design)

milo <- calcNhoodDistance(milo, d = 40, reduced.dim = "PCA")

da_results <- testNhoods(
  milo,
  design = design_mat,
  design.df = milo_design
)

da_results <- annotateNhoods(milo, da_results, coldata_col = "NewCellType")

milo <- buildNhoodGraph(milo)

save(da_results, file = 'miloR_NewCellType.RData')
