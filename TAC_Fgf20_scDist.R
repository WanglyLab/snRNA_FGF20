library(Seurat)
library(scDist)
library(dplyr)

rm(list=ls())
load('./seurat_combined.RData')

my_meta <- seurat_obj@meta.data %>%
  select(Group, Sample, NewCellType) %>%
  filter(Group %in% c('TAC', 'TAC+SLNP@Fgf20'))

my_norm <- GetAssayData(seurat_obj, assay = 'RNA', layer = 'data')
my_norm <- my_norm[, rownames(my_meta)]
my_meta <- my_meta[colnames(my_norm), ]

TAC_Fgf20_diff <- scDist(
  normalized_counts = my_norm,
  meta.data = my_meta,
  fixed.effects = "Group",
  random.effects = "Sample",
  d = 20,
  clusters = "NewCellType"
)

save(TAC_Fgf20_diff, file = 'TAC_Fgf20_diff_scDist.RData')
DistPlot(TAC_Fgf20_diff)
