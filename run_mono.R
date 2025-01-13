library(Seurat)
library(monocle)

data = readRDS('./data/scRNA/Tcell_atls/CD8.rds')
curves = readRDS('./data/scRNA/Tcell_atls/slingshot_CD8_curves.rds')
pt <- slingPseudotime(curves)

num.cells <- rowSums(data@assays$RNA@counts > 0)
genes.use <- names(num.cells[which(num.cells >= 100)])
pseudotime = pt[, 'Lineage2']
pseudotime = na.omit(pseudotime)
sub_data = data[genes.use,
                names(pseudotime)]
sub_data$Pseudotime = as.numeric(pseudotime)

pd <- new("AnnotatedDataFrame", data = sub_data@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame('gene' = rownames(sub_data), row.names = rownames(sub_data)))
mon <- newCellDataSet(sub_data@assays$RNA@counts,
    phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())
mon <- estimateSizeFactors(mon)mon <- estimateSizeFactors(mon)
mon <- estimateDispersions(mon)
saveRDS(mon,'./data/scRNA/Tcell_atls/monocle_CD8_lineage.rds')

