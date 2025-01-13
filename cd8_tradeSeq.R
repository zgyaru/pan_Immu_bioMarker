library(Seurat)
library(tradeSeq)
library(slingshot)

data = readRDS('./data/scRNA/Tcell_atls/CD8.rds')

num.cells <- rowSums(data@assays$RNA@counts > 0)
genes.use <- names(num.cells[which(num.cells >= 100)])

data = data[genes.use,]

counts = data@assays$RNA@counts

counts = counts[c('CD8A','PDCD1','TCF7','SLAMF6','CXCL13','SPRY1','HAVCR2','LAG3'),]

curves = readRDS('./scRNA/Tcell_atls/slingshot_CD8_curves.rds')

pseudotime <- slingPseudotime(curves, na = FALSE)
cellWeights <- slingCurveWeights(curves)

gamList <- fitGAM(counts = counts,
                  pseudotime = pseudotime,
                  cellWeights = cellWeights,
                  nknots = 6)

saveRDS(gamList, './data/scRNA/Tcell_atls/tradeSeq_CD8_smalGenes.rds')


