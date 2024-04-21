## Figure 3

## Fig3-ab
# similar to Fig1-bd

## Fig3-c-1
# similar to Fig2-h

## Fig3-c-2
gene <- "IGHG4"
GENE.EXP = scRNA@assays$RNA@data[which(rownames(scRNA) ==gene),]
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT$VALUE=GENE.EXP
OUT=vector.showValue(OUT)

## Fig3-d
library('circlize')
library('gatepoints')
library('stringr')
library('igraph')
library('gmodels')
VEC = scRNA@reductions$umap@cell.embeddings
rownames(VEC) = colnames(scRNA)
PCA =scRNA@reductions$pca@cell.embeddings
source('../Vector.R')
PCA=vector.rankPCA(PCA)
OUT=vector.buildGrid(VEC, N=30,SHOW=TRUE)
OUT=vector.buildNet(OUT, CUT=1, SHOW=TRUE)
OUT=vector.getValue(OUT, PCA, SHOW=TRUE)
OUT=vector.gridValue(OUT,SHOW=TRUE)
OUT=vector.autoCenter(OUT,UP=0.9,SHOW=TRUE)
OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)
vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)

## Fig3-ef
# similar to Fig2-ij

## Fig3-hij
# similar to Fig2-k

