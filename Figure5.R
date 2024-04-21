## Figure 5

## Fig5-abc
# similar to Fig1-bdc

## Fig5-d
# similar to Fig2-g

## Fig5-e
# similar to Fig1-e

## Fig5-f
# similar to Fig2-g

## Fig5-g
# similar to Fig2-i

## Fig5-h
library(Seurat)
library(pheatmap)
library(tidyr)
annotation_col = data.frame(
  celltype = c('Naive_CD8','Cytotoxic_CD8','Naive_CD4',
               'NAMPT_inflammatory_CD4','TH_CD4','Tfh','NK'),
  type = c(rep("CD8",2),rep("CD4",4),rep("NK",1)))
row.names(annotation_col) <- colnames(exp)
annotation_row = data.frame(
  celltype = c('Naive_CD8','Cytotoxic_CD8','Naive_CD4',
               'NAMPT_inflammatory_CD4','TH_CD4','Tfh','NK'),
  type = c(rep("CD8",2),rep("CD4",4),rep("NK",1))
)
row.names(annotation_row) <- rownames(exp)
allcolour <-  c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                         "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                         "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                         "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                         "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
mycol = allcolour[1:length(table(sce$cell_type))]      
ann_colors = list(celltype = c(Naive_CD8="#DC050C",Cytotoxic_CD8="#FB8072" ,
                               Naive_CD4="#1965B0" ,NAMPT_inflammatory_CD4="#7BAFDE" ,
                               TH_CD4="#882E72" ,Tfh= "#B17BA6",NK=  "#FF7F00"),
                               type = c(CD8="#612FDE",CD4="#397E92",NK="#00bfff")) 
library(psych)
library(pheatmap)
library(ggplot2)
library(gtable)
library(grid)
pheatmap::pheatmap(exp,annotation_col = annotation_col,annotation_row = annotation_row, 
                   color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),
                   cluster_rows = T,cluster_cols = T,annotation_colors = ann_colors,
                   fontsize_row = 15,fontsize_col = 15,fontsize = 11)
                        
## Fig5-i
# similar to Fig2-j

## Fig5-j
library(ggrepel)
library(ggrastr)
library(dplyr)
score <- AddModuleScore_UCell(uterus,
                              features=features,
                              name="_score")
ggplot(data2, aes(x=cellid,y=m2_score,fill=cellid,color=cellid)) +
theme_bw()+RotatedAxis()+
theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=12),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
labs(x=NULL,y=NULL,title = "M2_score")+ 
geom_jitter_rast(col="#00000033", pch=19,cex=2, position = position_jitter(0.2))+
geom_boxplot(position=position_dodge(0))+
scale_fill_manual(values = colors)+
geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
geom_text(data=cell_number2,aes(cellid,y,label=number), color='black') 
