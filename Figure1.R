options(stringsAsFactors = F)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(magrittr)
library(colorfindr)
library(dittoSeq)
library(ggrastr)
library(tidydr)
library(scatterpie)
library(tidyr)
library(DESeq2)
library(combinat)
library(clusterProfiler)
## Figure 1

## Fig1b-1
meta  = seuratObj@meta.data
meta$UMAP_1 <- seuratObj@reductions$umap@cell.embeddings[,1] 
meta$UMAP_2 <- seuratObj@reductions$umap@cell.embeddings[,2] 
col_df = data.frame(name=unique(meta$cell_type)) %>% 
  arrange(name)
mycol <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
mycol <- mycol[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,
                 2,4,6,8,10,12,14,16,18,20,22,24,26,28,30)]
mycol =setNames(mycol,col_df$name)
mid_coord_type = meta %>% 
dplyr::select(c('cell_type','UMAP_1','UMAP_2')) %>% 
group_by(cell_type) %>% 
summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
axis <- ggh4x::guide_axis_truncated(
trunc_lower = unit(0, "npc"),
trunc_upper = unit(3, "cm")
)
ggplot(meta,aes( UMAP_1 ,UMAP_2))+
geom_point(size=0.01,aes(color=cell_type))+
geom_text(data = mid_coord_type,size=5,
aes(UMAP_1,UMAP_2,label=cell_type)) +
theme_classic()+
theme(legend.title = element_text(size=24),
legend.text = element_text(size=21),
axis.title = element_text(size=24),
axis.text = element_text(size=21))+
guides(color=guide_legend(override.aes = list(size=4)))+
scale_color_manual(values = mycol)
## Fig1b-2
uterus <- sce
df <- uterus@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(cell_type = uterus@meta.data$cell_type)
Idents(uterus)<-"cell_type"
Cellratio <- prop.table(table(uterus$type,Idents(uterus)), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
freq <-spread(Cellratio, Var1, Freq)
colnames(freq)[1] <- 'celltype'
freq <- freq[sort(freq$celltype),]
label <- df %>%group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1),
            UMAP_2 = median(UMAP_2))%>%
  as.data.frame()
rownames(label) <- label$cell_type
cell_number <- as.data.frame(table(uterus$cell_type))
colnames(cell_number)[2]<-'cellnumber'
cell_number$cellnumber <- log2(cell_number$cellnumber)/9
data = cbind(freq,label[,c(2:3)], cell_number[,c(2)])
colnames(data)[10]<- 'cellnumber'
ggplot()+
  geom_point_rast(data=df, aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type),size = 2,shape=16) +
  scale_color_manual(values = alpha(dittoColors(),0.5))+ 
  theme_dr()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 21))+
  geom_scatterpie(data=data,
                  aes(x=UMAP_1,y=UMAP_2,
                      group=celltype,
                      r=cellnumber),
                  cols=names(freq)[2:7])+
  scale_fill_manual(values = c("#409079","#52a5c1","#c65341","#425785","#d6873b","#92b8da"),name='group')

## Fig1c
library(fmsb)
pdf('./youpath/radarchart.pdf', width=5, height=5)
radarchart(my.data, 
           pty = c(16,16,32),
           axistype = 1,
           pcol = c("#64299C", "#0439FD","black"), 
           pfcol = c(scales::alpha(c("#64299C", "#0439FD","black"), c(0.5, 0.5,0.5))),
           plwd = c(3,3,3),
           plty = 1,
           cglcol = "grey60", 
           cglty = 1, 
           cglwd = 1,
           axislabcol = "grey60",
           vlcex = 1.2, 
           vlabels = colnames(colnames(my.data)),
           caxislabels = c(0, 1000,2000,3000,4000),
           calcex=1.5)
legend(x = "bottomright", legend = c("blood","tissue"), horiz = F,
       bty = "n", pch = 15 , col = c("#64299C", "#0439FD"),
       text.col = "black", cex = 1, pt.cex = 3)
legend(x = "topright", legend = c("AIP-1"), horiz = TRUE,
       text.col = "black", cex = 3,bg=NULL,box.lty=0)
dev.off()

#Fig1d
library(scRNAtoolVis)
library(patchwork)
pdf1<-cellRatioPlot(object = sce1,
                    sample.name = "type",
                    celltype.name = "cell_type")+
  theme(axis.text.x=element_text(angle=45, hjust=1,size = 14),
        axis.text.y=element_text(size = 14),
        axis.title.y=element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14),
        legend.position = "none")
pdf2<-cellRatioPlot(object = sce2,
                    sample.name = "type",
                    celltype.name = "cell_type")+
  theme(axis.text.x=element_text(angle=45, hjust=1,size = 14),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14))
pdf <- pdf1+pdf2

#Fig1-ef
library(miloR)
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)
library(SeuratWrappers)
library(ggbeeswarm)
library(scater)
library(scales)
library(forcats)
library(data.table)
library(stringr)
library(dplyr)
plotNhoodGraphDA(sce1, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",
                       mid="lightgrey",
                       high="#910000",
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish)
da_results <- annotateNhoods(sce1, da_results, coldata_col = "cell_type")
plotDAbeeswarm(da_results, group.by = "cell_type") +
  scale_color_gradient2(low="#070091",
                        mid="lightgrey",
                        high="#910000",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black'))