options(stringsAsFactors = F)
## Figure 2

## Fig2-abcde 
# similar to Fig1-bcdef

## Fig2-f
plot_cells(mycds, 
           color_cells_by = 'cell_type',
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size=0.5,group_label_size=4)+
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=16))

## Fig2-g
plot_cells(mycds1, label_cell_groups = F, 
           color_cells_by = "pseudotime", 
           label_leaves = F, 
           label_branch_points = F, 
           graph_label_size = 2, 
           cell_size=0.5, 
           trajectory_graph_segment_size = 2)+
  theme(axis.text = element_text(size=28),
        axis.title = element_text(size=32))

## Fig2-h
gene_sce <- "IGHG4"
plot_cells(mycds1,
           genes=gene_sce,           
           label_cell_groups=F,
           show_trajectory_graph=T, 
           cell_size=1, trajectory_graph_color = "black", 
           label_branch_points = F, 
           label_roots = F, label_leaves = F)+
  scale_color_viridis(option="inferno")+
  theme(text = element_text(size = 14),
        legend.title =  element_text(size = 16),
        axis.title = element_text(size = 16)
  )

## Fig2-i
ggplot(plotdf, aes(x=pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8))+
  scale_x_continuous(position = 'top')

## Fig2-j
library(Hmisc)
library(dplyr)
library(scales)
library(dplyr)
library(future.apply)
library(viridis)
library(ComplexHeatmap)
library(clusterProfiler)
# library(org.Mm.eg.db)
library(org.Hs.eg.db)
fit_m_genes <- c("MS4A1",'HLA−DPB1','CD37','CD79B','MKI67','IGHD','TOP2A','IL4R','IGHG1',
                 'IGHG3','IGHV4−39','IGHG4','IGHG2','IGHV1−18','IGHM')
fit_m_genes <- as.data.frame(fit_m_genes)
pdf('pseudotime_hetamap.pdf', width=4, height=8)
ht_list <-Heatmap(module, col = module_col,
                  cluster_rows = FALSE,
                  show_row_names = FALSE,
                  show_column_names = FALSE,
                  row_title = NULL,
                  column_title = NULL,
                  heatmap_legend_param = list(title = "Module"))+
  Heatmap(
    as.matrix(ordered),
    show_column_names = FALSE, 
    show_row_names=FALSE,
    col = viridis(256),
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    use_raster = TRUE,
    heatmap_legend_param = list(title = "%Max")
  )+
  rowAnnotation(link = anno_mark(at = which(rownames(ordered) %in% fit_m_genes$fit_m_genes), 
                                 labels = fit_m_genes$fit_m_genes, labels_gp = gpar(fontsize = 8)))
draw(ht_list, row_km = 3, row_split = module,
     column_title = "Heatmap of pseudotime DEGs", 
     column_title_gp = gpar(fontsize = 12, fontface = "bold"), 
     merge_legends = TRUE, heatmap_legend_side = "right") 
dev.off()

## Fig2-k
library(cowplot)
library(ggpubr)
library(dplyr)
library(readr)
library(tidyr)
library(ggforce)
library(pals)
library(pheatmap)
library(scales)
library(ggthemes)
library(Seurat)
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(Hmisc)
ggplot(GSEA_df2)+
geom_col(aes(x = reorder(pathway, ES),y = ES, fill = group))+
scale_fill_manual(values = c("#2f73bb","#ae4531"))+
geom_segment(aes(y = 0, yend = 0,x = 0, xend = length))+
theme_classic()+
ylim(-1,1)+
coord_flip()+
theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
  )+
  ylab("Enrichment Score")+
  geom_text(data=GSEA_df2[which(GSEA_df2$ES > 0), ],aes(x = pathway, y = 0, label = pathway), 
            hjust = 1, size = 3)+
  geom_text(data=GSEA_df2[which(GSEA_df2$ES < 0), ],aes(x = pathway, y = 0, label = pathway), 
            hjust = -0, size = 4)+
  ggtitle(paste0(cell_type," Subtype GO-BP-ENRICH \n P < 0.1"))+
  scale_x_discrete(expand=expansion(add=c(0,1.5)))+
  geom_segment(aes(y = -0, yend = -0.8,x = length, xend = length),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
  geom_segment(aes(y = 0, yend = 0.8,x = length, xend = length),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), 
               size = 0.5)+
  annotate("text", x = length, y = -1, label = "down")+
  annotate("text", x = length, y = 1, label = "up")