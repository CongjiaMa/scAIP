## Figure 4

## Fig4-abc
# similar to Fig1-bdc

## Fig4-d
# similar to Fig2-g

## Fig4-e
# similar to Fig2-i

## Fig4-f
# similar to Fig2-j

## Fig4-g
genes_sig <- c('NKG7','CD8A','CXCR4','ETS1','S100A6','LEF1')
P1 <- plot_genes_in_pseudotime(mycds1[genes_sig,], color_cells_by="cell_type", 
                               min_expr=0.5, ncol = 1,
                               cell_size = 2)+
  scale_color_manual(values = mycol)+
  theme(text = element_text(size = 14),
        legend.title =  element_text(size = 16),
        axis.title = element_text(size = 16)
  )
P2 <- plot_genes_in_pseudotime(mycds1[genes_sig,], color_cells_by="pseudotime", 
                               min_expr=0.5, ncol = 1)+
  theme(text = element_text(size = 14),
        legend.title =  element_text(size = 16),
        axis.title = element_text(size = 16)
  )
p <- P1+P2

## Fig4-hi
# similar to Fig2-k