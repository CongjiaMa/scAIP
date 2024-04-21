## Figure 6

## Fig6-a
library(CellChat)
library(ggraph)
library(tidygraph)
library(dplyr)
ggraph(data,layout='igraph', algorithm = 'circle') +
  geom_edge_bend(mapping = aes(edge_width = inter_num),
                 strength = 0.2,alpha = 0.5,
                 flipped =T, edge_color = "#A9AAAA",
                 n=50, show.legend = F,
                 check_overlap =T)+
  geom_edge_loop(aes(edge_width = inter_num,
                     direction = (from - 1)*180 / length(data2)),
                 colour = "#A9AAAA",
                 alpha = 0.5, show.legend = F)+
  scale_edge_width_continuous(range = c(0,5))+
  geom_node_point(aes(size=inter_num,colour = inter_num)) +
  geom_node_point(aes(size=inter_num), show.legend = F,
                  shape=21,colour = 'black',stroke = 1.5)+
  geom_node_text(aes(x = x*1.06, y=y*1.06, label=name),
                 angle=0,hjust=0, fontface="bold",size=3) + 
  scale_color_gradientn(colors = colorRampPalette(c("#2166AC",'#478ABF','#90C0DC', "white",'#EF8C65','#CF4F45',"#B2182B"))(100))+
  scale_size_continuous(range = c(1, 15))+
  theme_graph()

## Fig6-b
library(patchwork)
pathways.show='VISFATIN'
netVisual_aggregate(group1_cellchat, signaling = pathways.show, layout = "circle",targets.use = c('Stellate_cell'))+
  title(main = "VISFATIN pathway \n Stellate_cell in Normal")