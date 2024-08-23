setwd("E:/Academic/SparCC/Gut_microbe_Protocol/新挑选微生物组22个/新挑选微生物组22个/data-for-run/total")
rm(list = ls())
library(igraph)
#读取分类信息
#####输入OTU物种分类数据 #####（提前把相对丰度信息加入到分类信息中）
tax_its=read.csv("micro_taxa.csv",row.names=1)
colnames(tax_its)=c("Kindom","Phylum","Class","Order","Family","Genus","Species")
tax_its[tax_its==""]="unassigned"

#读取fastspar计算的correlation 和p矩阵
cor_sparcc=read.delim('median_correlation.tsv',row.names=1,sep='\t',check.names=FALSE)
cor_sparcc[abs(cor_sparcc)<0.6]=0

pvals=read.delim('pvalues.tsv',row.names=1,sep='\t',check.names=FALSE)

#筛选p值
pvals[pvals>=0.05]= -1    # p值大于0.05的统一为-1
pvals[pvals>=0]= 1         #大于0的为1
pvals[pvals==-1]=0       #等于-1的为0
#两个矩阵相乘，获得共同matrix
adj=as.matrix(cor_sparcc)*as.matrix(pvals)
#输出矩阵
write.table(data.frame(adj,check.names=FALSE),'network.adj.txt',col.names=NA,sep='\t',quote=FALSE)
#读取网络矩阵
network_adj=read.delim('network.adj.txt',row.names=1,sep='\t',check.names=FALSE)
#（注意network.adj.txt中最后一列，如有FALSE，则删除network.adj.txt中最后一列FALSE）
#获取网络信息
g=graph_from_adjacency_matrix(as.matrix(network_adj),mode='undirected',weighted=TRUE,diag=FALSE)
#给边赋值weight，因为存在负值
E(g)$sparcc=E(g)$weight
E(g)$weight=abs(E(g)$weight)
E(g)$pattern=ifelse(E(g)$sparcc>0,"positive","negative")    #添加正负权重
edge=data.frame(as_edgelist(g))

#获得边信息
edge_list=data.frame(from=edge[[1]],to=edge[[2]],weight=E(g)$weight,sparcc=E(g)$sparcc, pattern=E(g)$pattern)
head(edge_list)
#获得点信息
nodeattrib_i_its=data.frame(node=union(edge_list$from,edge_list$to))
nodeattrib_i_its

#整合节点分类信息
nodeattrib_t_its=cbind(nodeattrib_i_its,tax_its[as.character(nodeattrib_i_its$node),])
t_net_its=graph_from_data_frame(edge_list,direct=FALSE, vertices = nodeattrib_t_its)


#输出网络：
write_graph(t_net_its,"66666666control_network.graphml","graphml")




occor.p=pvals
occor.r=cor_sparcc
#根据上述筛选的 r 值和 p 值保留数据
z <- occor.r * occor.p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
z[abs(z)>0]=1
z
adjacency_unweight <- z

####################################################################################
#计算节点度
V(g)$degree <- degree(g)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(g)$modularity <- membership(cluster_fast_greedy(g))

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
  nodes_id = V(g)$name, 
  degree = V(g)$degree, 
  modularity = V(g)$modularity
)
head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块

write.table(nodes_list, 'nodes_list.txt', sep = '\t', row.names = FALSE, quote = FALSE)






##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')#这里要用小白鱼大佬的代码，将代码放在自己R的工作目录下就行了，后台回复Zi-Pi可获得


#上述的邻接矩阵类型的网络文件
adjacency_unweight 

#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('nodes_list.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

write.table(zi_pi, 'zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
write.csv(zi_pi,"zipi结果.csv")
ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.6, size = 12,shape=17) +
  scale_y_continuous(limits=c(-3,3))+
  scale_color_manual(values = c("#443B84FF","#20A486FF","#FDE725FF", "#2C738EFF"), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.61,linetype=2,linewidth=1) +
  geom_hline(yintercept = 2.5,linetype=2,linewidth=1)+  
  theme_bw()+
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )







# Your ggplot code
p <- ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.6, size = 12, shape=17) +
  scale_y_continuous(limits=c(-3, 3)) +
  scale_color_manual(values = c("#443B84FF","#20A486FF","#FDE725FF", "#2C738EFF"), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs')) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.61, linetype=2, linewidth=1) +
  geom_hline(yintercept = 2.5, linetype=2, linewidth=1) +  
  theme_bw() +
  theme(axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")) )

# Save the plot to a PDF file
ggsave("plot.pdf", plot = p, width = 15, height = 14, device = 'pdf')
