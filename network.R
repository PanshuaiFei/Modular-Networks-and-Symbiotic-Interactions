setwd("C:/Users/feips/Desktop/新挑选微生物组22个/新挑选微生物组22个/data-for-run/network-cons")
library(igraph)
rm(list=ls())
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
write_graph(t_net_its,"control_network.graphml","graphml")




# 假设你已经有了图 't_net_its' 和分类学数据 'tax_its'
# 计算度和模块性类别
degrees <- degree(t_net_its)
modularity_classes <- membership(cluster_louvain(t_net_its))

# 创建包含节点信息的数据框
node_info <- data.frame(
  Id = V(t_net_its)$name,
  modularity_class = modularity_classes,
  degree = degrees,
  row.names = V(t_net_its)$name
)

# 添加丰度数据 (A1...B16) 到 node_info 数据框
abundance_data <- read.delim('total.tsv', row.names=1, sep='\t', check.names=FALSE)
node_info <- merge(node_info, abundance_data, by="row.names", all.x=TRUE)
# 添加分类学信息从 'tax_its' 到 node_info 数据框
# 请确保 'tax_its' 的 row.names 是节点的ID
node_info <- merge(node_info, tax_its, by.x="Id", by.y="row.names", all.x=TRUE)
# 获取 tax_its 数据框的列名
tax_columns <- colnames(tax_its)
# 现在 tax_columns 包含了正确的列名
# 使用 tax_columns 来替换 node_info 中的 NA 值
node_info[,tax_columns] <- lapply(node_info[,tax_columns], function(x) ifelse(is.na(x), "others", x))
# 写入节点文件
write.csv(node_info, file="nodes_info.csv", row.names=FALSE)#########################
r_value<-read.delim('median_correlation.tsv',row.names=1,sep='\t',check.names=FALSE) #提取相关系数
p_value<-read.delim('pvalues.tsv',row.names=1,sep='\t',check.names=FALSE) #提取p值
#筛选p值小于0.05的相关分析结果
r_value[p_value>0.05|abs(r_value)<0.6] = 0
#将r_value保存为csv文件，gephi可视化
write.csv(r_value,file="total_network.csv")
#gephi导入，计算模块性进行上色，导出nodes文件，合并nodes文件，OTU丰度文件。
#导入nodes文件
nodes<-read.csv("nodes_info.csv",row.names = 1)
#导出主要的4个模块


modu1<-subset(nodes,nodes$modularity_class=="1")
modu2<-subset(nodes,nodes$modularity_class=="2")
modu3<-subset(nodes,nodes$modularity_class=="5")
modu4<-subset(nodes,nodes$modularity_class=="7")
#计算每个模块标准化的平均丰度
zmodu1<-t(scale(t(modu1[,4:179])))
abundance1<-as.data.frame(colSums(zmodu1)/296)
zmodu2<-t(scale(t(modu2[,4:179])))
abundance2<-as.data.frame(colSums(zmodu2)/230)
zmodu3<-t(scale(t(modu3[,4:179])))
abundance3<-as.data.frame(colSums(zmodu3)/72)
zmodu4<-t(scale(t(modu4[,4:179])))
abundance4<-as.data.frame(colSums(zmodu4)/33)
#合并每个模块的平均丰度
df<-data.frame(abundance1$`colSums(zmodu1)/296`,abundance2$`colSums(zmodu2)/230`,abundance3$`colSums(zmodu3)/72`,abundance4$`colSums(zmodu4)/33`)
rownames(df)<-rownames(abundance1)
colnames(df)<-c("modu1","modu2","modu3","modu4")
#将模块丰度和固氮功能进行相关分析
# 't_net_its' is your igraph network object
# Make sure it contains the 'modularity_class' as a vertex attribute

# If 'modularity_class' is not yet part of your igraph object's vertex attributes, add it as follows:
V(t_net_its)$modularity_class <- modularity_classes # modularity_classes obtained from your previous analysis

# Now, write the graph to a GraphML file
write_graph(t_net_its, file="network_with_modules.graphml", format="graphml")
#将模块丰度和固氮功能进行相关分析
library(ggplot2)
library(ggpubr)
design<-read.csv("treatment.csv",row.names = 1,check.names = F)

aa<-cbind(df,design)

#计算模块1
ggplot(aa,aes(x=aa$modu1,y=aa$N.fixation))+
  geom_point(size=4,aes(color=aa$Treatment))+
  scale_color_manual(values=c("#0055AA","#C40003","#00C19B","#EAC862"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
  theme_classic() +
  labs(x="Abundance in module 1",y="N fixation",color="Treatment")+
  theme(axis.text=element_text(colour='black',size=9))
View(aa)
#计算模块1
ggplot(aa,aes(x=aa$modu1,y=aa$Treatment))+
  geom_point(size=4,aes(color=aa$Treatment))+
  scale_color_manual(values=c("#0055AA","#C40003","#00C19B","#EAC862"))+
  geom_smooth(method=lm,level=0.95,color="gray4")+
  stat_cor(method = "pearson",label.x.npc ="left",label.y.npc = 0.02)+
  theme_classic() +
  labs(x="Abundance in module 1",y="N fixation",color="Treatment")+
  theme(axis.text=element_text(colour='black',size=9))




##展示每个模块的主要微生物组成
#合并四个模块的微生物群落
microbe<-rbind(modu1,modu2)
microbe<-rbind(microbe,modu3)
microbe<-rbind(microbe,modu4)
#统计丰度组成
total<-aggregate(microbe[,4:179],by=list(microbe$Phylum,microbe$modularity_class),FUN=sum)
abundance<-data.frame(rowSums(total[,3:178])/176)
io<-cbind(total[,1:2],abundance)
colnames(io)<-c("taxa","modu","abundance")
#转化为因子型
io$modu<-as.factor(io$modu)
#每个模块的物种组成
ggplot(data=io,aes(x=modu,y=abundance,fill=taxa))+
  geom_bar(position = "fill",stat = "identity",width = 0.6)+
  ylab("Relative abundance")+
  scale_fill_manual(values=c("#0055AA","#7FD2FF","#007ED3","#C40003","#00C19B","#EAC862","#B2DF8A","#FFACAA","#FF9D1E","#C3EF00","#CAB2D6","#894FC6","skyblue","red","#df9dc0","#5c7ada","#312d50","#63599d","#86547a"))+  scale_y_reverse(expand = c(0,0),labels  = c("1","0.75","0.50","0.25","0"))+
  theme(axis.text.x = element_text(color = "black",size = 8,angle = 90))+
  theme(axis.text.y = element_text(color = "black",size = 8))+
  theme(legend.position = "right",legend.text = element_text(size = 7),
        panel.grid =element_blank())+scale_x_discrete(name='Module' #x轴坐标名称
        )+
  guides(fill=guide_legend(title="Phylum",color="black",reverse=TRUE))+theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=9))

