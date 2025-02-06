library(Seurat)
library(harmony)
library(dplyr)
library(limma)
library(R.utils)
library(scDblFinder)
library(BiocParallel)
library(ggplot2)

library(tidyverse)
library(CellChat)
setwd("D:/fx/GBMLGG/GSE182109")
seu_obj<-readRDS("seu_obj_celltype.RDS")
dim(seu_obj)
# 30210 122902
############下面进行cellchat分析
folders<-list.files("D:/fx/GBMLGG/DiscoveryOutput_LGGGBM_sample_cutoff0.99_p0.05")
folders<-folders[1:12]
folders<-folders[folders!="Ecotypes"]
data_all<-data.frame()

####读入所有被分配了细胞状态的数据
for (folder in folders) {
  file_path<-file.path("D:/fx/GBMLGG/DiscoveryOutput_LGGGBM_sample_cutoff0.99_p0.05",folder,"state_assignment.txt")
  data_state<-read.table(file_path,header = T,as.is = T,sep = '\t')
  row.names(data_state)<-NULL
  data_state$celltype<-folder
  data_all<-rbind(data_all,data_state)
}
data_all<-data_all[,-3]

##合并表型信息
metadata<-seu_obj@meta.data
metadata$ID<-row.names(metadata)
metadata<-merge(metadata,data_all)
metadata<-metadata%>%column_to_rownames(var="ID")
row.names(metadata)
seu_obj<-AddMetaData(seu_obj,metadata = metadata)
seu_obj_filter<-subset(seu_obj,cells=row.names(metadata))

##########根据提取到的细胞进行cellchat分析
seu_obj_filter$eco<-paste0(seu_obj_filter$celltype,"_",seu_obj_filter$State)
E4<-c("Astrocyte_S06","Fibroblast_S02","Macrophage_S02","Microglia_S03","Neuron_S03","NKT_S02","NPC_S04","Oligodendrocyte_S03","OPC_S01")
seu_obj_E4<-subset(seu_obj_filter,eco%in%E4)
#########创建cellchat对象
cellchat <- createCellChat(object = seu_obj_E4, meta = seu_obj_E4@meta.data, group.by = "eco")
cellchat@meta$samples<-cellchat@meta$orig.ident
####指定数据库
CellChatDB <- CellChatDB.human 
cellchat@DB <- CellChatDB   ####细胞通讯的数据库包含多个，如果不想提取子数据库直接用全部的就好

cellchat <- subsetData(cellchat) 
#识别过表达基因
cellchat <- identifyOverExpressedGenes(cellchat)
#识别过表达配体受体对
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)###这一步可以不做
cellchat@data.smooth[1:5,1:5]

###推断细胞通讯网络
cellchat <- computeCommunProb(cellchat, raw.use = F, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "./cellchat/net_lr.csv")


#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
df.netp<-subsetCommunication(cellchat,slot.name = "netP")
write.csv(df.netp,"./cellchatnet_pathway.csv")
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


cellchat@netP$pathways
#Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = "SPP1", layout = "chord")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = "SPP1", layout = "circle")

par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "SPP1", color.heatmap = "Reds")


netAnalysis_contribution(cellchat, signaling = "SPP1")

levels(cellchat@idents)   
vertex.receiver = c(1,2,3,5,6,7,8,9) 
#par(mar=c(5.1,4.1,4.1,2.1))
netVisual_aggregate(cellchat, signaling ="SPP1",  vertex.receiver = vertex.receiver,layout = "hierarchy")
# save as TIL/CXCL_hierarchy.pdf



netVisual_bubble(cellchat, sources.use = 3, targets.use = c(5:11), remove.isolate = FALSE)
cellchat@net$count


#######
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5,width = 6,height = 25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5,width = 6,height = 25)
ht1 + ht2
p = plotGeneExpression(cellchat, signaling = "SPP1")



###############################################              E2
E2<-c("Neuron_S01","Astrocyte_S02","Macrophage_S03","NKT_S03","Fibroblast_S07")
seu_obj_E2<-subset(seu_obj_filter,eco%in%E2)


#########创建cellchat对象
cellchat_E2 <- createCellChat(object = seu_obj_E2, meta = seu_obj_E2@meta.data, group.by = "eco")
cellchat_E2@meta$samples<-cellchat_E2@meta$orig.ident
####指定数据库
CellchatDB <- CellChatDB.human 
cellchat_E2@DB <- CellChatDB   ####细胞通讯的数据库包含多个，如果不想提取子数据库直接用全部的就好

cellchat_E2 <- subsetData(cellchat_E2) 
#识别过表达基因
cellchat_E2 <- identifyOverExpressedGenes(cellchat_E2)
#识别过表达配体受体对
cellchat_E2 <- identifyOverExpressedInteractions(cellchat_E2)
cellchat_E2 <- smoothData(cellchat_E2, adj = PPI.human)###这一步可以不做
cellchat_E2@data.smooth[1:5,1:5]

###推断细胞通讯网络
cellchat_E2 <- computeCommunProb(cellchat_E2, raw.use = F, population.size = TRUE) 
cellchat_E2 <- filterCommunication(cellchat_E2, min.cells = 10)
df.net_E2 <- subsetCommunication(cellchat_E2)
write.csv(df.net_E2, "D:/fx/GBMLGG/GSE182109/cellchat/net_lr_E2.csv")


#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat_E2 <- computeCommunProbPathway(cellchat_E2)
df.netp_E2<-subsetCommunication(cellchat_E2,slot.name = "netP")
write.csv(df.netp_E2,"D:/fx/GBMLGG/GSE182109/cellchat/cellchat_E2net_pathway.csv")
#计算整合的细胞类型之间通信结果
cellchat_E2 <- aggregateNet(cellchat_E2)

groupSize <- as.numeric(table(cellchat_E2@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_E2@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_E2@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")



cellchat_E2<- netAnalysis_computeCentrality(cellchat_E2, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat_E2, pattern = "outgoing", font.size = 5,width = 6,height = 25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_E2, pattern = "incoming", font.size = 5,width = 6,height = 25)
ht1 + ht2

par(mfrow=c(1,1))
netVisual_heatmap(cellchat_E2, signaling = "COLLAGEN", color.heatmap = "Reds")
