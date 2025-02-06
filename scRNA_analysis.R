library(Seurat)
library(harmony)
library(dplyr)
library(limma)
library(R.utils)
library(scDblFinder)
library(BiocParallel)
library(ggplot2)
library(scRNAtoolVis)
library(COSG)
setwd("D:/fx/GBMLGG/GSE182109")

samples <- list.files("./data_lgggbm")

# 创建一个空的列表来存储Seurat对象
seurat_list <- list()

# 读取每个样本的10x数据并创建Seurat对象
for (sample in samples) {
  # 拼接文件路径
  data.path <- paste0("./data_lgggbm/", sample)
  
  # 读取10x数据，data.dir参数指定存放文件的路径
  seurat_data <- Read10X(data.dir = data.path)
  
  # 创建Seurat对象，并指定项目名称为样本文件名
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)
  
  # 将Seurat对象添加到列表中
  seurat_list <- append(seurat_list, seurat_obj)
}


# 合并Seurat对象，将所有Seurat对象合并到一个对象中
seu_obj <- merge(seurat_list[[1]], 
                          y = seurat_list[-1],
                          add.cell.ids = samples)
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")

VlnPlot(seu_obj,features = "nFeature_RNA",group.by = "orig.ident")
VlnPlot(seu_obj,features = "nCount_RNA",group.by = "orig.ident")
VlnPlot(seu_obj,features = "pMT",group.by = "orig.ident")

dim(seu_obj)
nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 500
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 10

seu_obj <- subset(seu_obj, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & pMT < pMT_upper )

table(seu_obj$orig.ident)

dim(seu_obj)

seu_obj <- NormalizeData(seu_obj)
#高变基因
seu_obj <- FindVariableFeatures(seu_obj, selection.method = "vst", nfeatures = 2000)
#归一化
seu_obj <- ScaleData(seu_obj)
seu_obj<-RunPCA(seu_obj, npcs = 50)
seu_obj= seu_obj %>% RunHarmony("orig.ident", plot_convergence = F)

table(seu_obj$orig.ident)
#标准流程，参数不变
seu_obj <- seu_obj %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
seu_obj <- FindClusters(seu_obj ,resolution = 0.8)
DimPlot(seu_obj, reduction = "tsne", pt.size=0.5, label = T,repel = TRUE,group.by = "RNA_snn_res.0.8",raster = F)+theme(legend.position = "none")
DimPlot(seu_obj, reduction = "tsne", pt.size=0.5, label = T,repel = TRUE,group.by = "RNA_snn_res.0.8",raster = F)
DimPlot(seu_obj, reduction = "tsne", pt.size=0.5, label = F,repel = TRUE,group.by = "orig.ident",raster = F)
DimPlot(seu_obj, reduction = "tsne", pt.size=0.5, label = F,repel = TRUE,group.by = "orig.ident",raster = F)+theme(legend.position = "none")
DimPlot(seu_obj, reduction = "umap", pt.size=0.5, label = T,repel = TRUE,group.by = "RNA_snn_res.0.8",raster = F)
DimPlot(seu_obj, reduction = "umap", pt.size=0.5, label = F,repel = TRUE,group.by = "orig.ident",raster = F)
saveRDS(seu_obj,"GBMLGG_qc_data.rds")

########细胞注释
Idents(seu_obj)<-"RNA_snn_res.0.8"
cluster_0.8_markers<-FindAllMarkers(object = seu_obj, 
               only.pos = TRUE,
               min.pct = 0.25, 
               thresh.use = 0.25)
###saveRDS(cluster_0.8_markers,"cluster_0.8_markers.rds")
cluster_0.8_markers$filter<-substring(cluster_0.8_markers$gene,1,3)
cluster_0.8_markers<-subset(cluster_0.8_markers,subset = filter!="LOC")
cluster_0.8_markers<-subset(cluster_0.8_markers,subset = (p_val_adj<0.05))
cluster_0.8_markers<-subset(cluster_0.8_markers,subset = (avg_log2FC>1))
 
#####首先按照文章里面的进行注释
#######髓系细胞
VlnPlot(seu_obj,features = c("PTPRC","ITGAM","CD68"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)


#######gloima cell
VlnPlot(seu_obj,features = c("SOX2","OLIG1","GFAP","S100B"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)


#######t cell
VlnPlot(seu_obj,features = c("CD3D","CD3E","CD4","CD8A"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)


#####B cell
VlnPlot(seu_obj,features = c("CD79A","CD19"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)

#####Endothelial
VlnPlot(seu_obj,features = c("PECAM1","VWF"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)


####Fibrobalsts
VlnPlot(seu_obj,features = c("DCN","COL1A1"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)

table(seu_obj$RNA_snn_res.0.8)
#########Oligodendrocytes
VlnPlot(seu_obj,features = c("PLP1","TF","PPP1R14A"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)



VlnPlot(seu_obj,features = c("PDGFRA","OLIG1"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)

VlnPlot(seu_obj,features = c("BCAN","SCRG1"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)
VlnPlot(seu_obj,features = c("MAP2","STMN2"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)

VlnPlot(seu_obj,features = c("PCSK1N"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)
###NPC
VlnPlot(seu_obj,features = c("HIST1H4C","PCLAF"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)

###

##Astrocyte
VlnPlot(seu_obj,features = c("AQP4","GFAP"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)


####macrophage
VlnPlot(seu_obj,features = c("CD163","CD68"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)



#########Microglial
VlnPlot(seu_obj,features = c("P2RY12","CX3CR1","TMEM119"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)


VlnPlot(seu_obj,features = c("HMGB2","MKI67"),group.by ="RNA_snn_res.0.8",raster=F,pt.size = 0)

seu_obj = RenameIdents(seu_obj,"0"="Microglia",
                            "1"="Microglia",
                            "2"="OPC",
                            "3"="Astrocyte",
                            "4"="Macrophage",
                            "5"="T_cell",
                            "6"="NPC",
                            "7"="Astrocyte",
                            "8"="NPC",
                            "9"="Fibroblast",
                            "10"="Microglia",
                            "11"="Astrocyte",
                            "12"="Astrocyte",
                            "13"="T_cell",
                            "14"="Oligodendrocyte",
                            "15"="Microglia",
                            "16"="Neuron",
                            "17"="Endothelium",
                            "18"="NPC",
                            "19"="Astrocyte",
                            "20"="B_cell",
                            "21"="OPC",
                            "22"="Astrocyte",
                            "23"="Oligodendrocyte",
                            "24"="Astrocyte",
                            "25"="Astrocyte",
                            "26"="Astrocyte")
#保存细胞身份鉴定结果至Seurat对象的meta.data中
seu_obj$cell_type = Idents(seu_obj)
#可视化每个样品中的细胞占比情况
seu_obj@meta.data %>%
  ggplot(aes(x=orig.ident,fill=cell_type))+
  geom_bar(position = position_fill())+
  scale_fill_manual(values=c("#A2C8D9","#ACD486","#E49B9B","#F6BE73","#CCB5D3",
                             "#FED703","#FCE1D4","#B25A28","#D9D9D9","#99973B"
  ,"#D75B4E"))+theme_classic()+theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic"))+
  coord_flip() 
#保存最终数据
saveRDS(seu_obj, "seu_obj_celltype.RDS")

#### 6、可视化 ####
#降维图,展示细胞在低维空间的分布情况,group.by参数可以设置多个样品分组
DimPlot(seu_obj, label = TRUE, group.by="cell_type",reduction = "umap",repel = T,pt.size = 0.5,raster = F)

#####注释完发现有一部分T细胞和小胶质细胞发生了污染
T_cell<-subset(seu_obj,cell_type=="T cell")
T_cell<-T_cell@assays$RNA@counts
T_cell<-CreateSeuratObject(counts = T_cell)
T_cell <- NormalizeData(T_cell)
#高变基因
T_cell <- FindVariableFeatures(T_cell, selection.method = "vst", nfeatures = 2000)
#归一化
T_cell <- ScaleData(T_cell)
T_cell<-RunPCA(T_cell, npcs = 50)
T_cell= T_cell %>% RunHarmony("orig.ident", plot_convergence = F)

table(T_cell$orig.ident)
#标准流程，参数不变
T_cell <- T_cell %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  RunTSNE(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.2) %>% 
  identity()
DimPlot(T_cell, reduction = "tsne", pt.size=0.5, label = T,repel = TRUE,group.by = "RNA_snn_res.0.2",raster = F)

T_markers<-FindAllMarkers(object = T_cell, 
                          only.pos = TRUE,
                          min.pct = 0.25, 
                          thresh.use = 0.25)
T_markers$filter<-substring(T_markers$gene,1,3)
T_markers<-subset(T_markers,subset = filter!="LOC")
T_markers<-subset(T_markers,subset = (p_val_adj<0.05))
T_markers<-subset(T_markers,subset = (avg_log2FC>1))
VlnPlot(T_cell,features = c("CD3D","CD3E"),group.by ="RNA_snn_res.0.2",raster=F,pt.size = 0)
VlnPlot(T_cell,features = c("P2RY12","CX3CR1","TMEM119"),group.by ="RNA_snn_res.0.2",raster=F,pt.size = 0)
Idents(T_cell)<-"RNA_snn_res.0.2"
table(T_cell$RNA_snn_res.0.2)
T_cell = RenameIdents(T_cell,"0"="NKT",
                       "1"="Microglia",
                       "2"="NKT",
                       "3"="NKT",
                       "4"="NKT",
                       "5"="NKT",
                       "6"="NKT",
                       "7"="NKT")

#保存细胞身份鉴定结果至Seurat对象的meta.data中
T_cell$cell_type2 = Idents(T_cell)
metadata<-seu_obj@meta.data
metadata1<-T_cell@meta.data
metadata$cellname<-row.names(metadata)
metadata1$cellname<-row.names(metadata1)
metadata1<-metadata1[,-(1:5)]
metadata2<-merge(metadata,metadata1,by.x="cellname",by.y="cellname",all.x=T)
metadata2 <- metadata2 %>% 
  mutate(cell_type2 = coalesce(cell_type2, cell_type))
library(tidyverse)
metadata2<-metadata2%>%
  column_to_rownames(var="cellname")
seu_obj<-AddMetaData(seu_obj,metadata = metadata2)
DimPlot(T_cell, label = TRUE, group.by="cell_type2",reduction = "tsne",repel = T,pt.size = 0.5,raster = F)
DimPlot(seu_obj, label = F, group.by="cell_type2",reduction = "tsne",repel = T,pt.size = 0.5,raster = F)

saveRDS(seu_obj, "seu_obj_celltype.RDS")

seu_obj<-readRDS("D:/fx/GBMLGG/GSE182109/seu_obj_celltype.RDS")

######出现了一点问题
####需要将细胞名中的-换成.

colnames(seu_obj)<-gsub("-",".",colnames(seu_obj))


library(data.table)
#######接下来输出每种细胞类型提取500个细胞单细胞矩阵以及注释文件
protein_gene<-fread("D:/fx/protein_coding.tsv")
Idents(seu_obj)<-"cell_type2"
seu_obj_sample<-subset(seu_obj,downsample=10000)
table(seu_obj_sample$cell_type2)
metadata_sample<-seu_obj_sample@meta.data
metadata_sample<-metadata_sample[,c("orig.ident","cell_type2")]

colnames(metadata_sample)<-c("Sample","CellType")
celltype_data<-GetAssayData(seu_obj_sample, slot = "data")
celltype_data<-data.frame(Gene=rownames(celltype_data),celltype_data)
row.names(celltype_data)<-NULL
celltype_data<-celltype_data[celltype_data$Gene %in% protein_gene$name,]

metadata_sample$CellType<-as.character(metadata_sample$CellType)
metadata_sample$CellType[metadata_sample$CellType=="B cell"]<-"B_cell"
table(metadata_sample$CellType)
fwrite(data.frame(ID=rownames(metadata_sample),metadata_sample),"GBMLGG_sample_Annotation.txt", row.names=F,sep = "\t")
fwrite(celltype_data,"GBMLGG_sample.txt", row.names=F,sep = "\t")
table(metadata_sample$CellType)








#########绘制注释气泡图
levels(seu_obj$cell_type2)
celltype_markers<-c("CD3D","CD3E","P2RY12","CX3CR1","TMEM119","BCAN","SCRG1","AQP4","GFAP","CD163","CD68","CD3D","CD3E","HIST1H4C","PCLAF","DCN","COL1A1","PLP1","TF","PPP1R14A","STMN2","SOX4","PECAM1","VWF","CD79A","CD19")
jjDotPlot(object = seu_obj,
          gene = celltype_markers,
          xtree = F,
          rescale = T,
          rescale.min = 0,
          rescale.max = 1,
          point.geom = F,
          tile.geom = T,ytree=F,id="cell_type2")
Idents(seu_obj)<-seu_obj$cell_type2
averageHeatmap(object = seu_obj,
               markerGene = celltype_markers)

