####博士单细胞代码 质控部分导出为step1.1.rdata####
###安包
library(patchwork)
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(Augur)
version(Seurat)
update.packages()
##样本C1质控
dir='sample/'
samples <- list.files(dir)
C1.data <- Read10X(data.dir = "sample/C1K/")
colnames(C1.data) <- paste0("C1_", colnames(C1.data))
ncol(C1.data)
C1 <- CreateSeuratObject(counts = C1.data,
                         project = "C1",
                         min.cells = 3,
                         min.features = 200)
#计算线粒体比例
C1$mito.percent <- PercentageFeatureSet(object = C1, pattern = "^mt-")
C1$mito.percent <- C1@meta.data$mito.percent / 100
fivenum(C1@meta.data$mito.percent)#最小值，下四分位数，中位数，上四分位数，最大值
View(C1@meta.data)
#计算核糖体基因比例
ribo.percent=rownames(C1)[grep("^Rp[sl]", rownames(C1),ignore.case = T)]
ribo.percent
C1=PercentageFeatureSet(C1, "^Rp[sl]", col.name = "ribo.percent")
C1$ribo.percent <- C1@meta.data$ribo.percent / 100
fivenum(C1@meta.data$ribo.percent)
fivenum(C1@meta.data$nCount_RNA)
#画质控图
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "mito.percent","ribo.percent"), ncol = 4)
C1 <- subset(C1, subset = nFeature_RNA < 5000 & nCount_RNA < 9827 & mito.percent < 0.2 & ribo.percent < 0.183)
#标准化
C1 <- NormalizeData(C1, verbose = FALSE)
plot1 <- FeatureScatter(C1, feature1 = "nCount_RNA", feature2 = "mito.percent")
plot2 <- FeatureScatter(C1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#找高变基因
C1 <- FindVariableFeatures(C1, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
C1_top10 <- head(VariableFeatures(C1), 10)
C1_top10
plot1 <- VariableFeaturePlot(C1)
plot2 <- LabelPoints(plot = plot1, points = C1_top10)
CombinePlots(plots = list(plot1, plot2))


##样本C2质控
#计算线粒体比例
samples <- list.files(dir)
C2 <- Read10X(data.dir = "sample/C2T286/")
colnames(C2) <- paste0("C2_", colnames(C2))
ncol(C2)
C2 <- CreateSeuratObject(counts = C2,
                         project = "C2",
                         min.cells = 3,
                         min.features = 200)
C2$mito.percent <- PercentageFeatureSet(object = C2, pattern = "^mt-")
C2$mito.percent <- C2@meta.data$mito.percent / 100
fivenum(C2@meta.data$mito.percent)#最小值，下四分位数，中位数，上四分位数，最大值
head(C2@meta.data)
View(C2@meta.data)
#计算核糖体基因比例
ribo.percent=rownames(C2)[grep("^Rp[sl]", rownames(C2),ignore.case = T)]
ribo.percent
C2=PercentageFeatureSet(C2, "^Rp[sl]", col.name = "ribo.percent")
C2$ribo.percent <- C2@meta.data$ribo.percent / 100
fivenum(C2@meta.data$ribo.percent)
fivenum(C2@meta.data$nCount_RNA)
#画质控图
VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "mito.percent","ribo.percent"), ncol = 4)
C2 <- subset(C2, subset = nFeature_RNA < 5000 & nCount_RNA < 9421 & mito.percent < 0.2 & ribo.percent < 0.224)
#标准化
C2 <- NormalizeData(C2, verbose = FALSE)
plot1 <- FeatureScatter(C2, feature1 = "nCount_RNA", feature2 = "mito.percent")
plot2 <- FeatureScatter(C2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#找高变基因
C2 <- FindVariableFeatures(C2, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
C2_top10 <- head(VariableFeatures(C2), 10)
C2_top10
plot1 <- VariableFeaturePlot(C2)
plot2 <- LabelPoints(plot = plot1, points = C2_top10)
CombinePlots(plots = list(plot1, plot2))

##样本C3质控
C3 <- Read10X(data.dir = "sample/C3T279/")
colnames(C3) <- paste0("C3_", colnames(C3))
ncol(C3)
C3 <- CreateSeuratObject(counts = C3,
                         project = "C3",
                         min.cells = 3,
                         min.features = 200)
#计算线粒体比例
C3$mito.percent <- PercentageFeatureSet(object = C3, pattern = "^mt-")
C3$mito.percent <- C3@meta.data$mito.percent / 100
fivenum(C3@meta.data$mito.percent)#最小值，下四分位数，中位数，上四分位数，最大值
head(C3@meta.data)
View(C3@meta.data)
#计算核糖体基因比例
ribo.percent=rownames(C3)[grep("^Rp[sl]", rownames(C3),ignore.case = T)]
ribo.percent
C3=PercentageFeatureSet(C3, "^Rp[sl]", col.name = "ribo.percent")
C3$ribo.percent <- C3@meta.data$ribo.percent / 100
fivenum(C3@meta.data$ribo.percent)
fivenum(C3@meta.data$nCount_RNA)
#画质控图
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA", "mito.percent","ribo.percent"), ncol = 4)
#C3 <- subset(C3, subset = nFeature_RNA < 5000 & nCount_RNA < 10312 & mito.percent < 0.2 & ribo.percent < 0.209)
C3 <- subset(C3, subset = nFeature_RNA < 5000 & nCount_RNA < 10521 & mito.percent < 0.2 & ribo.percent < 0.209)
#标准化
C3 <- NormalizeData(C3, verbose = FALSE)
plot1 <- FeatureScatter(C3, feature1 = "nCount_RNA", feature2 = "mito.percent")
plot2 <- FeatureScatter(C3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#找高变基因
C3 <- FindVariableFeatures(C3, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
C3_top10 <- head(VariableFeatures(C3), 10)
C3_top10
plot1 <- VariableFeaturePlot(C3)
plot2 <- LabelPoints(plot = plot1, points = C3_top10)
CombinePlots(plots = list(plot1, plot2))



##样本Pb1质控
Pb1 <- Read10X(data.dir = "sample/Pb1T289/")
colnames(Pb1) <- paste0("Pb1_", colnames(Pb1))
ncol(Pb1)
Pb1 <- CreateSeuratObject(counts = Pb1,
                          project = "Pb1",
                          min.cells = 3,
                          min.features = 200)
Pb1$mito.percent <- PercentageFeatureSet(object = Pb1, pattern = "^mt-")
Pb1$mito.percent <- Pb1@meta.data$mito.percent / 100
fivenum(Pb1@meta.data$mito.percent)#最小值，下四分位数，中位数，上四分位数，最大值
head(Pb1@meta.data)
View(Pb1@meta.data)
#计算核糖体基因比例
ribo.percent=rownames(Pb1)[grep("^Rp[sl]", rownames(Pb1),ignore.case = T)]
ribo.percent
Pb1=PercentageFeatureSet(Pb1, "^Rp[sl]", col.name = "ribo.percent")
Pb1$ribo.percent <- Pb1@meta.data$ribo.percent / 100
fivenum(Pb1@meta.data$ribo.percent)
fivenum(Pb1@meta.data$nCount_RNA)
#画质控图
VlnPlot(Pb1, features = c("nFeature_RNA", "nCount_RNA", "mito.percent","ribo.percent"), ncol = 4)
#Pb1 <- subset(Pb1, subset = nFeature_RNA < 5000 & nCount_RNA < 10521 & mito.percent < 0.2 & ribo.percent < 0.209)
Pb1 <- subset(Pb1, subset = nFeature_RNA < 5000 & nCount_RNA < 10312 & mito.percent < 0.2 & ribo.percent < 0.209)
#标准化
Pb1 <- NormalizeData(Pb1, verbose = FALSE)
plot1 <- FeatureScatter(Pb1, feature1 = "nCount_RNA", feature2 = "mito.percent")
plot2 <- FeatureScatter(Pb1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#找高变基因
Pb1 <- FindVariableFeatures(Pb1, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
Pb1_top10 <- head(VariableFeatures(Pb1), 10)
Pb1_top10
plot1 <- VariableFeaturePlot(Pb1)
plot2 <- LabelPoints(plot = plot1, points = Pb1_top10)
CombinePlots(plots = list(plot1, plot2))


##样本Pb2质控
Pb2 <- Read10X(data.dir = "sample/Pb2T282/")
colnames(Pb2) <- paste0("Pb2_", colnames(Pb2))
ncol(Pb2)
Pb2 <- CreateSeuratObject(counts = Pb2,
                          project = "Pb2",
                          min.cells = 3,
                          min.features = 200)
View(Pb2@meta.data)
Pb2$mito.percent <- PercentageFeatureSet(object = Pb2, pattern = "^mt-")
Pb2$mito.percent <- Pb2@meta.data$mito.percent / 100
fivenum(Pb2@meta.data$mito.percent)#最小值，下四分位数，中位数，上四分位数，最大值
head(Pb2@meta.data)
View(Pb2@meta.data)
#计算核糖体基因比例
ribo.percent=rownames(Pb2)[grep("^Rp[sl]", rownames(Pb2),ignore.case = T)]
ribo.percent
Pb2=PercentageFeatureSet(Pb2, "^Rp[sl]", col.name = "ribo.percent")
Pb2$ribo.percent <- Pb2@meta.data$ribo.percent / 100
fivenum(Pb2@meta.data$ribo.percent)
fivenum(Pb2@meta.data$nCount_RNA)
#画质控图
VlnPlot(Pb2, features = c("nFeature_RNA", "nCount_RNA", "mito.percent","ribo.percent"), ncol = 4)
Pb2 <- subset(Pb2, subset = nFeature_RNA < 5000 & nCount_RNA < 8870 & mito.percent < 0.2 & ribo.percent < 0.206)
#标准化
Pb2 <- NormalizeData(Pb2, verbose = FALSE)
plot1 <- FeatureScatter(Pb2, feature1 = "nCount_RNA", feature2 = "mito.percent")
plot2 <- FeatureScatter(Pb2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#找高变基因
Pb2 <- FindVariableFeatures(Pb2, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
Pb2_top10 <- head(VariableFeatures(Pb2), 10)
Pb2_top10
plot1 <- VariableFeaturePlot(Pb2)
plot2 <- LabelPoints(plot = plot1, points = Pb2_top10)
CombinePlots(plots = list(plot1, plot2))


##样本Pb3质控
Pb3 <- Read10X(data.dir = "sample/Pb3323/")
colnames(Pb3) <- paste0("Pb3_", colnames(Pb3))
ncol(Pb3)
Pb3 <- CreateSeuratObject(counts = Pb3,
                          project = "Pb3",
                          min.cells = 3,
                          min.features = 200)
View(Pb3@meta.data)
Pb3$mito.percent <- PercentageFeatureSet(object = Pb3, pattern = "^mt-")
Pb3$mito.percent <- Pb3@meta.data$mito.percent / 100
fivenum(Pb3@meta.data$mito.percent)#最小值，下四分位数，中位数，上四分位数，最大值
head(Pb3@meta.data)
View(Pb3@meta.data)
#计算核糖体基因比例
ribo.percent=rownames(Pb3)[grep("^Rp[sl]", rownames(Pb3),ignore.case = T)]
ribo.percent
Pb3=PercentageFeatureSet(Pb3, "^Rp[sl]", col.name = "ribo.percent")
Pb3$ribo.percent <- Pb3@meta.data$ribo.percent / 100
fivenum(Pb3@meta.data$ribo.percent)
fivenum(Pb3@meta.data$nCount_RNA)
#画质控图
VlnPlot(Pb3, features = c("nFeature_RNA", "nCount_RNA", "mito.percent","ribo.percent"), ncol = 4)
Pb3 <- subset(Pb3, subset = nFeature_RNA < 5000 & nCount_RNA < 12966 & mito.percent < 0.2 & ribo.percent < 0.199)
#标准化
Pb3 <- NormalizeData(Pb3, verbose = FALSE)
plot1 <- FeatureScatter(Pb3, feature1 = "nCount_RNA", feature2 = "mito.percent")
plot2 <- FeatureScatter(Pb3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#找高变基因
Pb3 <- FindVariableFeatures(Pb3, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
Pb3_top10 <- head(VariableFeatures(Pb3), 10)
Pb3_top10
plot1 <- VariableFeaturePlot(Pb3)
plot2 <- LabelPoints(plot = plot1, points = Pb3_top10)
CombinePlots(plots = list(plot1, plot2))
#保存数据
save(Pb1,Pb2,Pb3, C1, C2, C3, file = "step1.1.rdata")

####博士单细胞代码 数据合并及降为去批次和不去批次版 导出为dm2qpc.Rdata####
##包的安装
library(patchwork)
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(rJava)
library(devtools)
library(usethis)
library(rJava)
library(harmony)
library(celldex)
library(SingleR)
library(matrixStats)
library(IRanges)
library(devtools)
library(DelayedArray)
library(clustree)
library(tidyverse)
library(RColorBrewer)
library(DO.db)
##数据分组及降为 保存为

View(C1@meta.data)
scelist <- list(C1,C2,C3,Pb1,Pb2,Pb3)
sce.all <- merge(scelist[[1]],
                 y= scelist[ -1 ]) 

View(sce.all@meta.data)
dim(sce.all)
metadata <- sce.all@meta.data
metadata$cells <- rownames(metadata)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^C"))] <- "Control"
metadata$sample[which(str_detect(metadata$cells, "^P"))] <- "Pb"
sce.all@meta.data <- metadata
sce.all <- sce.all %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
plot1 <- VariableFeaturePlot(sce.all)
plot2 <- LabelPoints(plot = plot1, points = Pb2_top10)
CombinePlots(plots = list(plot1, plot2))
##runpca不去批次效应
sce.all <- RunPCA(sce.all, npcs = 30, verbose = FALSE)
sce.all <- FindNeighbors(sce.all, dims = 1:20)
sce.all <-FindClusters(sce.all, resolution = 0.3)
sce.all <- RunUMAP(sce.all,reduction = "pca",dims = 1:20)
DimPlot(sce.all,reduction = "umap",pt.size = 0.3 )
DimPlot(sce.all,reduction = "umap",group.by="sample",pt.size = 0.3) 
##harmony去批次效应
sce.all <- RunPCA(sce.all, npcs = 30, verbose = FALSE)
sce.all <- RunHarmony(sce.all, "orig.ident")
names(sce.all@reductions)
sce.all <- RunUMAP(sce.all,dims = 1:20, 
                   reduction = "harmony")
sce.all <- RunTSNE(sce.all,  dims = 1:20,
                   reduction = "harmony")

sce.all <- FindNeighbors(sce.all, reduction = "harmony", dims = 1:20)
sce.all <-FindClusters(sce.all, resolution = 0.3)

DimPlot(sce.all,reduction = "umap",group.by="sample",pt.size = 0.2 ) 
DimPlot(sce.all,reduction = "umap",pt.size = 0.3 ) 

DimPlot(sce.all,reduction = "tsne",group.by="sample",pt.size = 0.2 ) 
DimPlot(sce.all,reduction = "tsne",pt.size = 0.3 ) 
DimPlot(sce.all, reduction = "pca")
#clustree找最合适的分辨率点
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  sce.all <- FindClusters(sce.all, resolution = res)
}
p1 <- clustree(sce.all, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(sce.all, group.by = 'RNA_snn_res.0.3', label = T)
p1 + p2 + plot_layout(widths = c(3, 1))
#看一下三维及标准化后的数据
DimHeatmap(sce.all, dims = 1:20,cell = 500)
ElbowPlot(sce.all, ndims = 30)
View(sce.all[["RNA"]]@scale.data[1:30,1:30])
#cell cycle看一下是否有细胞周期效应
sce.all<- CellCycleScoring(object = sce.all, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = sce.all@meta.data)
DimPlot(sce.all,reduction = "umap",label = F,group.by="Phase",pt.size = 0.2)
ggsave('figure/cellcycle.tiff',height = 7,width = 11)
##找亚群基因
sce.all.markers <- FindAllMarkers(sce.all, 
                                  only.pos = TRUE,  # 只返回positive基因
                                  min.pct = 0.25)
sce.top.markers = sce.all.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC) 
write.table(sce.all.markers,file="allmarker1.txt",sep="\t")
write.table(sce.top.markers,file="topmarker1.txt",sep="\t")
save(metadata,sce.all,sce.all.markers,file = "dm2qpc.Rdata")
####亚群注释
####细胞分群注释 bs_totalcell.Rdata####
###安包加载数据

load("dm2qpc.Rdata")
load("bs_code/bs_totalcell.Rdata")
##细胞注释
View(sce.all@meta.data)
cell.used <- rownames(sce.all@meta.data[which(!(sce.all@meta.data$seurat_clusters%in% c(8))),])
sce.all <- subset(sce.all, cells = cell.used)
DimPlot(sce.all, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
DimPlot(sce.all,reduction = "tSNE",pt.size = 0.3,label = T) 
DimPlot(sce.all,reduction = "umap",pt.size = 0.3,label = T,split.by = "sample") 
levels(sce.all)
Idents(sce.all)
#小胶质细胞(0,3,7)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Cx3cr1","Hexb","Tmem119","C1qa"))
#星形胶质细胞(1,2,5,11)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Atp1b2","Aldoc","Gfap"))
#Ependymal cell(13)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Tmem212"))
#少突胶质细胞(4,16,17,6,15)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Mbp","Apod","Plp1"))
#内皮细胞(12)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Cldn5"))
#上皮细胞(6,15)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Ttr"))
#血管平滑肌细胞(10)SMC
FeaturePlot(sce.all, pt.size = 0.5,features = c("Acta2"))
#神经元(9)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Dcx","Sox11"))
#Neutrophils (18)
FeaturePlot(sce.all, pt.size = 0.5,features = c("S100a8"))
# neural stem cell (14)
FeaturePlot(sce.all, pt.size = 0.5,features = c("Ube2c"))
p <- FeaturePlot(sce.all, pt.size = 0.5,features = c("Hexb","Aldoc","Mbp","Ttr","Tmem212","Acta2",
                                                     "Ube2c","Sox11","Vtn","Cldn5","S100a8"),ncol = 3)
p
ggsave('figure/bsall_marker.tiff',height = 11,width = 11)
##marker重新定群
celltype=data.frame(ClusterID=0:18 ,
                    celltype= 0:18) 
#名字和群匹配  
celltype[celltype$ClusterID %in% c(0,3,7),2]='Microglia' 
celltype[celltype$ClusterID %in% c(2,5,1,11),2]='Astrocyte'  
celltype[celltype$ClusterID %in% c(4,16,17),2]='Oligodendrocyte'  
celltype[celltype$ClusterID %in% c(12),2]='Endothelial cells'  
celltype[celltype$ClusterID %in% c(6,15),2]='Epithelial cells'   
celltype[celltype$ClusterID %in% c(9),2]='Neuron'  
celltype[celltype$ClusterID %in% c(18),2]='Neutrophils'  
celltype[celltype$ClusterID %in% c(10),2]='SMC'  
celltype[celltype$ClusterID %in% c(14),2]='Neural stem cell' 
celltype[celltype$ClusterID %in% c(13),2]='Ependymal cell' 
celltype[celltype$ClusterID %in% c(8),2]='Mic_Ast' 
head(celltype)
celltype
table(celltype$celltype)
sce.all@meta.data$celltype = "NA"

for(i in 1:nrow(celltype)){
  sce.all@meta.data[which(sce.all@meta.data$RNA_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce.all@meta.data$celltype)
library(ggsci)
library(ggplot2)
p <- DimPlot(sce.all, reduction = "umap",group.by="celltype",label = TRUE, pt.size = 0.3) + NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Celltype")
plot4 = plot3 + scale_color_npg()

cell_type_cols <- c(brewer.pal(9, "Set2"), "#FF34B3", "#BC8F8F", "#20B2AA", "#00F5FF", 
                    "#FFA500", "#ADFF2F", "#FF6A6A", "#7FFFD4", "#AB82FF", "#90EE90", 
                    "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B",
                    "#FF3030", "#7CFC00", "#000000", "#708090")
DimPlot(sce.all, label = T, cols = cell_type_cols,  pt.size = 0.3, repel = T)
ggsave('figure/bsall_cell.tiff',height = 7,width = 11)
ggsave('figure/bsall_cell.pdf',height = 7,width = 11)
DimPlot(sce.all, group.by="sample",label = F, pt.size = 0.2)
ggsave('figure/bsall_cell_group.tiff',height = 7,width = 11)
#计算细胞比例
sce.all@active.ident <- plyr::mapvalues(x = sce.all@active.ident, 
                                        from = celltype$ClusterID, to = celltype$celltype)
sce.all@active.ident<-factor(sce.all@active.ident,levels = unique(celltype$celltype))
table(sce.all$sample)
prop.table(table(Idents(sce.all)))
table(Idents(sce.all), sce.all$sample)
Cellratio <- prop.table(table(Idents(sce.all), sce.all$sample), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Celltype","Group","Percentage")
colourCount = length(unique(Cellratio$Celltype))
ggplot(Cellratio) + 
  geom_bar(aes(x =Group, y= Percentage, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Percentage')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('bs_figure/radio_allcell.pdf',height = 4,width = 11)
save(celltype,sce.all,file = "bs_totalcell.Rdata")
#细胞比例画图
p1 = plot.clusters.group(data = sce.all,clusters =  "celltype", 
                         xlab = "Celltype", log = TRUE,
                         group = "orig.ident",legend.title = "Sample",
                         widths = c(3,1),color = 1)
p1 = plot.clusters.group(data = sce.all,clusters =  "celltype", 
                         xlab = "Celltype", log = T,
                         group = "sample",legend.title = "Sample",
                         widths = c(3,1),color = 1)
p1
####差异基因及富集分析####
##安包及加载数据
load("bs_code/bs_totalcell.Rdata")
load("bs_code/total_marker.Rdata")
library(dplyr)
library(Seurat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(msigdbr)
##总差异基因比较
meta.data <- sce.all@meta.data
all_marker<-FindMarkers(sce.all, group.by = "sample",
                        ident.1 = "Pb", 
                        ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.15)

all_marker = row.names(all_marker)
all_marker <- bitr(all_marker,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
all_marker = all_marker[,2]

go_total<- enrichGO(
    all_marker,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 3,
    maxGSSize = 20,
    readable = TRUE,
    pool = FALSE
  )
kegg_gene <- enrichKEGG(
  all_marker,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 3,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE)
write.csv(go_total,file = "table/go_total.csv")
write.csv(kegg_gene,file = "table/kegg_gene_total.csv")
##找差异基因
sce.all@meta.data$cell_group <- paste(sce.all@active.ident, sce.all@meta.data$sample, sep = "_")
mig_marker<-FindMarkers(sce.all, group.by="cell_group",
                        ident.1 = "Microglia_Pb", 
                        ident.2 = "Microglia_Control", logfc.threshold = 0.1, min.pct = 0.1)
ast_marker<-FindMarkers(sce.all, group.by="cell_group",
                        ident.1 = "Astrocyte_Pb", 
                        ident.2 = "Astrocyte_Control", logfc.threshold = 0.17, min.pct = 0.1)
Oli_marker<-FindMarkers(sce.all, group.by="cell_group",
                        ident.1 = "Oligodendrocyte_Pb", 
                        ident.2 = "Oligodendrocyte_Control", logfc.threshold = 0.15, min.pct = 0.1)
Neuron_marker <- FindMarkers(sce.all, group.by="cell_group",
                             ident.1 = "Neuron_Pb", 
                             ident.2 = "Neuron_Control", logfc.threshold = 0.22, min.pct = 0.1)
Epithelial_cells_marker <- FindMarkers(sce.all, group.by="cell_group",
                                      ident.1 = "Epithelial cells_Pb", 
                                      ident.2 = "Epithelial cells_Control", logfc.threshold = 0.17, min.pct = 0.1)
Neutrophils_marker <- FindMarkers(sce.all, group.by="cell_group",
                                  ident.1 = "Neutrophils_Pb", 
                                  ident.2 = "Neutrophils_Control", logfc.threshold = 0.7, min.pct = 0.1)
SMC_marker <- FindMarkers(sce.all, group.by="cell_group",
                          ident.1 = "SMC_Pb", 
                          ident.2 = "SMC_Control", logfc.threshold = 0.4, min.pct = 0.1)
Endothelial_cells_marker <- FindMarkers(sce.all, group.by="cell_group",
                                        ident.1 = "Endothelial cells_Pb", 
                                        ident.2 = "Endothelial cells_Control", 
                                        logfc.threshold = 0.5, min.pct = 0.1)
Ependymal_cell_marker<- FindMarkers(sce.all, group.by="cell_group",
                                    ident.1 = "Ependymal cell_Pb", 
                                    ident.2 = "Ependymal cell_Control", 
                                    logfc.threshold = 0.6, min.pct = 0.1)
Neural_stem_cell_marker <- FindMarkers(sce.all, group.by="cell_group",
                                       ident.1 = "Neural stem cell_Pb", 
                                       ident.2 = "Neural stem cell_Control", 
                                       logfc.threshold = 0.8, min.pct = 0.1)

write.csv(mig_marker,file = "bs_table/mig_marker.csv")
write.csv(ast_marker,file = "bs_table/ast_marker.csv")
write.csv(Oli_marker,file = "bs_table/Oli_marker.csv")
write.csv(Neuron_marker,file = "bs_table/Neuron_marker.csv")
write.csv(Epithelial_cells_marker,file = "bs_table/Epithelial_cells_marker.csv")
write.csv(Neutrophils_marker,file = "bs_table/Neutrophils_marker.csv")
write.csv(SMC_marker,file = "bs_table/SMC_marker.csv")
write.csv(Endothelial_cells_marker,file = "bs_table/Endothelial_cells_marker.csv")
write.csv(Ependymal_cell_marker,file = "bs_table/Ependymal_cell_marker.csv")
write.csv(Neural_stem_cell_marker,file = "bs_table/Neural_stem_cell_marker.csv")
#保存marker
save(ast_marker,Endothelial_cells_marker,Ependymal_cell_marker,
     Epithelial_cells_marker,mig_marker,Neural_stem_cell_marker,
     Neuron_marker,Neutrophils_marker,Oli_marker,SMC_marker,file = "bscode/total_marker.Rdata")

total <- list(mig=mig_marker,ast=ast_marker,Oli=Oli_marker,Neuron=Neuron_marker,
              Epi=Epithelial_cells_marker,End=Endothelial_cells_marker,SMC=SMC_marker,
              Neu=Neutrophils_marker,Epe=Ependymal_cell_marker,Ns=Endothelial_cells_marker)
i=1
for (i in 1:10) {
                  total[[i]]= row.names(total[[i]])
                  total[[i]] <- bitr(total[[i]],
                                   fromType = "SYMBOL",
                                   toType = "ENTREZID",
                                   OrgDb = org.Mm.eg.db)
                  total[[i]] = total[[i]][,2]
                  }
lapply(total, head)
#小胶质细胞差异基因作图
mig_marker<-FindMarkers(sce.all, group.by="cell_group",
                        ident.1 = "Microglia_Pb", 
                        ident.2 = "Microglia_Control", logfc.threshold = 0.1, min.pct = 0.1)
mig_marker$gene <- row.names(mig_marker)
logFC_t=0.20
P.Value_t = 0.05
k1 = (mig_marker$p_val_adj < P.Value_t)&(mig_marker$avg_log2FC < -logFC_t)
k2 = (mig_marker$p_val_adj < P.Value_t)&(mig_marker$avg_log2FC > logFC_t)
deg <- mutate(mig_marker,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)
dat  = deg[!duplicated(deg$gene),]
dat <- dat[-c(1,2),]
dat <- dat[-nrow(dat),]
dat
p <- ggplot(data = dat, 
            aes(x = avg_log2FC, 
                y = -log10(p_val_adj))) +
  geom_point(alpha=0.8, size=3, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("Blue", "Grey","Red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()+
  theme(panel.grid=element_blank())
#小胶质细胞差异富集分析

go_result <- list(mig_go=NA,ast_go=NA,Oli_go=NA,Neuron_go=NA,Epi_go=NA,End_go=NA,SMC_go=NA,
                  Neu_go=NA,Epe_go=NA,End_go=NA)
length(total)
i = 1
for (i in 1:length(total)){
  go_result[[i]]<- enrichGO(
  total[[i]],
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 3,
  maxGSSize = 20,
  readable = TRUE,
  pool = FALSE
  )
}

#批量导出go结果
filename<-c("mig_go","ast_go","Oli_go","Neuron_go","Epi_go","End_go","SMC_go",
      "Neu_go","Epe_go","End_go")

for(i in 1:length(go_result)){
  write.csv(go_result[[i]],file = paste0("bs_table/",filename[i],".csv"))
}

mig_gene_kegg <- enrichKEGG(
  mig_gene,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 3,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE
)
write.csv(mig_gene_kegg,file = "table/mig_gene_kegg.csv")
p <- do.call(cbind,total)
save(total,go_result, file = "bscode/total_go.Rdata")
df <- data.frame(matrix(unlist(total), nrow=length(total)))

####细胞通讯分析####
load("bs_code/bs_totalcell.Rdata")

library(CellChat)
library(patchwork)
library(Seurat)
###去掉混杂群
cell.used <- rownames(sce.all@meta.data[-which(sce.all$celltype %in% c('Mic_Ast',"Neural stem cell")),])
length(cell.used)
sce.all <- subset(sce.all, cells = cell.used)
DimPlot(sce.all, label = F, pt.size = 0.5)
###选取500个细胞
allCells=names(Idents(sce.all))
allType = levels(Idents(sce.all))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,500,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sce.all[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(cg_sce))
DimPlot(cg_sce, label = F, pt.size = 1)
alldata <- SplitObject(cg_sce,
                       split.by = "sample")
alldata
table(alldata[["Control"]]@meta.data$celltype)
table(alldata[["Pb"]]@meta.data$celltype)

##对照组通讯
object_C = alldata$Control
C.str <- GetAssayData(object = object_C,
                       assay = "RNA", slot = "data")

C.meta.data <- alldata[["Control"]]@meta.data
C.cellchat <- createCellChat(object = C.str,
                              meta = C.meta.data,
                              group.by = "celltype")
levels(C.cellchat@idents)
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
str(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB,
                           search = "Secreted Signaling")
C.cellchat@DB <- CellChatDB.use
future::plan("multiprocess", workers = 6)
C.cellchat <- subsetData(C.cellchat)
C.cellchat <- identifyOverExpressedGenes(C.cellchat)
C.cellchat <- identifyOverExpressedInteractions(C.cellchat)
C.cellchat <- projectData(C.cellchat, PPI.mouse)

C.cellchat <- computeCommunProb(C.cellchat, raw.use = TRUE)
C.cellchat <- filterCommunication(C.cellchat, min.cells = 10)
C.cellchat <- computeCommunProbPathway(C.cellchat)
C.cellchat <- aggregateNet(C.cellchat)



groupSize <- as.numeric(table(C.cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(C.cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Number of interactions")

netVisual_circle(C.cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Interaction weights/strength")

mat <- C.cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}



library(tidyverse)
df.net <- subsetCommunication(AT.cellchat,thresh = 0.05)
df<-df.net%>%filter(source=="TSK")
df<-df.net%>%filter(target=="Fibroblast")
pathways.show <- c("MIF") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(AT.cellchat, signaling =pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(AT.cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(AT.cellchat, signaling = pathways.show, layout = "chord")
unique(labels)
netVisual_bubble(AT.cellchat, 
                 sources.use = 15, 
                 targets.use = 6, 
                 remove.isolate = FALSE)

netVisual_bubble(AT.cellchat, 
                 sources.use = 2, 
                 targets.use = c(1:7), 
                 remove.isolate = FALSE)

netVisual_chord_gene(AT.cellchat, 
                     sources.use = c(1:7), 
                     targets.use =2, 
                     lab.cex =0.3,
                     legend.pos.y = 20)
##Pb组通讯
object_Pb = alldata$Pb
Pb.str <- GetAssayData(object = object_Pb,
                      assay = "RNA", slot = "data")

Pb.meta.data <- alldata[["Pb"]]@meta.data
Pb.cellchat <- createCellChat(object = Pb.str,
                             meta = Pb.meta.data,
                             group.by = "celltype")
levels(Pb.cellchat@idents)
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB,
                           search = "Secreted Signaling")
Pb.cellchat@DB <- CellChatDB.use
future::plan("multiprocess", workers = 6)
Pb.cellchat <- subsetData(Pb.cellchat)
Pb.cellchat <- identifyOverExpressedGenes(Pb.cellchat)
Pb.cellchat <- identifyOverExpressedInteractions(Pb.cellchat)
Pb.cellchat <- projectData(Pb.cellchat, PPI.mouse)

Pb.cellchat <- computeCommunProb(Pb.cellchat, raw.use = TRUE)
Pb.cellchat <- filterCommunication(Pb.cellchat, min.cells = 10)
Pb.cellchat <- computeCommunProbPathway(Pb.cellchat)
Pb.cellchat <- aggregateNet(Pb.cellchat)
groupSize <- as.numeric(table(Pb.cellchat@idents))

par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(Pb.cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "CNumber of interactions")

netVisual_circle(C.cellchat@net$count,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "PNumber of interactions")

netVisual_circle(Pb.cellchat@net$weight,
                 vertex.weight = groupSize,
                 weight.scale = T,
                 label.edge= F,
                 title.name = "Interaction weights/strength")

mat <- Pb.cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#合并比较
object.list <- list(Control = C.cellchat, Pb = Pb.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2],
                   edge.width.max = 12, title.name = paste0("Number of interactions - ",names(object.list)[i]))}

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison",stacked = T, do.stat = T)
gg2 <- rankNet(cellchat, mode = "comparison",stacked = F, do.stat = T)

netVisual_bubble(cellchat,sources.use = 1, targets.use = c(5:11), comparison = c(1,2), angle.x = 45)

library(tidyverse)
df.net <- subsetCommunication(C.cellchat,thresh = 0.05)
df<-df.net%>%filter(source=="TSK")
df<-df.net%>%filter(target=="Microglia")
df1<-df.net%>%filter(target=="Astrocyte")
pathways.show <- c("MK") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(C.cellchat, signaling =pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(AT.cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(AT.cellchat, signaling = pathways.show, layout = "chord")
unique(labels)
netVisual_bubble(Pb.cellchat, 
                 sources.use = 15, 
                 targets.use = 6, 
                 remove.isolate = FALSE)

netVisual_bubble(Pb.cellchat, 
                 sources.use = 4, 
                 targets.use = c(1:7), 
                 remove.isolate = FALSE)

netVisual_chord_gene(Pb.cellchat, 
                     sources.use = c(1:5), 
                     targets.use =2, 
                     lab.cex =0.3,
                     legend.pos.y = 20)

####小胶质细胞亚群分析####
###包的安装及数据读入
library(ggplot2) 
library(cowplot) 
library(paletteer)  
library(gplots)
library(ggpubr)    
library(ggsci) 
library(stringr)
library(rio)
source("bs_figure//custom_seurat_functions.R")
load("bs_code/bs_totalcell.Rdata")
load("bs_code/MG/bs_microglia.Rdata")
load("bs_code/MG/fe_microglia.Rdata")
load("bs_code/MG/fe-MG.Rdata")

###选取小胶质细胞亚群
View(sce.all@meta.data)
cell.used <- rownames(sce.all@meta.data[which(sce.all@meta.data$seurat_clusters%in% c(0,3,7)),])
length(cell.used)
###降为分析
sub_microglia <- subset(sce.all, cells = cell.used)
sub_microglia <- NormalizeData(sub_microglia, normalization.method = "LogNormalize", scale.factor = 1e4) 
sub_microglia <- FindVariableFeatures(sub_microglia, selection.method = 'vst', nfeatures = 2000)
sub_microglia <- ScaleData(sub_microglia, vars.to.regress = "mito.percent")
sub_microglia <- RunPCA(sub_microglia, features = VariableFeatures(object = sub_microglia))
sub_microglia <- FindNeighbors(sub_microglia, dims = 1:10)
sub_microglia <- FindClusters(sub_microglia, resolution = 0.5)
sub_microglia <- FindClusters(sub_microglia, resolution = 0.7)
sub_microglia <- FindClusters(sub_microglia, resolution = 0.6)
head(Idents(sub_microglia), 5)
table(sub_microglia$seurat_clusters)
sub_microglia <- RunUMAP(sub_microglia, dims = 1:10)
sub_microglia <- RunTSNE(sub_microglia, dims = 1:10)
DimPlot(sub_microglia, reduction = 'tsne', label = TRUE, pt.size = 0.6)
DimPlot(sub_microglia, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
DimPlot(sub_microglia, reduction = 'umap', label = TRUE, split.by = "sample", pt.size = 0.8)
##harmony去批次效应
sub_microglia <- RunHarmony(sub_microglia, "orig.ident")
names(sub_microglia@reductions)
sub_microglia <- RunUMAP(sub_microglia,dims = 1:10, 
                   reduction = "harmony")
sub_microglia <- RunTSNE(sub_microglia, dims = 1:10,reduction = "harmony")
##计算比例
table(sub_microglia$sample)
prop.table(table(Idents(sub_microglia)))
table(Idents(sub_microglia), sub_microglia$sample)
Cellratio <- prop.table(table(Idents(sub_microglia), sub_microglia$sample), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
##find marker
#总体差异基因
sub_microglia.markers0vs1 <- FindMarkers(sub_microglia, group.by="mg_celltype",
                                         ident.1 = "Lps_MG", 
                                         ident.2 = "Hom_MG", logfc.threshold = 0.15, min.pct = 0.25)
genes_to_check = c('Cx3cr1', 'C1qb', 'Nfkbia', 'Ifitm3',"Mt1","Tyrobp","Trem2","Atf3")
genes_to_check = c('Gpx4','Slc40a1', 'Fth1',"Hmox1","Tfrc","Nfe2l2", "Hspb1")
DotPlot(sub_microglia, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()
VlnPlot(subset(sub_microglia1, downsample = 1500),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
VlnPlot(subset(sub_microglia_6, downsample = 50),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)

#亚群合并及拆分
sub_microglia1 <- sub_microglia
cell.used <- rownames(sub_microglia1@meta.data[which(sub_microglia1@meta.data$seurat_clusters%in% c(0,1,2,3,4,5,6,7,8)),])
length(cell.used)
sub_microglia1 <- subset(sub_microglia1, cells = cell.used)
celltype=data.frame(ClusterID=0:8 ,
                    mgcelltype= 0:8) 
levels(sub_microglia1)
#名字和群匹配  
celltype[celltype$ClusterID %in% c(2,0),2]='0' 
celltype[celltype$ClusterID %in% c(1),2]='1'  
celltype[celltype$ClusterID %in% c(3),2]='2'  
celltype[celltype$ClusterID %in% c(4),2]='3'  
celltype[celltype$ClusterID %in% c(5),2]='4'   
celltype[celltype$ClusterID %in% c(6),2]='5'  
celltype[celltype$ClusterID %in% c(7),2]='6'  
celltype[celltype$ClusterID %in% c(8),2]='7'  
head(celltype)
celltype
table(celltype$mgcelltype)
sub_microglia1@meta.data$mgcelltype = "NA"

for(i in 1:nrow(celltype)){
  sub_microglia1@meta.data[which(sub_microglia1@meta.data$seurat_clusters == celltype$ClusterID[i]),'mgcelltype'] <- celltype$mgcelltype[i]}
table(sub_microglia1@meta.data$mgcelltype)
p <- DimPlot(sub_microglia1, reduction = "tsne",group.by="mgcelltype",
             split.by = "sample",label = TRUE, pt.size = 0.3)
p
DimPlot(sub_microglia1, reduction = 'tsne', group.by="mgcelltype",
        label = TRUE, split.by = "sample", pt.size = 0.6)
DimPlot(sub_microglia1, reduction = 'tsne', group.by="mgcelltype",
        label = TRUE,  pt.size = 0.6)
DimPlot(sub_microglia1, reduction = 'umap', group.by="mgcelltype",
        label = TRUE, split.by = "sample", pt.size = 0.6)
cell.used <- rownames(sub_microglia1@meta.data[which(sub_microglia1@meta.data$mgcelltype%in% c(0,1,2,3,4,5,6)),])
length(cell.used)
cell.used <- rownames(sub_microglia1@meta.data[which(sub_microglia1@meta.data$mgcelltype%in% c(0)),])
length(cell.used)
sub1_microglia1 <- subset(sub_microglia1, cells = cell.used)
sub1_microglia1$mgcelltype = ifelse(sub1_microglia1@assays$RNA@counts['Gpx4',]>1,'pos','neg')
sub_microglia1$mgcelltype[match(colnames(sub1_microglia1),colnames(sub_microglia1))] =  sub1_microglia1$mgcelltype
table(sub_microglia1$mgcelltype)
##计算比例
sub_microglia1@active.ident <- plyr::mapvalues(x = sub_microglia1@active.ident, 
                                        from = celltype$ClusterID, to = celltype$mgcelltype)
sub_microglia1@active.ident<-factor(sub_microglia1@active.ident,levels = unique(celltype$mgcelltype))
table(sub_microglia1$sample)
prop.table(table(Idents(sub_microglia1)))
table(Idents(sub_microglia1), sub_microglia1$sample)
Cellratio <- prop.table(table(Idents(sub_microglia1), sub_microglia1$sample), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
###降为分析
sub_microglia1 <- subset(sub_microglia1, cells = cell.used)
save(sub_microglia1,file = "bs_code/MG/fe_microglia.Rdata")
#铁死亡打分
ferroptosis_all <- list(c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                          "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                          "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                          "Cdo1","Me1","Coq2","Gch1","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                          'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                          "Phgdh","Hmgb1","Pparg"))
WNT_features <- ferroptosis_all
sub_microglia1 <- AddModuleScore(sub_microglia1,features = WNT_features,ctrl = 100,name = "WNT_features")
p <- ggboxplot(sub_microglia1@meta.data, x="mgcelltype", y="WNT_features1", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#1E90FF", "#FF6347"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")

compare_means(WNT_features1 ~ sample, data = sub_microglia1@meta.data,group.by = "mgcelltype")
p + stat_compare_means(aes(sample = sample),label = "p.signif")

#换颜色
p <- DimPlot(sub_microglia, reduction = "umap",group.by="celltype",label = TRUE, pt.size = 0.3) + NoLegend()+labs(x = "UMAP1", y = "UMAP2",title = "Celltype")
plot4 = plot3 + scale_color_npg()

cell_type_cols <- c(brewer.pal(9, "Set2"), "#FF34B3", "#BC8F8F", "#20B2AA", "#00F5FF", 
                    "#FFA500", "#ADFF2F", "#FF6A6A", "#7FFFD4", "#AB82FF", "#90EE90", 
                    "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B",
                    "#FF3030", "#7CFC00", "#000000", "#708090")
DimPlot(sub_microglia1, label = T, cols = cell_type_cols,  pt.size = 0.7, repel = T,split.by = "sample",reduction = "tsne")
#亚群差异基因比较
sub_microglia.markers <- FindAllMarkers(sub_microglia,   # 只返回positive基因
                                  min.pct = 0.25)
sub_microglia.markers <- FindAllMarkers(sub_microglia1,
                                         only.pos = TRUE,
                                         logfc.threshold = 0.2,# 只返回positive基因
                                         min.pct = 0.2)
sub_microglia.markers <- FindMarkers(sub_microglia, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.2)

sub_microglia.markers0vs1 <- FindMarkers(sub_microglia, group.by="seurat_clusters",
                        ident.1 = "1", 
                        ident.2 = "0", logfc.threshold = 0.15, min.pct = 0.25)
sub_microglia.markers0vs2 <- FindMarkers(sub_microglia, group.by="seurat_clusters",
                                         ident.1 = "2", 
                                         ident.2 = "0", logfc.threshold = 0.15, min.pct = 0.25)
sub_microglia.markers0vs3 <- FindMarkers(sub_microglia, group.by="seurat_clusters",
                                         ident.1 = "3", 
                                         ident.2 = "0", logfc.threshold = 0.15, min.pct = 0.25)
sub_microglia.markers0vs4 <- FindMarkers(sub_microglia, group.by="seurat_clusters",
                                         ident.1 = "4", 
                                         ident.2 = "0", logfc.threshold = 0.15, min.pct = 0.25)
sub_microglia.markers0vs5 <- FindMarkers(sub_microglia, group.by="seurat_clusters",
                                         ident.1 = "5", 
                                         ident.2 = "0", logfc.threshold = 0.15, min.pct = 0.25)
sub_microglia.markers0vs6 <- FindMarkers(sub_microglia, group.by="seurat_clusters",
                                         ident.1 = "6", 
                                         ident.2 = "0", logfc.threshold = 0.15, min.pct = 0.25)
write.csv(sub_microglia.markers,file = "bs_table/yqmarker_MG/sub__microgliamarkers.csv")
write.csv(sub_microglia.markers1,file = "bs_table/femarker_MG/sub__microgliamarkers.csv")
write.csv(sub_microglia.markers0vs1,file = "bs_table/yqmarker_MG/sub_microglia.markers0vs1.csv")
write.csv(sub_microglia.markers0vs2,file = "bs_table/yqmarker_MG/sub_microglia.markers0vs2.csv")
write.csv(sub_microglia.markers0vs3,file = "bs_table/yqmarker_MG/sub_microglia.markers0vs3.csv")
write.csv(sub_microglia.markers0vs4,file = "bs_table/yqmarker_MG/sub_microglia.markers0vs4.csv")
write.csv(sub_microglia.markers0vs5,file = "bs_table/yqmarker_MG/sub_microglia.markers0vs5.csv")
write.csv(sub_microglia.markers0vs6,file = "bs_table/yqmarker_MG/sub_microglia.markers0vs6.csv")
meta.data <- sub_microglia@meta.data
table(meta.data$mg_celltype)
sub_microglia.markers0vs3 <- FindMarkers(sub_microglia, group.by="mg_celltype",
                                         ident.1 = "Mid_MG1", 
                                         ident.2 = "Hom_MG", logfc.threshold = 0.5, min.pct = 0.25)
markers0vs5 = row.names(sub_microglia.markers0vs5)
markers0vs5 <- bitr(markers0vs5,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
markers0vs5 = markers0vs5[,2]
markers0vs5_kegg <- enrichKEGG(
  markers0vs5,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE
)

markers0vs6 = row.names(sub_microglia.markers0vs6)
markers0vs6 <- bitr(markers0vs6,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
markers0vs6 = markers0vs6[,2]
markers0vs6_kegg <- enrichKEGG(
  markers0vs6,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE
)
write.csv(markers0vs5_kegg,file = "bs_table/MG/markers0vs5_kegg.csv")
write.csv(markers0vs6_kegg,file = "bs_table/MG/markers0vs6_kegg.csv")
write.csv(sub_microglia.markers0vs5,file = "bs_table/MG/sub_microglia.markers0vs5.csv")
write.csv(sub_microglia.markers0vs6,file = "bs_table/MG/sub_microglia.markers0vs6.csv.")
write.csv(sub_microglia.markers,file = "bs_table/MG/sub__microgliamarkers.")
write.csv(sub_microglia.markers,file = "bs_table/MG/sub__microgliamarkers1.csv")
save(sub_microglia1,sub_microglia.markers, file = "bs_code/MG/fe-MG.Rdata")
##富集分析
dif_micgene <- import("bs_table/MG/dif_gene1.csv")
dif_micgene[dif_micgene == ""]<-NA
total <- list(Hom_MG = dif_micgene$`Hom-MG`,M1_like = dif_micgene$`M1-like`,M2_like = dif_micgene$`M2-like`,
              Pha_MG = dif_micgene$`Pha-MG`, Mid_MG = dif_micgene$`Mid-MG`, Pro_DAM = dif_micgene$`Pro-DAM`,
              Pha_DAM = dif_micgene$`Pha-DAM`)
i=1
for (i in 1:length(total)) {
  total[[i]] <- bitr(total[[i]],
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)
  total[[i]] = total[[i]][,2]
}
lapply(total, head)

go_result <- list(Hom_MGgo=NA,M1_likego=NA,M2_likego=NA,Pha_MGgo=NA,Mid_MGgo=NA,Pro_DAMgo=NA,Pha_DAMgo=NA)
length(total)
i = 1
for (i in 1:length(total)){
  go_result[[i]]<- enrichGO(
    total[[i]],
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 3,
    maxGSSize = 20,
    readable = TRUE,
    pool = FALSE
  )
}

#批量导出go结果
filename<-c("Hom_MGgo","M1_likego","M2_likego","Pha_MGgo","Mid_MGgo","Pro_DAMgo","Pha_DAMgo")

for(i in 1:length(go_result)){
  write.csv(go_result[[i]],file = paste0("bs_table/MG/",filename[i],".csv"))
}
##画图
genes_to_check = c('Cx3cr1', 'Lpl', 'Nfkbia', 'Ifitm3',"Apoe","Tyrobp","Trem2","Maf")
DotPlot(sub_microglia, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()
VlnPlot(subset(sub_microglia, downsample = 500),group.by = "seurat_clusters",
        features = unique(genes_to_check),ncol = 2)
VlnPlot(subset(sub_microglia_6, downsample = 50),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(6)),])
length(cell.used)
sub_microglia_6 <- subset(sub_microglia, cells = cell.used)
##亚群命名
mg_celltype=data.frame(ClusterID=0:6 ,
                       mg_celltype= 0:6) 
mg_celltype[mg_celltype$ClusterID %in% c(0),2]='Hom-MG' 
mg_celltype[mg_celltype$ClusterID %in% c(1),2]='M1-like'  
mg_celltype[mg_celltype$ClusterID %in% c(2),2]='M2-like'  
mg_celltype[mg_celltype$ClusterID %in% c(3),2]='M2-like'  
mg_celltype[mg_celltype$ClusterID %in% c(4),2]='Mid-MG'   
mg_celltype[mg_celltype$ClusterID %in% c(5),2]='DAM1'  
mg_celltype[mg_celltype$ClusterID %in% c(6),2]='DAM2'
mg_celltype
table(mg_celltype$mg_celltype)
sub_microglia@meta.data$mg_celltype = "NA"
for(i in 1:nrow(mg_celltype)){
  sub_microglia@meta.data[which(sub_microglia@meta.data$RNA_snn_res.0.3 == mg_celltype$ClusterID[i]),
                          'mg_celltype'] <- mg_celltype$mg_celltype[i]}
table(sub_microglia@meta.data$mg_celltype)
DimPlot(sub_microglia, reduction = "umap",group.by="mg_celltype",label = TRUE, pt.size = 0.5)
save(sub_microglia,sub_microglia.markers,mg_celltype, file = "bs_microglia.Rdata")
sub_microglia@active.ident <- plyr::mapvalues(x = sub_microglia@active.ident, 
                                        from = mg_celltype$ClusterID, to = mg_celltype$mg_celltype)
sub_microglia@active.ident<-factor(sub_microglia@active.ident,levels = unique(mg_celltype$mg_celltype))
##计算亚群比例
p1 = plot.clusters.group(data = sub_microglia,clusters =  "mg_celltype", 
                         xlab = "Celltype", log = TRUE,
                         group = "orig.ident",legend.title = "Sample",
                         widths = c(3,1),color = 1)
p1 = plot.clusters.group(data = sub_microglia,clusters =  "mg_celltype", 
                         xlab = "Celltype", log = F,
                         group = "sample",legend.title = "Sample",
                         widths = c(3,1),color = 1)
p1
table(sub_microglia$sample)
prop.table(table(Idents(sub_microglia)))
table(Idents(sub_microglia), sub_microglia$sample)
Cellratio <- prop.table(table(Idents(sub_microglia), sub_microglia$sample), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("Celltype","Group","Percentage")
colourCount = length(unique(Cellratio$Celltype))
ggplot(Cellratio) + 
  geom_bar(aes(x =Group, y= Percentage, fill = Celltype),stat = "identity", width = 0.5,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Percentage')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave('bs_figure/mg.pdf',height = 5,width = 11)
save(sub_microglia,sub_microglia.markers, file = "bs_code/MG.Rdata")
##拟时序分析
rm(list=ls())
load("bs_code/MG.Rdata")
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(devtools)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(mclust, quietly = TRUE)
devtools::install_github('cole-trapnell-lab/monocle3')
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github("cole-trapnell-lab/monocle3")
#抽取部分细胞进行拟时序分析
allCells=names(Idents(sub_microglia))
allType = levels(Idents(sub_microglia))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sub_microglia)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,3000,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sub_microglia[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(cg_sce))
DimPlot(cg_sce, label = F, pt.size = 0.5)

#拟时序分析步骤slingshot包
cluster.ident<-FetchData(cg_sce,c('mg_celltype','PC_1','PC_2','UMAP_1','UMAP_2'))
sce <- as.SingleCellExperiment(cg_sce)

sce <- slingshot(sce, clusterLabels = 'seurat',reducedDim = "UMAP", start.clus = c("Hom_MG"))
cell.used <- cluster.ident
sce <- SingleCellExperiment(assays=list(counts=cg_sce@assays$RNA@counts[,],
                                        logcounts=cg_sce@assays$RNA@data[,]),
                            reducedDims=SimpleList(PCA=cluster.ident[,c('PC_1','PC_2')], 
                                                   UMAP=cluster.ident[,c('UMAP_1','UMAP_2')]))
colData(sce)$seurat <- cluster.ident$mg_celltype
sce <- slingshot(sce, start.clus = c("Hom-MG"), 
                 end.clus =c("M2-like","M1-like" ,"DAM1","DAM2"), 
                 clusterLabels = sce$seurat, reducedDim = 'UMAP') 
summary(sce$slingPseudotime_4)

colnames(colData(sce))
colData(sce)
coul = brewer.pal(nlevels(as.factor(cluster.ident$mg_celltype)), "Set2")
colors.type=coul[as.numeric(as.factor(cluster.ident$mg_celltype))]
plot(reducedDims(sce)$UMAP, col = colors.type, pch=16, asp = 0.5,cex = 0.5)
lines(SlingshotDataSet(sce), lwd=2, col='black')

library(gam)
sling.pseu<-data.frame(sce$seurat,sce$slingPseudotime_1,sce$slingPseudotime_2,
                       sce$slingPseudotime_3,sce$slingPseudotime_4)
rownames(sling.pseu)<-colnames(sce)
head(sling.pseu)
t <- sling.pseu$sce.slingPseudotime_4
sling.cells <- rownames(sling.pseu)[!is.na(t)] #去na
head(sling.cells)
length(sling.cells)
table(sling.pseu[sling.cells,]$sce.seurat)
test<-FetchData(cg_sce,c("mg_celltype","Ctsb"))
dim(test)
test<-test[sling.cells,]
test$pseu<-sling.pseu[sling.cells,]$sce.slingPseudotime_4
ggplot(test,aes(pseu,Ctsb))+ geom_point(aes(colour=factor(mg_celltype)),alpha = 0.8)+ 
  geom_smooth(colour="gray30",size=1.2) +
  theme_bw()
#拟时序热图

#差异基因图
head(sub_microglia.markers)
diff_cell <- sub_microglia.markers
#显著性
diff_cell$label <- ifelse(diff_cell$p_val_adj<0.01,"adjust P-val<0.01", "adjust P-val>=0.01")
top10_cell = diff_cell %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
diff_cell$size <- case_when(!(diff_cell$gene %in% top10_cell$gene)~1,
                            diff_cell$gene %in% top10_cell$gene ~2)
dt <- filter(diff_cell,size ==1)

dfbar <- data.frame(x=c("0","1","2","3","4","5","6"),
                    y=c(3.5,3,3,3,3,3,3))
p <- ggplot() + 
  geom_col(data = dfbar,
           mapping = aes(x = x, y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width = 0.4)+
  geom_jitter(data = top10_cell,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width = 0.4)

#添加显著基因
dfcol <- data.frame(x=c("0","1","2","3","4","5","6"),
                  y=0,
                  label=c("0","1","2","3","4","5","6"))
mycol <- c("#E64B357F","#00A0877F","#3C54887F","#F39B7F7F","#E64B357F","#00A0877F","#3C54887F")
library(ggrepel) 
p1 <-  p + geom_tile(data = dfcol,
                                aes(x=x,y=y),
                                height=0.5,
                                color = "black",
                                fill = mycol,
                                alpha = 0.6,
                                show.legend = F) +
  geom_text_repel(
    data=top10_cell,
    aes(x = cluster, y = avg_log2FC, label = gene),
    size = 2,
    arrow = arrow(length = unit(0.008,"npc"),
                  type = "open",ends = "last"))
#对图进行修饰
p2 <- p1+
  scale_color_manual(name = NULL,
                     values = c("red","grey"))+
  labs(x="cluster",y = "avg_log2FC")+
  geom_text(data = dfcol,
            aes(x = x, y = y, label = label),
            size = 4,
            color = "black")+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 12,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  )


####小胶质与疾病相关的分析####
library(rio)
library(Hmisc)
mg_marker <- import("bs_table/MG/sub__microgliamarkers.csv")
select_mg_marker <- mg_marker[c(mg_marker$cluster == "6")|c(mg_marker$cluster == "5"),]
select_mg_marker <- select_mg_marker[,c("gene","avg_log2FC")]
DIS_GENE <- import("bs_table/zys-lysome.xlsx")
Lysome <- tolower(DIS_GENE$lysome)
Lysome <- capitalize(Lysome)
DEGs <- DIS_GENE$DEGs
Inter_lysome_DEGs <- intersect(DEGs, Lysome)
Inter_lysome_DEGs <- as.data.frame(Inter_lysome_DEGs)
write.csv(Inter_lysome_DEGs,file = "bs_table/zys-DEGslysome.csv")

DIS_GENE <- DIS_GENE[,c("AD","Aging","ALS","PD","LPS")]
DIS_GENE <- na.omit(DIS_GENE)
AD <- capitalize(DIS_GENE$AD)
Aging <- capitalize(DIS_GENE$Aging)
PD <- capitalize(DIS_GENE$PD)
ALS <- capitalize(DIS_GENE$ALS)
LPS <- capitalize(DIS_GENE$LPS)
AD_MG <- select_mg_marker[c(select_mg_marker$gene) %in% AD,]
Aging_MG <- select_mg_marker[c(select_mg_marker$gene) %in% Aging,]
PD_MG <- select_mg_marker[c(select_mg_marker$gene) %in% PD,]
ALS_MG <- select_mg_marker[c(select_mg_marker$gene) %in% ALS,]
LPS_MG <- select_mg_marker[c(select_mg_marker$gene) %in% LPS,]
total_dismg <- rbind(AD_MG,Aging_MG,PD_MG,ALS_MG,LPS_MG)
write.csv(total_dismg,file = "bs_table/MG/total_dismg.csv")
###AD差异基因富集GSEA富集分析
AD_dat <- import("bs_table/MG/total_dismg1.csv")
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(msigdbr)
gene <- bitr(AD_dat$Gene,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Mm.eg.db)

AD_dat <- merge(gene,AD_GSEA,by = 1,all = F)
geneList = AD_dat[,1]
names(geneList) = as.character(AD_dat[,2])
head(geneList)
geneList <- geneList[!duplicated(names(geneList))]
geneList = sort(geneList, decreasing = TRUE)
AD_gsego <- gseGO(geneList = geneList,#排序后的基因列表，一般根据logFC进行排序#有100个基因的数量要求
               OrgDb = org.Mm.eg.db,
               ont = "BP",#可选择bp.MF,CC,ALL
               nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
               minGSSize = 1,#最小基因集的基因数
               maxGSSize = 500,#最大基因集的基因数
               pvalueCutoff = 10,
               pAdjustMethod = "BH",#p值的阈???
               verbose      = FALSE)#是否输出提示信息，默认为false

geneList = gene[,2]
AD_go <- enrichGO(
  geneList,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 2,
  maxGSSize = 20,
  readable = TRUE,
  pool = FALSE
)
write.csv(AD_go,file = "bs_table/MG/dismg_go.csv")

AD_go_CC <- enrichGO(
  geneList,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1,
  minGSSize = 2,
  maxGSSize = 20,
  readable = TRUE,
  pool = FALSE
)

write.csv(AD_go_CC,file = "bs_table/MG/dismg_go_CC.csv")

dis_kegg <- enrichKEGG(
  geneList,
  organism = "mmu",
  keyType = "kegg",
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 2,
  maxGSSize = 500,
  qvalueCutoff = 1,
  use_internal_data = FALSE
)

write.csv(dis_kegg,file = "bs_table/MG/dismg_kegg.csv")

p = dotplot(AD_go, showCategory=10)
p + scale_y_discrete(labels=function(x) str_wrap(x,width = 50))
dotplot(dis_kegg, showCategory=10)
dotplot(gsegmt, showCategory=30)
cnetplot(AD_go, colorEdge = TRUE, circular = TRUE) + 
  ggraph::scale_edge_color_discrete(labels = function(x) str_wrap(x, 20) )


####TREM2调控小胶质细胞亚群分析####
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
sub_microglia.markers0vs6 <- import("bs_table/MG/sub_microglia.markers0vs6.csv")
sub_microglia.markers0vs6 <- sub_microglia.markers0vs6[,c(1,3)]
colnames(sub_microglia.markers0vs6) <- c("Gene","LogFC")
#ID转换
all_marker = sub_microglia.markers0vs6[,1]
all_marker <- bitr(all_marker,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
all_marker <- merge(all_marker,sub_microglia.markers0vs6,by = 1,all = F)
# 1: 提取logFC值，并储存在一个向量中
geneList = all_marker[,3]
# 2: 对geneList进行命名
names(geneList) = as.character(all_marker[,2])
head(geneList)
# 3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)
sub_microglia.markers0vs6_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
      OrgDb = org.Mm.eg.db,
      ont = "BP",#可选择bp.MF,CC,ALL
      nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
      minGSSize = 1,#最小基因集的基因数
      maxGSSize = 500,#最大基因集的基因数
      pvalueCutoff = 10,
      pAdjustMethod = "BH",#p值的阈???
      verbose      = FALSE)#是否输出提示信息，默认为false
sub_microglia.markers0vs6_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                               organism     = 'mmu',
                                               nPerm        = 1000,
                                               minGSSize    = 3,
                                               pvalueCutoff = 1,
                                               verbose      = FALSE)
write.csv(sub_microglia.markers0vs6_GSEAGO,file = "bs_table/MG/sub_microglia.markers0vs6_GSEAGO.csv")
write.csv(sub_microglia.markers0vs6_GSEAKEGG,file = "bs_table/MG/sub_microglia.markers0vs6_GSEAKEGG.csv")
gseaplot2(sub_microglia.markers0vs6_GSEAKEGG, geneSetID = 16, title = sub_microglia.markers0vs6_GSEAKEGG$Description[16],
          subplots = 1:3,
          color = "#0AFF99")
gseaplot2(sub_microglia.markers0vs6_GSEAKEGG, geneSetID = 71, title = sub_microglia.markers0vs6_GSEAKEGG$Description[71],
          subplots = 1:3,
          color = "#0AFF99")
gseaplot2(sub_microglia.markers0vs6_GSEAKEGG, geneSetID = 63, title = sub_microglia.markers0vs6_GSEAKEGG$Description[63],
          subplots = 1:3,
          color = "#0AFF99")
gseaplot2(sub_microglia.markers0vs6_GSEAKEGG, geneSetID = 195, title = sub_microglia.markers0vs6_GSEAKEGG$Description[195],
          subplots = 1:3,
          color = "#0AFF99")
gseaplot2(sub_microglia.markers0vs6_GSEAKEGG, geneSetID = 227, title = sub_microglia.markers0vs6_GSEAKEGG$Description[227],
          subplots = 1:3,
          color = "#0AFF99")
#选取
genes_to_check = c("Trem2",	"Tyrobp",	"Syk","Nfam1", 'Cd33', 'Lgals3', "Nlrp3","Il1b","Il1rn","Tlr4",
                   "Tpm1","Cd86","Tnf","Apoe","Itgb2",	"Manf",	"Irf3","Akt1", "Lpl",
                   "Fxyd5",	"Cd53",	"Il10rb",	"Lilrb4a","Axl","Stat6","Icam2","Adam10","Adam17",
                   "Mapk1",	"Nfe2l2",	"Grb2",	"Adam15",	"Cst7",	"Cd52",	"Ctsd",	"Itgax",
                   "Cd68","Ccl6",	"Ccl2",	"Lyz2",	"C1qc",	"Itgb2","Ctsb","Irak3","Tlr2","Myd88",
                   "Plcg1","Ras","Itpr3","Rasgrp1","Itpr1","Pten")
#6号亚群Trem2值
group6 <- AverageExpression(sub_microglia,features = unique(genes_to_check),group.by = "sample")
group6 <- as.data.frame(group6)
group6$FC <- group6$RNA.Pb/group6$RNA.Control
group6$logFC <- log2(group6$FC)
FC_t1=1.2
FC_t2=0.85
k1 = (group6$FC < FC_t2)
k2 = (group6$FC > FC_t1)
group6 <- mutate(group6, change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(group6$change)
write.csv(group6,file = "bs_table/different/group6.csv")
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(6)),])
length(cell.used)
sub_microglia_6 <- subset(sub_microglia, cells = cell.used)
AverageExpression(sub_microglia_6,features = unique(genes_to_check),group.by = "sample")
VlnPlot(subset(sub_microglia_6, downsample = 50),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
sub_microglia_groupmarker <- FindMarkers(sub_microglia_6, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.2)
#5号亚群Trem2
group5 <- AverageExpression(sub_microglia_5,features = unique(genes_to_check),group.by = "sample")
group5 <- as.data.frame(group5)
group5$FC <- group5$RNA.Pb/group5$RNA.Control
group5$logFC <- log2(group5$FC)
FC_t1=1.2
FC_t2=0.85
k1 = (group5$FC < FC_t2)
k2 = (group5$FC > FC_t1)
group5 <- mutate(group5, change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(group5$change)
write.csv(group5,file = "bs_table/different/group5.csv")
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(5)),])
length(cell.used)
sub_microglia_5 <- subset(sub_microglia, cells = cell.used)
AverageExpression(sub_microglia_5,features = unique(genes_to_check),group.by = "sample")
VlnPlot(subset(sub_microglia_5),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
#4号亚群Trem2
group4 <- AverageExpression(sub_microglia_4,features = unique(genes_to_check),group.by = "sample")
group4 <- as.data.frame(group4)
group4$FC <- group4$RNA.Pb/group4$RNA.Control
group4$logFC <- log2(group4$FC)
FC_t1=1.2
FC_t2=0.85
k1 = (group4$FC < FC_t2)
k2 = (group4$FC > FC_t1)
group4 <- mutate(group4, change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(group4$change)
write.csv(group4,file = "bs_table/different/group4.csv")
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(4)),])
length(cell.used)
sub_microglia_4 <- subset(sub_microglia, cells = cell.used)
AverageExpression(sub_microglia_4,features = unique(genes_to_check),group.by = "sample")
VlnPlot(subset(sub_microglia_4),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
#3号亚群Trem2
group3 <- AverageExpression(sub_microglia_3,features = unique(genes_to_check),group.by = "sample")
group3 <- as.data.frame(group3)
group3$FC <- group3$RNA.Pb/group3$RNA.Control
group3$logFC <- log2(group3$FC)
FC_t1=1.2
FC_t2=0.85
k1 = (group3$FC < FC_t2)
k2 = (group3$FC > FC_t1)
group3 <- mutate(group3, change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(group3$change)
write.csv(group3,file = "bs_table/different/group3.csv")
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(3)),])
length(cell.used)
sub_microglia_3 <- subset(sub_microglia, cells = cell.used)
AverageExpression(sub_microglia_3,features = unique(genes_to_check),group.by = "sample")
VlnPlot(subset(sub_microglia_3, downsample = 500),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
#2号亚群Trem2
group2 <- AverageExpression(sub_microglia_2,features = unique(genes_to_check),group.by = "sample")
group2 <- as.data.frame(group2)
group2$FC <- group2$RNA.Pb/group2$RNA.Control
group2$logFC <- log2(group2$FC)
FC_t1=1.2
FC_t2=0.85
k1 = (group2$FC < FC_t2)
k2 = (group2$FC > FC_t1)
group2 <- mutate(group2, change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(group2$change)
write.csv(group2,file = "bs_table/different/group2.csv")
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(2)),])
length(cell.used)
sub_microglia_2 <- subset(sub_microglia, cells = cell.used)
AverageExpression(sub_microglia_2,features = unique(genes_to_check),group.by = "sample")
VlnPlot(subset(sub_microglia_2),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
#1号亚群Trem2
group1 <- AverageExpression(sub_microglia_1,features = unique(genes_to_check),group.by = "sample")
group1 <- as.data.frame(group1)
group1$FC <- group1$RNA.Pb/group1$RNA.Control
group1$logFC <- log2(group1$FC)
FC_t1=1.2
FC_t2=0.85
k1 = (group1$FC < FC_t2)
k2 = (group1$FC > FC_t1)
group1 <- mutate(group1, change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(group1$change)
write.csv(group1,file = "bs_table/different/group1.csv")
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(1)),])
length(cell.used)
sub_microglia_1 <- subset(sub_microglia, cells = cell.used)
AverageExpression(sub_microglia_1,features = unique(genes_to_check),group.by = "sample")
VlnPlot(subset(sub_microglia_1),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
#0号亚群Trem2
group0 <- AverageExpression(sub_microglia_0,features = unique(genes_to_check),group.by = "sample")
group0 <- as.data.frame(group0)
group0$FC <- group0$RNA.Pb/group0$RNA.Control
group0$logFC <- log2(group0$FC)
FC_t1=1.2
FC_t2=0.85
k1 = (group0$FC < FC_t2)
k2 = (group0$FC > FC_t1)
group0 <- mutate(group0, change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(group0$change)
write.csv(group0,file = "bs_table/different/group0.csv")
cell.used <- rownames(sub_microglia@meta.data[which(sub_microglia@meta.data$seurat_clusters%in% c(0)),])
length(cell.used)
sub_microglia_0 <- subset(sub_microglia, cells = cell.used)
AverageExpression(sub_microglia_1,features = unique(genes_to_check),group.by = "sample")
VlnPlot(subset(sub_microglia_0),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
#保存所有亚群
save(sub_microglia,sub_microglia.markers, file = "bs_code/MG.Rdata")

####内皮细胞铁死亡分析####
cell.used <- rownames(sce.all@meta.data[which(sce.all@meta.data$seurat_clusters%in% c(12)),])
length(cell.used)
end.all <- subset(sce.all, cells = cell.used)
DimPlot(end.all, reduction = 'umap', label = TRUE,  pt.size = 0.8)
diff.end.gene <- FindMarkers(end.all, group.by="sample",
                        ident.1 = "Pb", ident.2 = "Control", 
                        logfc.threshold = 0.1, min.pct = 0.1)
ferroptosis_all <- c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                     "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                     "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                     "Cdo1","Me1","Coq2","Gch1","Slc25a39","Sqstm1","Hif1a",
                     'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                     "Phgdh","Hmgb1","Pparg")
genes_to_check = c("Gch1",	"Pebp1",	"Slc40a1","Tfrc")
VlnPlot(subset(end.all, downsample = 50),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
write.csv(diff.end.gene,file = "bs_table/end/diff.end.gene.csv")
##1号亚群的富集分析
library(rio)
Fe_marker1 <- import("bs_table/end/diff.end.gene.csv")
Fe_marker1 <- Fe_marker1[,c(1:3)]
colnames(Fe_marker1) <- c("Gene","p.value","LogFC")
ferroptosis_1 <- intersect(Fe_marker1$Gene,ferroptosis_all)
ferroptosis_1 <- as.data.frame(ferroptosis_1)
names(ferroptosis_1) <- "Gene"
ferroptosis_1_dat <- merge(ferroptosis_1,Fe_marker1,by = 1,all = F)
VlnPlot(subset(end.all, downsample = 50),group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
write.csv(ferroptosis_1_dat,file = "bs_table/end/ferroptosis_1_dat.csv")
####所有小胶质细胞亚群不同组别差异基因及功能富集####
##找差异基因
mig1_marker<-FindMarkers(sub_microglia_1, group.by="sample",
                        ident.1 = "Pb", 
                        ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.1)
mig2_marker<-FindMarkers(sub_microglia_2, group.by="sample",
                         ident.1 = "Pb", 
                         ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.1)
mig3_marker<-FindMarkers(sub_microglia_3, group.by="sample",
                         ident.1 = "Pb", 
                         ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.1)
mig4_marker<-FindMarkers(sub_microglia_4, group.by="sample",
                         ident.1 = "Pb", 
                         ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.1)
mig5_marker<-FindMarkers(sub_microglia_5, group.by="sample",
                         ident.1 = "Pb", 
                         ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.1)
mig6_marker<-FindMarkers(sub_microglia_6, group.by="sample",
                         ident.1 = "Pb", 
                         ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.1)
mig0_marker<-FindMarkers(sub_microglia_0, group.by="sample",
                         ident.1 = "Pb", 
                         ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.1)
write.csv(mig1_marker,file = "bs_table/mig1_marker.csv")
write.csv(mig2_marker,file = "bs_table/mig2_marker.csv")
write.csv(mig3_marker,file = "bs_table/mig3_marker.csv")
write.csv(mig4_marker,file = "bs_table/mig4_marker.csv")
write.csv(mig5_marker,file = "bs_table/mig5_marker.csv")
write.csv(mig6_marker,file = "bs_table/mig6_marker.csv")
write.csv(mig0_marker,file = "bs_table/mig0_marker.csv")
#保存marker
save(mig1_marker,mig2_marker,mig3_marker,
     mig4_marker,mig5_marker,mig6_marker,mig0_marker,file = "bs_code/mg_marker_sample.Rdata")
all_marker = sub_microglia.markers0vs6[,1]
all_marker <- bitr(all_marker,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
all_marker <- merge(all_marker,sub_microglia.markers0vs6,by = 1,all = F)
##GSEA分析
#0号群
mig0_marker <- import("bs_table/wfx/mig0_marker.csv")
mig0_marker <- mig0_marker[,c(1,3)]
colnames(mig0_marker) <- c("Gene","LogFC")
mig_marker0 = mig0_marker[,1]
mig_marker0 <- bitr(mig_marker0,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
mig0_marker <- merge(mig_marker0,mig0_marker,by = 1,all = F)
geneList = mig0_marker[,3]
names(geneList) = as.character(mig0_marker[,2])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
sub_microglia.markers0group_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                          OrgDb = org.Mm.eg.db,
                                          ont = "BP",#可选择bp.MF,CC,ALL
                                          nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                          minGSSize = 1,#最小基因集的基因数
                                          maxGSSize = 500,#最大基因集的基因数
                                          pvalueCutoff = 10,
                                          pAdjustMethod = "BH",#p值的阈???
                                          verbose      = FALSE)#是否输出提示信息，默认为false
sub_microglia.markers0group_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                               organism     = 'mmu',
                                               nPerm        = 1000,
                                               minGSSize    = 3,
                                               pvalueCutoff = 1,
                                               verbose      = FALSE)
write.csv(sub_microglia.markers0group_GSEAGO,file = "bs_table/wfx/sub_microglia.markers0group_GSEAGO.csv")
write.csv(sub_microglia.markers0group_GSEAKEGG,file = "bs_table/wfx/sub_microglia.markers0group_KEGG.csv")
#1号群
mig1_marker <- import("bs_table/wfx/mig1_marker.csv")
mig1_marker <- mig1_marker[,c(1,3)]
colnames(mig1_marker) <- c("Gene","LogFC")
mig_marker1 = mig1_marker[,1]
mig_marker1 <- bitr(mig_marker1,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
mig1_marker <- merge(mig_marker1,mig1_marker,by = 1,all = F)
geneList = mig1_marker[,3]
names(geneList) = as.character(mig1_marker[,2])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
sub_microglia.markers1group_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                            OrgDb = org.Mm.eg.db,
                                            ont = "BP",#可选择bp.MF,CC,ALL
                                            nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                            minGSSize = 1,#最小基因集的基因数
                                            maxGSSize = 500,#最大基因集的基因数
                                            pvalueCutoff = 10,
                                            pAdjustMethod = "BH",#p值的阈???
                                            verbose      = FALSE)#是否输出提示信息，默认为false
sub_microglia.markers1group_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                                 organism     = 'mmu',
                                                 nPerm        = 1000,
                                                 minGSSize    = 3,
                                                 pvalueCutoff = 1,
                                                 verbose      = FALSE)
write.csv(sub_microglia.markers1group_GSEAGO,file = "bs_table/wfx/sub_microglia.markers1group_GSEAGO.csv")
write.csv(sub_microglia.markers1group_GSEAKEGG,file = "bs_table/wfx/sub_microglia.markers1group_KEGG.csv")
#2号群
mig2_marker <- import("bs_table/wfx/mig2_marker.csv")
mig2_marker <- mig2_marker[,c(1,3)]
colnames(mig2_marker) <- c("Gene","LogFC")
mig_marker2 = mig2_marker[,1]
mig_marker2 <- bitr(mig_marker2,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
mig2_marker <- merge(mig_marker2,mig2_marker,by = 1,all = F)
geneList = mig2_marker[,3]
names(geneList) = as.character(mig2_marker[,2])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
sub_microglia.markers2group_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                            OrgDb = org.Mm.eg.db,
                                            ont = "BP",#可选择bp.MF,CC,ALL
                                            nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                            minGSSize = 1,#最小基因集的基因数
                                            maxGSSize = 500,#最大基因集的基因数
                                            pvalueCutoff = 10,
                                            pAdjustMethod = "BH",#p值的阈???
                                            verbose      = FALSE)#是否输出提示信息，默认为false
sub_microglia.markers2group_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                                 organism     = 'mmu',
                                                 nPerm        = 1000,
                                                 minGSSize    = 3,
                                                 pvalueCutoff = 1,
                                                 verbose      = FALSE)
write.csv(sub_microglia.markers2group_GSEAGO,file = "bs_table/wfx/sub_microglia.markers2group_GSEAGO.csv")
write.csv(sub_microglia.markers2group_GSEAKEGG,file = "bs_table/wfx/sub_microglia.markers2group_KEGG.csv")
#4号群
mig4_marker <- import("bs_table/wfx/mig4_marker.csv")
mig4_marker <- mig4_marker[,c(1,3)]
colnames(mig4_marker) <- c("Gene","LogFC")
mig_marker4 = mig4_marker[,1]
mig_marker4 <- bitr(mig_marker4,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
mig4_marker <- merge(mig_marker4,mig4_marker,by = 1,all = F)
geneList = mig4_marker[,3]
names(geneList) = as.character(mig4_marker[,2])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
sub_microglia.markers4group_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                            OrgDb = org.Mm.eg.db,
                                            ont = "BP",#可选择bp.MF,CC,ALL
                                            nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                            minGSSize = 1,#最小基因集的基因数
                                            maxGSSize = 500,#最大基因集的基因数
                                            pvalueCutoff = 10,
                                            pAdjustMethod = "BH",#p值的阈???
                                            verbose      = FALSE)#是否输出提示信息，默认为false
sub_microglia.markers4group_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                                 organism     = 'mmu',
                                                 nPerm        = 1000,
                                                 minGSSize    = 3,
                                                 pvalueCutoff = 1,
                                                 verbose      = FALSE)
write.csv(sub_microglia.markers4group_GSEAGO,file = "bs_table/wfx/sub_microglia.markers4group_GSEAGO.csv")
write.csv(sub_microglia.markers4group_GSEAKEGG,file = "bs_table/wfx/sub_microglia.markers4group_KEGG.csv")
#5号群
mig5_marker <- import("bs_table/wfx/mig5_marker.csv")
mig5_marker <- mig5_marker[,c(1,3)]
colnames(mig5_marker) <- c("Gene","LogFC")
mig_marker5 = mig4_marker[,1]
mig_marker5 <- bitr(mig_marker5,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
mig5_marker <- merge(mig_marker5,mig5_marker,by = 1,all = F)
geneList = mig5_marker[,3]
names(geneList) = as.character(mig5_marker[,2])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
sub_microglia.markers5group_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                            OrgDb = org.Mm.eg.db,
                                            ont = "BP",#可选择bp.MF,CC,ALL
                                            nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                            minGSSize = 1,#最小基因集的基因数
                                            maxGSSize = 500,#最大基因集的基因数
                                            pvalueCutoff = 10,
                                            pAdjustMethod = "BH",#p值的阈???
                                            verbose      = FALSE)#是否输出提示信息，默认为false
sub_microglia.markers5group_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                                 organism     = 'mmu',
                                                 nPerm        = 1000,
                                                 minGSSize    = 3,
                                                 pvalueCutoff = 1,
                                                 verbose      = FALSE)
write.csv(sub_microglia.markers5group_GSEAGO,file = "bs_table/wfx/sub_microglia.markers5group_GSEAGO.csv")
write.csv(sub_microglia.markers5group_GSEAKEGG,file = "bs_table/wfx/sub_microglia.markers5group_KEGG.csv")
#6号群
mig6_marker <- import("bs_table/wfx/mig6_marker.csv")
mig6_marker <- mig6_marker[,c(1,3)]
colnames(mig6_marker) <- c("Gene","LogFC")
mig_marker6 = mig6_marker[,1]
mig_marker6 <- bitr(mig_marker6,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
mig6_marker <- merge(mig_marker6,mig6_marker,by = 1,all = F)
geneList = mig6_marker[,3]
names(geneList) = as.character(mig6_marker[,2])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
sub_microglia.markers6group_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                            OrgDb = org.Mm.eg.db,
                                            ont = "BP",#可选择bp.MF,CC,ALL
                                            nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                            minGSSize = 1,#最小基因集的基因数
                                            maxGSSize = 500,#最大基因集的基因数
                                            pvalueCutoff = 10,
                                            pAdjustMethod = "BH",#p值的阈???
                                            verbose      = TRUE)#是否输出提示信息，默认为false
sub_microglia.markers6group_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                                 organism     = 'mmu',
                                                 nPerm        = 1000,
                                                 minGSSize    = 3,
                                                 pvalueCutoff = 1,
                                                 verbose      = FALSE)
write.csv(sub_microglia.markers6group_GSEAGO,file = "bs_table/wfx/sub_microglia.markers6group_GSEAGO.csv")
write.csv(sub_microglia.markers6group_GSEAKEGG,file = "bs_table/wfx/sub_microglia.markers6group_KEGG.csv")
#总体
sub_microglia.markers <- FindMarkers(sub_microglia, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.2)
write.csv(sub_microglia.markers,file = "bs_table/wfx/sub_microglia.markers.csv")
genes_to_check = c('Nfam1', 'Cd33', 'Gal3', 'Tyrobp',"Syk","Itam","Trem2","Nlrp3","Il1b","Il1rn","Tlr4","Tpm1","Cd86","Tnf")
AverageExpression(sub_microglia,features = unique(genes_to_check),group.by = "sample")
####铁死亡亚群的相关分析####
save(sub_microglia1,sub_microglia.markers1,file = "bs_code/MG/fe_microglia.Rdata")
all_microglia.difgene <- FindMarkers(sub_microglia1, group.by="sample",
                                    ident.1 = "Pb", 
                                    ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.1)
cell.used <- rownames(sub_microglia1@meta.data[which(sub_microglia1@meta.data$mgcelltype%in% c(1)),])
cell.used <- rownames(sub_microglia1@meta.data[which(sub_microglia1@meta.data$mgcelltype%in% c(5)),])
cell.used <- rownames(sub_microglia1@meta.data[which(sub_microglia1@meta.data$mgcelltype%in% c(2)),])
length(cell.used)
fe_microglia1 <- subset(sub_microglia1, cells = cell.used)
fe_microglia5 <- subset(sub_microglia1, cells = cell.used)
fe_microglia2 <- subset(sub_microglia1, cells = cell.used)
DimPlot(fe_microglia1, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
DimPlot(fe_microglia5, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
DimPlot(fe_microglia2, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
fe_microglia.markers1 <- FindMarkers(fe_microglia1, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.1)
fe_microglia.markers5 <- FindMarkers(fe_microglia5, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.1)
ferroptosis_all <- c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                          "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                          "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                          "Cdo1","Me1","Coq2","Gch1","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                          'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                          "Phgdh","Hmgb1","Pparg")

fe_microglia.markers2 <- FindMarkers(fe_microglia2, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(fe_microglia.markers,file = "bs_table/femarker_MG/femarker_MG.csv")
write.csv(fe_microglia.markers1,file = "bs_table/femarker_MG/fe_li/femarker_MG1.csv")
write.csv(fe_microglia.markers5,file = "bs_table/femarker_MG/fe_li/femarker_MG5.csv")
##1号亚群的富集分析
library(rio)
Fe_marker1 <- import("bs_table/femarker_MG/fe_li/femarker_MG1.csv")
Fe_marker5 <- import("bs_table/femarker_MG/fe_li/femarker_MG5.csv")
Fe_marker1 <- Fe_marker1[,c(1:3)]
Fe_marker5 <- Fe_marker5[,c(1:3)]
colnames(Fe_marker1) <- c("Gene","p.value","LogFC")
colnames(Fe_marker5) <- c("Gene","p.value","LogFC")
Fe_gene = Fe_marker[,1]
Fe_gene <- bitr(Fe_gene,
                fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Mm.eg.db)
Fe_marker <- merge(Fe_marker,Fe_gene,by = 1,all = F)
geneList = Fe_marker[,2]
names(geneList) = as.character(Fe_marker[,3])
head(geneList)
geneList = sort(geneList, decreasing = TRUE)
fe_microglia_GSEAGO <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                            OrgDb = org.Mm.eg.db,
                                            ont = "BP",#可选择bp.MF,CC,ALL
                                            nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                            minGSSize = 1,#最小基因集的基因数
                                            maxGSSize = 500,#最大基因集的基因数
                                            pvalueCutoff = 10,
                                            pAdjustMethod = "BH",#p值的阈???
                                            verbose      = TRUE)#是否输出提示信息，默认为false
fe_microglia_GSEAKEGG <-  gseKEGG(geneList     = geneList,
                                                 organism     = 'mmu',
                                                 nPerm        = 1000,
                                                 minGSSize    = 3,
                                                 pvalueCutoff = 1,
                                                 verbose      = TRUE)
write.csv(fe_microglia_GSEAGO,file = "bs_table/femarker_MG/fe_microglia_GSEAGO.csv")
write.csv(fe_microglia_GSEAKEGG,file = "bs_table/femarker_MG/fe_microglia_GSEAKEGG.csv")
##铁死亡数据库
library(Hmisc)
ferroptosis_all <- c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
"Tf","Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Flvcr","Abcb7","Abcb8","Fxn",
"Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Mrp1",
"Cdo1","Me1","Coq2","Vkh2","Gch1","Bh4","Slc25a39",
'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Lpcat3",
"Chac1")
ferroptosis_driver <- import("bs_table/femarker_MG/ferroptosis_driver.csv")
ferroptosis_marker <- import("bs_table/femarker_MG/ferroptosis_marker.csv")
ferroptosis_suppressor <- import("bs_table/femarker_MG/ferroptosis_suppressor.csv")
ferroptosis_unclassified <- import("bs_table/femarker_MG/ferroptosis_unclassified.csv")
ferroptosis_driver <- ferroptosis_driver$symbol
ferroptosis_marker <- ferroptosis_marker$symbol
ferroptosis_suppressor <- ferroptosis_suppressor$symbol
ferroptosis_unclassified <- ferroptosis_unclassified$symbol
ferroptosis_all <- c(ferroptosis_driver,ferroptosis_marker,ferroptosis_suppressor)
ferroptosis_all <- ferroptosis_all[!duplicated(ferroptosis_all)]
ferroptosis_all <- tolower(ferroptosis_all)
ferroptosis_all <- capitalize(ferroptosis_all)
ferroptosis_all <- as.data.frame(ferroptosis_all)
write.csv(ferroptosis_all,file = "bs_table/femarker_MG/ferroptosis_all.csv")
ferroptosis_all <- import("bs_table/femarker_MG/ferroptosis_all.csv")
diff_as_gene <- import("bs_table/AS/sub_astrocyte.diffgeneall1.csv")
ferroptosis_as_gene <- merge(diff_as_gene,ferroptosis_all,by.x = "Gene",by.y = "symbol")
write.csv(ferroptosis_as_gene,file = "bs_table/AS/ferroptosis_as_gene.csv")
library(dplyr)
fe_used <- inner_join(Fe_marker,ferroptosis_all,by = c("Gene"= "symbol"))
#与数据集铁死亡基因进行对比
load("bs_code/MG/fe_microglia.Rdata")
ferroptosis_1 <- intersect(Fe_marker1$Gene,ferroptosis_all)
ferroptosis_5 <- intersect(Fe_marker5$Gene,ferroptosis_all)
ferroptosis_5 <- as.data.frame(ferroptosis_5)
names(ferroptosis_5) <- "Gene"
ferroptosis_5_dat <- merge(ferroptosis_5,Fe_marker5,by = 1,all = F)
ferroptosis_1 <- as.data.frame(ferroptosis_1)
names(ferroptosis_1) <- "Gene"
ferroptosis_1_dat <- merge(ferroptosis_1,Fe_marker1,by = 1,all = F)
ferroptosis_1vs7 <- intersect(ferroptosis_1$Gene,ferroptosis_7$Gene )
ferroptosis_1vs7 <- as.data.frame(ferroptosis_1vs7)
names(ferroptosis_1vs7) <- "Gene"
ferroptosis_1vs7_dat1 <- merge(ferroptosis_1vs7,ferroptosis_1_dat,by = 1,all = F)
ferroptosis_1vs7_dat1 <- ferroptosis_1vs7_dat1[,c(1:3)]
names(ferroptosis_1vs7_dat1) <- c("Gene","p.vaule1","LogFC1")
ferroptosis_1vs7_dat7 <- merge(ferroptosis_1vs7,ferroptosis_7_dat,by = 1,all = F)
ferroptosis_1vs7_dat7 <- ferroptosis_1vs7_dat7[,c(1:3)]
names(ferroptosis_1vs7_dat7) <- c("Gene","p.vaule7","LogFC7")
ferroptosis_1vs7_dat <- merge(ferroptosis_1vs7_dat7,ferroptosis_1vs7_dat1,by = 1,all = F)
write.csv(ferroptosis_1vs7_dat,file = "bs_table/femarker_MG/ferroptosis_1vs7_dat.csv")
write.csv(ferroptosis_1_dat,file = "bs_table/femarker_MG/fe_li/ferroptosis_1_dat.csv")
write.csv(ferroptosis_5_dat,file = "bs_table/femarker_MG/fe_li/ferroptosis_5_dat.csv")
genes_to_check = c("Slc7a11",	"Fth1",	"Hmox1","Ftmt", 'Slc25a28',"Gpx4")
genes_to_check = c("Dhodh",	"Coq2",	"Coq4","Atl1", 'Acsl4',"Acsl1","Slc7a11","Trem2","Syk","Tyrobp")
genes_to_check = c("Dhodh",	"Coq2",	"Fps1","mt-Co2","Uqcrh","mt-Cytb","mt-Nd1")
AverageExpression(fe_microglia1,features = unique(genes_to_check),group.by = "sample")
AverageExpression(sub_microglia1,features = unique(genes_to_check),group.by = "sample")
VlnPlot(sub_microglia1,group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
genes_to_check = c("Dhodh",	"Coq2",	"Acsl4","Pebp1","Trem2","mt-Cytb","mt-Nd1")
genes_to_check = c("Fam210b", "Fth1", "Ftmt", "Slc25a28", "Slc25a37","Nfs1","Coq2","Dhodh",
                  "Coq10h2","Aifm2", "Aldh3a2", "Aldh7a1", "Aldh9a1", "Cat", "Glrx2", "Gpx1", "Gpx4", 
                  "Gsr", "Hagh", "Mgst1", "Msra", "Msrb2",  "Nit1", "Oxr1", "Prdx2", "Prdx3", "Prdx4", 
                  "Prdx5", "Prdx6",  "Txn2", "Txnrd1", "Txnrd2","Hmox1","Cisd1","Cisd2",
                  "Nox4","Nox3","Nox1","Nox5","Cdsd1","Flvcr1b","Vdac2","Vdac3")

genes_to_check = c("Slc40a1","Slc11a2","Tfrc","Dhodh","Trf","Cisd1","Vdac3","Tmem173","Nox4","Aco1","Vdac2","Idh2")
m_Fegene <- AverageExpression(sub_astrocyte,features = unique(genes_to_check),group.by = "sample")
m_Fegene <- m_Fegene$RNA
write.csv(m_Fegene,file = "bs_table/femarker_MG/mferroptosis_all.csv")
genes_to_check = c("Aifm2", "Aldh3a2", "Aldh7a1", "Aldh9a1", "Cat", "Glrx2", "Gpx1", "Gpx4", "Gsr", "Hagh", 
"Mgst1", "Msra", "Msrb2",  "Nit1", "Oxr1", "Prdx2", "Prdx3", "Prdx4", "Prdx5", "Prdx6", "Sod1", 
"Sod2", "Txn2", "Txnrd1", "Txnrd2","Nox4")

####mental response####
library(Hmisc)
mental <- import("bs_table/AS/mental response.csv")
mental  = mental[!duplicated(mental$V1),]
mental <- mental$V1
mental <- toupper(mental)
mental <- tolower(mental)
mental <- capitalize(mental)
cell.used <- rownames(sce.all@meta.data[which(sce.all@meta.data$celltype%in% c("Astrocyte","Microglia","Oligodendrocyte",
                                                                               "Neuron","Epithelial cells","Endothelial cells",
                                                                               "SMC","Neural stem cell")),])
sce.all <- subset(sce.all, cells = cell.used)
mental_socre <- list(mental)

sce.all <- AddModuleScore(sce.all,features = mental_socre,ctrl = 100,name = "mental_socre")
p <- ggboxplot(sce.all@meta.data, x="celltype", y="mental_socre1", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#E7B800", "#00AFBB"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")
compare_means(mental_socre1 ~ sample, data = sce.all@meta.data,group.by = "celltype")
p + stat_compare_means(aes(sample = sample),label = "p.signif")
a <- GetAssayData(object = sce.all,slot = "counts")
####铁死亡打分####
cell.used <- rownames(sce.all@meta.data[which(sce.all@meta.data$celltype%in% c("Astrocyte","Microglia","Oligodendrocyte",
                                                                               "Neuron","Epithelial cells","Endothelial cells",
                                                                               "SMC","Neural stem cell")),])
sce.all <- subset(sce.all, cells = cell.used)
ferroptosis_all <- list(c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                     "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                     "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                     "Cdo1","Me1","Coq2","Gch1","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                     'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                     "Phgdh","Hmgb1","Pparg"))
WNT_features <- ferroptosis_all
sce.all <- AddModuleScore(sce.all,features = WNT_features,ctrl = 100,name = "WNT_features")
p <- ggboxplot(sce.all@meta.data, x="celltype", y="WNT_features1", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#E7B800", "#00AFBB"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")
compare_means(WNT_features1 ~ sample, data = sce.all@meta.data,group.by = "celltype")
p + stat_compare_means(aes(sample = sample),label = "p.signif")
BiocManager::install("AUCell",force = TRUE)
install.packages("Matrix")
remove.packages("Matrix")
library(AUCell)
library(Matrix)
WNT_features <- list(c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                       "Tf","Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Flvcr","Abcb7","Abcb8","Fxn",
                       "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Mrp1",
                       "Cdo1","Me1","Coq2","Vkh2","Gch1","Bh4","Slc25a39",
                       'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Lpcat3",
                       "Chac1" ))

sub_microglia1 <- AddModuleScore(sub_microglia1,
                       features = WNT_features,
                       ctrl = 100,
                       name = "WNT_features")
head(sub_microglia1@meta.data)

sce.all <- AddModuleScore(sce.all,
                                 features = WNT_features,
                                 ctrl = 100,
                                 name = "Ferroptosis_Score")
#这里就得到了基因集评分结果，但是注意列名为 WNT_features1
colnames(sub_microglia1@meta.data)[15] <- 'mgcelltype'
VlnPlot(sub_microglia1,features = 'WNT_features1', 
        pt.size = 0.2,group.by = "mgcelltype")
ggboxplot(sub_microglia1@meta.data, x="mgcelltype", y="WNT_features1", width = 0.6, 
          color = "black",#轮廓颜色
          fill = "mgcelltype",
          palette = colour,
          #palette =c("#E7B800", "#00AFBB"),#分组着色
          xlab = F, #不显示x轴的标签
          x.text.angle = 90,
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right")
colnames(sce.all@meta.data)[11] <- 'Ferroptosis_Score'
library(ggpubr)
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p=VlnPlot(sub_microglia1,'WNT_features1')
df = aggregate(sub_microglia1$WNT_features1,list(sub_microglia1$mgcelltype),median)
ggbarplot(df,'Group.1','x') + coord_flip()
p=VlnPlot(sce.all,'Ferroptosis_Score1')

df = aggregate(sce.all$sample,list(sce.all$Ferroptosis_Score1),median)

ggbarplot(df,'Group.1','x') + coord_flip()
p <- ggboxplot(sub_microglia1@meta.data, x="mgcelltype", y="WNT_features1", width = 0.6, 
          color = "sample",#轮廓颜色
          palette =c("#E7B800", "#00AFBB"),#分组着色
          xlab = F, #不显示x轴的标签
          x.text.angle = 90,
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right")
ggboxplot(sce.all@meta.data, x="celltype", y="WNT_features1", width = 0.6, 
               color = "black",#轮廓颜色
               fill = "celltype",
               palette = colour,
               #palette =c("#E7B800", "#00AFBB"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")

compare_means(WNT_features1 ~ sample, data = sub_microglia1@meta.data,group.by = "mgcelltype")
p + stat_compare_means(aes(sample = sample),label = "p.signif")

mean(sub_microglia1$Ferroptosis_Score1,)
#提取做分析
mg = sce.all@meta.data[sce.all$celltype%in% c("Microglia"),]
fe_socure <- data.frame( group = mg$sample,
                  socure = mg$Ferroptosis_Score1)
df = aggregate(sub_microglia1$WNT_features1,list(sub_microglia1$sample),median)
str(fe_socure)
fe_socure$group <- as.factor(fe_socure$group)
oneway <- aov(socure~group, data = fe_socure)
anova(oneway)
install.packages("agricolae")
library(agricolae)
TukeyHSD(oneway)
out <- LSD.test(oneway,"fe_socure")
out
####铁死亡1号群和7号群合并分析####
#名字和群匹配 
sub_microglia2 <- sub_microglia1
celltype=data.frame(ClusterID=0:7 ,
                    mgcelltype1= 0:7) 
celltype[celltype$ClusterID %in% c(0),2]='0' 
celltype[celltype$ClusterID %in% c(1,7),2]='1'  
celltype[celltype$ClusterID %in% c(2),2]='2'  
celltype[celltype$ClusterID %in% c(3),2]='3'  
celltype[celltype$ClusterID %in% c(4),2]='4'   
celltype[celltype$ClusterID %in% c(5),2]='5'
celltype[celltype$ClusterID %in% c(6),2]='6'
head(celltype)
celltype
table(celltype$mgcelltype1)
sub_microglia2@meta.data$mgcelltype1 = "NA"

for(i in 1:nrow(celltype)){
  sub_microglia2@meta.data[which(sub_microglia2@meta.data$mgcelltype == celltype$ClusterID[i]),'mgcelltype1'] <- celltype$mgcelltype1[i]}
table(sub_microglia2@meta.data$mgcelltype1)
p <- DimPlot(sub_microglia2, reduction = "tsne",group.by="mgcelltype1",
             split.by = "sample",label = TRUE, pt.size = 0.3)
p
##计算比例
sub_microglia2@active.ident <- plyr::mapvalues(x = sub_microglia2@active.ident, 
                                               from = celltype$ClusterID, to = celltype$mgcelltype1)
sub_microglia2@active.ident<-factor(sub_microglia2@active.ident,levels = unique(celltype$mgcelltype1))
table(sub_microglia2$sample)
prop.table(table(Idents(sub_microglia2)))
table(Idents(sub_microglia2), sub_microglia2$sample)
Cellratio <- prop.table(table(Idents(sub_microglia2), sub_microglia2$sample), margin = 2)
Cellratio
ferroptosis_all <- c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                     "Tf","Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Flvcr","Abcb7","Abcb8","Fxn",
                     "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Mrp1",
                     "Cdo1","Me1","Coq2","Vkh2","Gch1","Bh4","Slc25a39",
                     'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Lpcat3",
                     "Chac1")
ferroptosis_all <- as.data.frame(ferroptosis_all)
cell.used <- rownames(sub_microglia2@meta.data[which(sub_microglia2@meta.data$mgcelltype1%in% c(1)),])
length(cell.used)
fe_microglia1 <- subset(sub_microglia2, cells = cell.used)
DimPlot(fe_microglia1, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
fe_microglia.markers1 <- FindMarkers(fe_microglia1, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(fe_microglia.markers1,file = "bs_table/femarker_MG/femarker_MG1.1.csv")

##1号亚群的富集分析
library(rio)
Fe_marker1 <- import("bs_table/femarker_MG/femarker_MG1.1.csv")
Fe_marker1 <- Fe_marker1[,c(1:3)]
colnames(Fe_marker1) <- c("Gene","p.value","LogFC")
#与数据集铁死亡基因进行对比
ferroptosis_1 <- intersect(Fe_marker1$Gene,ferroptosis_all$ferroptosis_all)
ferroptosis_7 <- as.data.frame(ferroptosis_7)
names(ferroptosis_7) <- "Gene"
ferroptosis_7_dat <- merge(ferroptosis_7,Fe_marker7,by = 1,all = F)
ferroptosis_1 <- as.data.frame(ferroptosis_1)
names(ferroptosis_1) <- "Gene"
ferroptosis_1_dat <- merge(ferroptosis_1,Fe_marker1,by = 1,all = F)
ferroptosis_1vs7 <- intersect(ferroptosis_1$Gene,ferroptosis_7$Gene )
ferroptosis_1vs7 <- as.data.frame(ferroptosis_1vs7)
names(ferroptosis_1vs7) <- "Gene"
ferroptosis_1vs7_dat1 <- merge(ferroptosis_1vs7,ferroptosis_1_dat,by = 1,all = F)
ferroptosis_1vs7_dat1 <- ferroptosis_1vs7_dat1[,c(1:3)]
names(ferroptosis_1vs7_dat1) <- c("Gene","p.vaule1","LogFC1")
ferroptosis_1vs7_dat7 <- merge(ferroptosis_1vs7,ferroptosis_7_dat,by = 1,all = F)
ferroptosis_1vs7_dat7 <- ferroptosis_1vs7_dat7[,c(1:3)]
names(ferroptosis_1vs7_dat7) <- c("Gene","p.vaule7","LogFC7")
ferroptosis_1vs7_dat <- merge(ferroptosis_1vs7_dat7,ferroptosis_1vs7_dat1,by = 1,all = F)
write.csv(ferroptosis_1vs7_dat,file = "bs_table/femarker_MG/ferroptosis_1vs7_dat.csv")
write.csv(ferroptosis_1_dat,file = "bs_table/femarker_MG/ferroptosis_1_dat.csv")
write.csv(ferroptosis_7_dat,file = "bs_table/femarker_MG/ferroptosis_7_dat.csv")
genes_to_check = c("Slc7a11",	"Fth1",	"Hmox1","Ftmt", 'Slc25a28',"Gpx4")
genes_to_check = c("Dhodh",	"Coq2",	"Coq4","Atl1", 'Acsl4',"Acsl1","Slc7a11","Trem2","Syk","Tyrobp")
genes_to_check = c("Dhodh",	"Coq2",	"Fps1","mt-Co2","Uqcrh","mt-Cytb","mt-Nd1")
AverageExpression(fe_microglia1,features = unique(genes_to_check),group.by = "sample")
AverageExpression(sub_microglia1,features = unique(genes_to_check),group.by = "sample")
VlnPlot(sub_microglia1,group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
genes_to_check = c("Dhodh",	"Coq2",	"Acsl4","Pebp1","Trem2","mt-Cytb","mt-Nd1")
genes_to_check = c("Fam210b", "Fth1", "Ftmt", "Slc25a28", "Slc25a37","Nfs1","Coq2","Dhodh",
                   "Coq10h2","Aifm2", "Aldh3a2", "Aldh7a1", "Aldh9a1", "Cat", "Glrx2", "Gpx1", "Gpx4", 
                   "Gsr", "Hagh", "Mgst1", "Msra", "Msrb2",  "Nit1", "Oxr1", "Prdx2", "Prdx3", "Prdx4", 
                   "Prdx5", "Prdx6",  "Txn2", "Txnrd1", "Txnrd2","Hmox1","Cisd1","Cisd2",
                   "Nox4","Nox3","Nox1","Nox5","Cdsd1","Flvcr1b","Vdac2","Vdac3")

m_Fegene <- AverageExpression(fe_microglia7,features = unique(genes_to_check),group.by = "sample")
m_Fegene <- m_Fegene$RNA
write.csv(m_Fegene,file = "bs_table/femarker_MG/mferroptosis_all.csv")
genes_to_check = c("Aifm2", "Aldh3a2", "Aldh7a1", "Aldh9a1", "Cat", "Glrx2", "Gpx1", "Gpx4", "Gsr", "Hagh", 
                   "Mgst1", "Msra", "Msrb2",  "Nit1", "Oxr1", "Prdx2", "Prdx3", "Prdx4", "Prdx5", "Prdx6", "Sod1", 
                   "Sod2", "Txn2", "Txnrd1", "Txnrd2","Nox4")
####星胶####
library(ggplot2) 
library(cowplot) 
library(paletteer)  
library(gplots)
library(ggpubr)    
library(ggsci) 
library(stringr)
library(rio)
library(irGSEA)
# install packages from Bioconductor
bioconductor.packages <- c("AUCell", "BiocParallel", "ComplexHeatmap", 
                           "decoupleR", "fgsea", "ggtree", "GSEABase", 
                           "GSVA", "Nebulosa", "scde", "singscore",
                           "SummarizedExperiment", "UCell",
                           "viper","sparseMatrixStats")

for (i in bioconductor.packages) {
  if (!requireNamespace(i, quietly = TRUE)) {
    BiocManager::install(i, ask = F, update = F)
  }
}

# install packages from Github
if (!requireNamespace("irGSEA", quietly = TRUE)) { 
  devtools::install_github("chuiqin/irGSEA", force =T)
}
install.packages('ks')
sub_astrocyte1 <- irGSEA.score(object =sub_astrocyte, assay = "RNA", 
                             slot = "data", seeds = 123, ncores = 4,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T, 
                             species = "Mus musculus", category = "H",  
                             subcategory = NULL, geneid = "symbol",
                             method = c("singscore"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                             kcdf = 'Gaussian')

result.dge <- irGSEA.integrate(object = sub_astrocyte1,
                               group.by = "seurat_clusters",
                               method = c("singscore"))
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "RRA",
                                      top = 50,
                                      show.geneset = NULL)

scatterplot <- irGSEA.density.scatterplot(object = sub_astrocyte1,
                                          method = "singscore",
                                          show.geneset = "HALLMARK-E2F-TARGETS",
                                          reduction = "umap")

ferroptosis_all <- list(c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                          "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                          "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                          "Cdo1","Me1","Coq2","Gch1","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                          'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                          "Phgdh","Hmgb1","Pparg","Prdx3"))

sub_astrocyte1 <- irGSEA.score(object = sub_astrocyte, assay = "RNA", slot = "data", 
                             seeds = 123, ncores = 1,msigdb=F, custom = T, geneset = markers, method = c( "singscore"), kcdf = 'Gaussian')
result.dge <- irGSEA.integrate(object = sub_astrocyte1,
                               group.by = "seurat_clusters",
                               method = c("singscore"))

scatterplot <- irGSEA.density.scatterplot(object = sub_astrocyte1,
                                          method = "singscore",
                                          show.geneset = "ferroptosis-all",
                                          reduction = "umap")

densityheatmap <- irGSEA.densityheatmap(object = sub_astrocyte1,
                                        method = "singscore",
                                        show.geneset = "ferroptosis-all")
densityheatmap


ridgeplot <- irGSEA.ridgeplot(object = sub_astrocyte1,
                              method = "singscore",
                              show.geneset = "ferroptosis-all")
ridgeplot
markers <- list()
ferroptosis_all <- c("Fth1","Slc40a1","Tfrc","Slc11a2","Vdac2","Vdac3","Cisd1","Slc7a11","Gpx4","Slc3a2","Tfr","Gclc","Atf3","Atf4","Gch1","Slc3a2","Hif1a",
                     "Acsl1","Acsl3","Cs","Hmgcr","Nupr1")




p <- ggboxplot(sub_astrocyte1@meta.data, x="seurat_clusters", y="nCount_singscore", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#1E90FF", "#FF6347"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")

compare_means(nCount_singscore ~ sample, data = sub_astrocyte1@meta.data,group.by = "seurat_clusters")
p + stat_compare_means(aes(sample = sample),label = "p.signif")

load("bs_code/bs_totalcell.Rdata")
cell.used <- rownames(sce.all@meta.data[which(sce.all@meta.data$celltypes%in% c(1,2,5,11)),])
load(file = "bs_code/AS/bs_AS.Rdata")
load(file = "bs_code/AS/bs_AS_gj.Rdata")
load(file = "bs_code/AS/bs_AS_cx.Rdata")
##计算比例
table(sub_astrocyte$sample)
prop.table(table(Idents(sub_astrocyte)))
table(Idents(sub_astrocyte), sub_astrocyte$sample)
Cellratio <- prop.table(table(Idents(sub_astrocyte), sub_astrocyte$sample), margin = 2)
Cellratio
Cellratio <- prop.table(table(mean(sub_astrocyte$WNT_features1), sub_astrocyte$sample), margin = 2)
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$sample%in% c("Control")),])
sub_astrocyte_control <- subset(sub_astrocyte, cells = cell.used)
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$sample%in% c("Pb")),])
sub_astrocyte_Pb <- subset(sub_astrocyte, cells = cell.used)
aggregate(sub_astrocyte$WNT_features1,list(sub_astrocyte$sample),mean)
aggregate(sub_astrocyte_Pb$WNT_features1,list(sub_astrocyte_Pb$seurat_clusters),mean)

Cellratio
Cellratio <- as.data.frame(Cellratio)
sub_astrocyte.markers <- FindAllMarkers(sub_astrocyte,
                                        only.pos = TRUE,
                                        logfc.threshold = 0.25)# 只返回positive基因)
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(0,1,2,3,4,5,8)),])
sub_astrocyte <- subset(sub_astrocyte, cells = cell.used)
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(0,1,2,3,4,5,6,8)),])
write.csv(sub_astrocyte.markers,file = "bs_table/AS/sub_astrocyte.marker1.csv")
save(sub_astrocyte, sub_astrocyte.markers,file = "bs_code/AS/bs_AS_cx.Rdata")

top_astrocyte.markers = sub_astrocyte.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 

sub_astrocyte$seurat_clusters <- factor(x = PBMC$seurat_clusters, levels = c('0','1','2','3','4','5','6'))
library(RColorBrewer)
library(paletteer)
DoHeatmap(sub_astrocyte,
            features = as.character(unique(top_astrocyte.markers$gene)),
            group.by = "seurat_clusters",
            assay = "RNA",
            group.colors = c("#C77CFF","#7CAE00","#00BFC4","#F8766D","#AB82FF","#90EE90","#00CD00","#008B8B","#FFA500"))+
scale_fill_gradientn(colors = c("navy","white","firebrick3"))

sub_astrocyte <- subset(sub_astrocyte, cells = cell.used)
p <- DimPlot(sub_astrocyte, reduction = 'umap', label = TRUE, split.by = "sample", pt.size = 0.6)
p <- DimPlot(sce.all, reduction = "umap",group.by="celltype",label = TRUE, pt.size = 0.3)
plot4 = p + scale_color_npg()
plot4
astro.top.markers = sub_astrocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 


cell_type_cols <- c(brewer.pal(9, "Set2"), "#FF34B3", "#BC8F8F", "#20B2AA", "#00F5FF", 
                    "#FFA500", "#ADFF2F", "#FF6A6A", "#7FFFD4", "#AB82FF", "#90EE90", 
                    "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B",
                    "#FF3030", "#7CFC00", "#000000", "#708090")
DimPlot(sub_astrocyte, label = T, cols = cell_type_cols,  pt.size = 0.6, repel = T, split.by = "sample")
astro.top.markers = sub_astrocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 



##保存分析
save(sub_astrocyte, file = "bs_code/AS/bs_AS_gj.Rdata")
genes_to_check = c("Fth1","Vdac2")
VlnPlot(sub_astrocyte,group.by = "sample",
        features = unique(genes_to_check),ncol = 2)
group6 <- AverageExpression(sub_microglia,features = unique(genes_to_check),group.by = "sample")
###降为分析
sub_astrocyte <- subset(sce.all, cells = cell.used)
sub_astrocyte <- NormalizeData(sub_astrocyte, normalization.method = "LogNormalize", scale.factor = 1e4) 
sub_astrocyte <- FindVariableFeatures(sub_astrocyte, selection.method = 'vst', nfeatures = 2000)
sub_astrocyte <- ScaleData(sub_astrocyte, vars.to.regress = "mito.percent")
sub_astrocyte <- RunPCA(sub_astrocyte, features = VariableFeatures(object = sub_astrocyte))
sub_astrocyte <- FindNeighbors(sub_astrocyte, dims = 1:10)
sub_astrocyte <- FindClusters(sub_astrocyte, resolution = 0.5)
sub_astrocyte <- FindClusters(sub_astrocyte, resolution = 0.7)
sub_astrocyte <- FindClusters(sub_astrocyte, resolution = 0.6)
head(Idents(sub_astrocyte), 5)
table(sub_astrocyte$seurat_clusters)
sub_astrocyte <- RunUMAP(sub_astrocyte, dims = 1:10)
sub_astrocyte <- RunTSNE(sub_astrocyte, dims = 1:10)
DimPlot(sub_astrocyte, reduction = 'tsne', label = TRUE, pt.size = 0.6)
DimPlot(sub_astrocyte, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
DimPlot(sub_astrocyte, reduction = 'umap', label = TRUE, split.by = "sample", pt.size = 0.6)
genes_to_check = c('Gpx4','Slc40a1', 'Fth1',"Hmox1","Tfrc","Nfe2l2", "Vdac2","Ptgds","Mdk","Acsl4")
genes_to_check = c('Slc1a3',"Slc1a2","Cadm2","Ltih3","Mgst1","Slc6a11","Tfrc","Mdk","Ftl1","B2m",'Gfap',"Vim","App","Nupr1","Agt","Aqp4")
FeaturePlot(sub_astrocyte, pt.size = 0.5,features = c("Mdk"),split.by = "sample")
DotPlot(sub_astrocyte, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()
genes_to_check = c('Gfap',"Cd44","Osmr","Chi3l1","C3","S100a10","Steap4","Lcn2","Cxcl10","Fth1","Nfkb1","Il10")
DotPlot(sub_astrocyte, group.by = 'seurat_clusters',
        features = unique(genes_to_check)) + RotatedAxis()
##保存
save(sub_astrocyte,sub_astrocyte.markers, file = "bs_code/AS/bs_AS.Rdata")
load(file = "bs_code/AS/bs_AS.Rdata")
##all
sub_astrocyte.diffgene <- FindMarkers(sub_astrocyte, group.by="sample",
                                       ident.1 = "Pb", 
                                       ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(sub_astrocyte.diffgene,file = "bs_table/AS/sub_astrocyte.diffgene1.csv")
Fe_marker <- import(file = "bs_table/AS/sub_astrocyte.diffgene1.csv")
ferroptosis_all <- c("Fth1","Ftl","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                     "Tf","Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Flvcr","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                     "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Mrp1","Atf3","Atf4",
                     "Cdo1","Me1","Coq2","Vkh2","Gch1","Bh4","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                     'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                     "Phgdh","Hmgb1","Pparg")
AS_ferroptosis_gene <- intersect(Fe_marker$Gene,ferroptosis_all)
AS_ferroptosis_gene <- as.data.frame(AS_ferroptosis_gene)
names(AS_ferroptosis_gene) <- "Gene"
AS_ferroptosis_dat <- merge(AS_ferroptosis_gene,Fe_marker,by = 1,all = F)
yhlsh <- c("mt-Nd1","mt-Nd2","mt-Nd3","mt-Co1","mt-Co3","mt-Atp6","mt-Cytb","mt-Co2","mt-Atp8")
VlnPlot(subset(sub_astrocyte, downsample = 100),group.by = "sample",
        features = unique(yhlsh), ncol = 5,assay = "RNA",cols = c("#6495ED","#FF6A6A"))
##星胶3号群
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(3)),])
length(cell.used)
sub_astrocyte3 <- subset(sub_astrocyte, cells = cell.used)
sub_astrocyte.diffgene3 <- FindMarkers(sub_astrocyte3, group.by="sample",
                                      ident.1 = "Pb", 
                                      ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.15)
write.csv(sub_astrocyte.diffgene3,file = "bs_table/AS/sub_astrocyte.diffgene3.csv")
Fe_marker3 <- import(file = "bs_table/AS/sub_astrocyte.diffgene3.csv")
AS_ferroptosis_gene <- intersect(Fe_marker3$V1,ferroptosis_all)
AS_ferroptosis_gene <- as.data.frame(AS_ferroptosis_gene)
names(AS_ferroptosis_gene) <- "Gene"
AS_ferroptosis_dat3 <- merge(AS_ferroptosis_gene,Fe_marker3,by = 1,all = F)
##星胶5号群
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(4)),])
length(cell.used)
sub_astrocyte4 <- subset(sub_astrocyte, cells = cell.used)
sub_astrocyte.diffgene4 <- FindMarkers(sub_astrocyte4, group.by="sample",
                                       ident.1 = "Pb", 
                                       ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.15)
write.csv(sub_astrocyte.diffgene4,file = "bs_table/AS/sub_astrocyte.diffgene4.csv")
Fe_marker3 <- import(file = "bs_table/AS/sub_astrocyte.diffgene3.csv")
AS_ferroptosis_gene <- intersect(Fe_marker3$V1,ferroptosis_all)
AS_ferroptosis_gene <- as.data.frame(AS_ferroptosis_gene)
names(AS_ferroptosis_gene) <- "Gene"
AS_ferroptosis_dat3 <- merge(AS_ferroptosis_gene,Fe_marker3,by = 1,all = F)
##星胶8号群
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(8)),])
length(cell.used)
sub_astrocyte8 <- subset(sub_astrocyte, cells = cell.used)
sub_astrocyte.diffgene8 <- FindMarkers(sub_astrocyte8, group.by="sample",
                                      ident.1 = "Pb", 
                                      ident.2 = "Control", logfc.threshold = 0.2, min.pct = 0.15)
write.csv(sub_astrocyte.diffgene8,file = "bs_table/AS/sub_astrocyte.diffgene8.csv")
Fe_marker8 <- import(file = "bs_table/AS/sub_astrocyte.diffgene8.csv")
AS_ferroptosis_gene <- intersect(Fe_marker8$V1,ferroptosis_all)
AS_ferroptosis_gene <- as.data.frame(AS_ferroptosis_gene)
names(AS_ferroptosis_gene) <- "Gene"
AS_ferroptosis_dat8 <- merge(AS_ferroptosis_gene,Fe_marker8,by = 1,all = F)
##星胶5号群
ferroptosis_all <- import("bs_table/femarker_MG/ferroptosis_all.csv")
ferroptosis_all <- ferroptosis_all$Symbols
ferroptosis_all <- list(ferroptosis_all)
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(5)),])
length(cell.used)
sub_astrocyte5 <- subset(sub_astrocyte, cells = cell.used)
DimPlot(sub_astrocyte5, reduction = 'tsne', label = TRUE, split.by = "sample", pt.size = 0.6)
sub_astrocyte.diffgene <- FindMarkers(sub_astrocyte, group.by="sample",
                                     ident.1 = "Pb", 
                                     ident.2 = "Control", logfc.threshold = 0.1, min.pct = 0.1)
ferroptosis_all <- list(c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                     "Tf","Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Flvcr","Abcb7","Abcb8","Fxn",
                     "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Mrp1",
                     "Cdo1","Me1","Coq2","Vkh2","Gch1","Bh4","Slc25a39",
                     'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Lpcat3",
                     "Chac1"))

ferroptosis_all <- list(c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                      "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                     "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                     "Cdo1","Me1","Coq2","Gch1","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                     'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                     "Phgdh","Hmgb1","Pparg","Prdx3"))
ferroptosis_all <- c("Fth1","Slc40a1","Tfrc","Slc11a2","Vdac2","Vdac1","Cisd1","Slc7a11","Gpx4","Slc3a2","Gclc","Atf3","Atf4","Gch1","Slc3a2","Hif1a",
                     "Acsl1","Acsl3","Cs","Hspa5","Nupr1","Iscu","Prdx1","Cbs","Gstm1")
WNT_features <- ferroptosis_all
sub_astrocyte <- AddModuleScore(sub_astrocyte,
                                 features = WNT_features,
                                 ctrl = 100,
                                 name = "WNT_features")
p <- ggboxplot(sub_astrocyte@meta.data, x="seurat_clusters", y="WNT_features1", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#1E90FF", "#FF6347"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")

compare_means(WNT_features1 ~ sample, data = sub_astrocyte@meta.data,group.by = "seurat_clusters")
p + stat_compare_means(aes(sample = sample),label = "p.signif")
library(viridis)
library(irGSEA)
FeaturePlot(sub_astrocyte, features = "WNT_features1",order= T, cols = viridis(256), split.by = "sample")



p <- ggboxplot(sub_astrocyte@meta.data, x="sample", y="WNT_features1", width = 0.6, 
               palette =c("#1E90FF", "#FF6347"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")

compare_means(WNT_features1 ~ sample, data = sub_astrocyte@meta.data)
p + stat_compare_means(aes(sample = sample),label = "p.signif")


df = aggregate(sub_astrocyte$sample,list(sub_astrocyte$WNT_features1),mean)
p <- ggboxplot(sub_astrocyte@meta.data, x="seurat_clusters", y="WNT_features1", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#E7B800", "#00AFBB"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")

colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",  
         "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
         "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
         "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
VlnPlot(sub_astrocyte,'WNT_features1')
ggboxplot(sub_astrocyte@meta.data, x="seurat_clusters", y="WNT_features1", width = 0.6,
          color = "black",#轮廓颜色
          fill = "seurat_clusters",
          palette = colour,
          #palette =c("#E7B800", "#00AFBB"),#分组着色
          xlab = F, #不显示x轴的标签
          x.text.angle = 90,
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=0.5, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          legend = "right")
df = aggregate(sub_astrocyte$WNT_features1,list(sub_astrocyte$seurat_clusters),mean)
ggbarplot(df,'Group.1','x') + coord_flip()
df = compare_means(WNT_features1 ~ seurat_clusters, data = sub_astrocyte@meta.data)

p + stat_compare_means(aes(sample = sample),label = "p.signif")
write.csv(sub_astrocyte.markers5,file = "bs_table/AS/sub_astrocyte.diff5.csv")
write.csv(sub_astrocyte.markers,file = "bs_table/AS/sub_astrocyte.markers.csv")
Fe_marker5 <- import("bs_table/AS/sub_astrocyte.diff5.csv")
ferroptosis_all <- c("Fth1","Ftl","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                          "Tf","Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Flvcr","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                          "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Mrp1","Atf3","Atf4",
                          "Cdo1","Me1","Coq2","Vkh2","Gch1","Bh4","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                          'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                          "Phgdh","Hmgb1","Pparg")
AS_ferroptosis_gene <- intersect(Fe_marker5$V1,ferroptosis_all)
AS_ferroptosis_gene <- as.data.frame(AS_ferroptosis_gene)
names(AS_ferroptosis_gene) <- "Gene"
AS_ferroptosis_dat <- merge(AS_ferroptosis_gene,Fe_marker5,by = 1,all = F)
write.csv(AS_ferroptosis_dat,file = "bs_table/AS/AS_ferroptosis_dat.csv")
p <- ggboxplot(sub_astrocyte@meta.data, x="seurat_clusters", y="WNT_features1", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#E7B800", "#00AFBB"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")
compare_means(WNT_features1 ~ sample, data = sub_astrocyte@meta.data,group.by = "seurat_clusters")
p + stat_compare_means(aes(sample = sample),label = "p.signif")
ferroptosis_expression <- as.data.frame(AverageExpression(sub_astrocyte,features = unique(ferroptosis_all),
                                                          group.by = "seurat_clusters"))
As = sub_astrocyte@meta.data[sce.all$celltype%in% c("Microglia"),]
fe_socure <- data.frame( group = mg$sample,
                         socure = mg$Ferroptosis_Score1)
df = aggregate(sub_astrocyte$WNT_features1,list(sub_astrocyte$sample),median)
##火山图差异基因
diff5 <- import(file = "bs_table/AS/sub_astrocyte.diff5.csv")
logFC_t=0.20
P.Value_t = 0.05
k1 = (diff5$p_val < P.Value_t)&(diff5$avg_log2FC < -logFC_t)
k2 = (diff5$p_val < P.Value_t)&(diff5$avg_log2FC > logFC_t)
diff5 <- mutate(diff5,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(diff5$change)
dat  = diff5[!duplicated(diff5$gene),]
dat = mutate(dat,row.names(dat))
p <- ggplot(data = dat, 
            aes(x = avg_log2FC, 
                y = -log10(p_val_adj))) +
  geom_point(alpha=0.8, size=3, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("Blue", "Grey","Red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()+
  theme(panel.grid=element_blank())

library(dplyr)
library(ggplot2)

dat <- diff
p <- ggplot(data = dat, 
            aes(x = avg_log2FC,y = row.names(dat)) )+
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("ID")+  xlab("LogFC") +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8)+ 
  coord_flip()
p


library(ggrepel)

sub_astrocyte.diffgene <- FindMarkers(sub_astrocyte, group.by="sample",
                                      ident.1 = "Pb", 
                                      ident.2 = "Control", logfc.threshold = 0, min.pct = 0,pseaudocount.use = 0.01)
diff <- import("bs_table/AS/sub_astrocyte.diffgene.csv")
logFC_t=0.13
P.Value_t = 0.05
k1 = (diff$p_val < P.Value_t)&(diff$avg_log2FC < -logFC_t)
k2 = (diff$p_val < P.Value_t)&(diff$avg_log2FC > logFC_t)
diff <- mutate(diff,change = ifelse(k1,"down",ifelse(k2,"up","stable")))


diff$Difference = diff$pct.1 - diff$pct.2
ggplot(diff, aes(x=Difference, y=avg_log2FC)) + 
  geom_point(size=3, aes(color = change)) + 
  scale_color_manual(values=c( "#8C9EFF","grey","#FF5722") ) +
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()


####亚群相关性####
#比例与基因表达
rm(list = ls())
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(data.table)
library(tidyverse)
library(paletteer) 
library(vcd)
library(limma)
library(gplots)
load(file = 'first_sce3.Rdata') ####Seurat对象
table(sub_astrocyte$orig.ident)
sce.all=sub_astrocyte
table(sce.all$orig.ident,sce.all$seurat_clusters)
library(stringr)
sce.all@meta.data$orig.ident2=sce.all@meta.data$orig.ident
table(sce.all@meta.data$orig.ident2)
sce.all$orig.ident2<-sub("-",'_',sce.all$orig.ident2)

phe=str_split(sce.all$orig.ident2,'-',simplify = T)###将样本分为两组
table(sce.all$orig.ident2)
sce.all$orig.ident2=phe[,1]
sce.all$orig.ident2 <- sce.all$sample 
sce.all=sce.all[,sce.all$orig.ident2%in% c("AdjNorm_TISSUE","PDAC_TISSUE")]

#差异表达箱线图
gene=c("Atf4","Cisd1","Cisd2","Cs","Fdx1","Fth1","Gpx4","Hmgb1","Hmgcr","Nfe2l2","Phgdh",
       "Slc3a2","Slc7a11","Sqstm1","Tfrc","Vdac2","Vdac3")
sce.sub.input <- GetAssayData(sce.all, assay = "RNA", slot = "data") # normalized data matrix
sce.sub.input=sce.sub.input[rownames(sce.sub.input)%in%gene,]
sce.sub.input=as.data.frame(sce.sub.input)

sce_sub.meta <- sce.all@meta.data[,c(1,10)]#细胞类型和分组
head(sce_sub.meta)
library(devtools)  
install_github("JanCoUnchained/ggunchained")
library(ggunchained)
library(reshape2)
library(plyr)
library(reshape2)  
library(ggsci)
library(ggpubr)
library(scales)
library(tidyverse)
library(rstatix)
sce.sub.input2 <- data.frame(t(sce.sub.input))
sce.sub.input2$sample = row.names(sce.sub.input2)
sce_sub.meta$sample=row.names(sce_sub.meta)
sce_sub.meta=sce_sub.meta[rownames(sce_sub.meta)%in%rownames(sce.sub.input2),]
table(sce_sub.meta$seurat_clusters)
data_new <- merge(sce.sub.input2, sce_sub.meta,by.x = "sample",by.y ="sample")
data_new = melt(data_new)
data_new$group <- NA
table(data_new$seurat_clusters,data_new$orig.ident2)
data_new$group[which(str_detect(data_new$sample, "^C"))] <- "Control"
data_new$group[which(str_detect(data_new$sample, "^P"))] <- "Pb"
colnames(data_new) = c("sample","celltype","gene","expression", "group")

table(data_new$celltype)
cell=unique(data_new$celltype)
cell=as.data.frame(cell)

mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12)
)
table(data_new$group)
setwd('ggplot-boxplot/')
for(i in c(1:length(cell$cell))){
  #my_comparisons <- list(c("-L", "-NL"),c("-L", "-WT"),c("-NL", "-WT"))
  my_comparisons <- list(c("Control", "Pb"))
  a=cell[c(i),]
  data=data_new[data_new$celltype%in%a,]
  box <- ggplot(data, aes(x = group, y = expression))+ 
    #scale_fill_manual(values = c("#56B4E9","#E69F00"))+
    labs(y="gene expression",title = "gene")+  
    #geom_violin(aes(fill = group),trim = T,scale = "width")+
    geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
    theme_classic()+
    facet_wrap(~gene,scales = "free")+
    stat_compare_means(comparisons=my_comparisons,label = "p.signif",method = "t.test", hide.ns = F)+
    mytheme
  ggsave(filename=paste0(a,'_boxPlot-by-gene.pdf'),
         height = 25,width = 6)
}
setwd('..')


tb <- data.frame(table(sce.all@meta.data$seurat_clusters,sce.all@meta.data$orig.ident2))
tb$Var3=tb$Var2
#tb$Var3=gsub("[-,C,E,-,1,2,3,4,5,6,7,8,9,0]", "", tb$Var3)
#tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var1 == x[1],3]))

tb$Total <- apply(tb,1,function(x)sum(tb[tb$Var2 == x[2],3]))
tb<- tb %>% mutate(Percentage = round(Freq/Total,3) * 100)
table(tb$Var3,tb$Var1)
tb=tb[,c(1,4,6)]
tb$Var1=as.factor(tb$Var1)
tb$Var3=as.factor(tb$Var3)

#tb$Var3<-sub("-",'_',tb$Var3)
#phe=str_split(tb$Var3,'-',simplify = T)
#tb$Var3=phe[,1]

ggplot(tb) + 
  geom_bar(aes(x =Percentage, y=Var3 , fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Ratio',y = 'Sample')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(angle=45,hjust = 1))
ggsave(filename='percent_celltype.pdf',
       height = 5,width = 10)

tb$Var3<-sub("-",'_',tb$Var3)
phe=str_split(tb$Var3,'-',simplify = T)
tb$Var3=phe[,1]

p <- ggboxplot(tb, x = "Var3", y = "Percentage",
               color = "Var1", 
               add = "jitter",
               facet.by = "Var1", short.panel.labs = FALSE)
my_comparisons <- list(c("AdjNorm_TISSUE", "PDAC_TISSUE"))
p + stat_compare_means(comparisons=my_comparisons,method = "t.test")+
  theme(axis.text.x=element_text(angle=45,hjust = 1)) + 
  labs(fill = "Cluster")+facet_wrap(~Var1,scales = "free")+mytheme

ggsave(filename='percent_celltype2.pdf',
       height = 18,width = 10)

DimPlot(sce.all, reduction = "umap", group.by = "celltype",
        split.by = 'orig.ident2',
        label = T,pt.size = 0.1,label.size = 3,
        repel = T,label.box = T) + 
  scale_colour_manual(values = pal_d3("category20")(20),
                      aesthetics = c("colour", "fill"))
ggsave('umap-by-celltype-ggsci.pdf',height = 4,width=8)
####拟时序分析星胶####
load(file = "bs_code/AS/bs_AS_cx.Rdata")
diff_as3 <- import("bs_table/AS/sub_astrocyte.diffgene3.csv")
diff_as5 <- import("bs_table/AS/sub_astrocyte.diff5.csv")
diff_as8 <- import("bs_table/AS/sub_astrocyte.diffgene8.csv")
a <- intersect(diff_as3$V1,intersect(diff_as5$gene,diff_as8$gene))
install.packages("BiocManager")
BiocManager::install("monocle")
library(monocle)
BEAM_res <- readRDS("bs_code/BEAM_res.rds")
mycds <- readRDS("bs_code/mycds.rds")
#抽取部分细胞进行拟时序分析
allCells=names(Idents(sub_astrocyte))
allType = levels(Idents(sub_astrocyte))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sub_astrocyte)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,1200,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sub_astrocyte[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(cg_sce))
DimPlot(cg_sce, label = T, pt.size = 0.5)
cg_sce <- subset(cg_sce, idents=c(3,5))
cg_sce <- subset(cg_sce, idents=c(3,8))
cg_sce <- subset(cg_sce, idents=c(3,5,8))
scRNAsub <- cg_sce
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
# meta表转成特定格式
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# 基因名表转成特定格式
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
expressed_genes <- row.names(subset(fData(mycds)))
top10 <- sub_astrocyte.markers[sub_astrocyte.markers$cluster %in% c(3,5),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- sub_astrocyte.markers[sub_astrocyte.markers$cluster %in% c(3,8),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- sub_astrocyte.markers[sub_astrocyte.markers$cluster %in% c(3,5,8),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

diff.wilcox = FindAllMarkers(scRNAsub)
all.markers = diff.wilcox %>%dplyr::select(gene, everything()) %>% subset(p_val<0.05)
diff.genes <- subset(all.markers,p_val_adj<0.01)$gene

ordering_genes<-as.matrix(top10)
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
p1
mycds <- setOrderingFilter(mycds, ordering_genes)
p1 <- plot_ordering_genes(mycds)
p1
mycds <- reduceDimension(mycds, max_components = 4, method = 'DDRTree')
#排序
mycds <- orderCells(mycds, root_state = 4)
mycds <- orderCells(mycds)
mycds <- orderCells(mycds,reverse = T)
save()
p1 = plot_cell_trajectory(mycds, color_by = "seurat_clusters")
p1 + facet_wrap(~sample, nrow = 1)
tsw = c("Fth1")
p1 = plot_genes_in_pseudotime(mycds[tsw,],color_by = "seurat_clusters") 
p1 + facet_wrap(~sample, nrow = 1)

plot_cell_trajectory(mycds,color_by = "Pseudotime")
p1 = plot_cell_trajectory(mycds, color_by = "Pseudotime") 
p1 + scale_color_viridis_c()

cg=as.character(head(top10$gene)) 
plot_genes_in_pseudotime(mycds[diff.genes,],color_by = "seurat_clusters")
library(ggpubr)
df <- pData(mycds)
ggplot(df,aes(Pseudotime, colour = seurat_clusters, fill = seurat_clusters))+geom_density(bw=0.5,size = 1, alpha = 0.5) + theme_classic2()
ferroptosis_all <- c("Atf3","Cbs","Cirbp","Gstm1","Mgst1","Gsk3b","Hif1a","Vdac3","Vdac2",
                     "Hspa5","Sat1","Slc3a2","Fth1","Cisd1","Trf","Tfrc","Iscu","Slc40a1","Slc7a11","Ftl1",
                     "Gpx4","Prdx1","Atf4","Nupr1")
tsw <- c("Tfrc","Vdac1","Fth1","Ftl1","Aco1","Prdx1","Gpx4","Slc7a11","Slc3a2","Cisd1")
test <- top10$gene
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[ferroptosis_all,], num_clusters = c(3,5,8), show_rownames = TRUE,
                                                 return_heatmap = TRUE)
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[ferroptosis_all,], num_clusters = c(3), show_rownames = TRUE,
                                                 return_heatmap = TRUE,
                                                 hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
plot_pseudotime_heatmap(mycds[ferroptosis_all,], num_clusters = c(3), show_rownames = TRUE,
                                                 return_heatmap = TRUE,
                                                 hmcols = colorRampPalette(c("#6666FF","white","#FF9999"))(100))

BEAM_res=BEAM(mycds,branch_point = 1,progenitor_method = "duplicate",cores = 4)
saveRDS(mycds, file = "bs_code/mycds.rds")
BEAM_branch1 = BEAM_res
head(BEAM_branch1)
colnames(BEAM_branch1)

BEAM_branch1 <- BEAM_branch1[order(BEAM_branch1$qval),]

BEAM_branch1 <- BEAM_branch1[,c("gene_short_name", "pval", "qval")]
head(BEAM_branch1) 

BEAM_res = BEAM_branch1
ferroptosis_all_fe <- c("Vdac1","Vdac2","Fth1","Cisd1","Trf","Tfrc","Iscu","Slc40a1","Ftl1","Slc11a2")
ferroptosis_all_gugt <- c("Cbs","Gstm1","Mgst1","Gsk3b","Slc3a2","Slc7a11","Gpx4")
ferroptosis_all_gyh <- c("Hspa5","Prdx1","Sat1")
ferroptosis_all_tf <- c("Atf3","Hif1a","Atf4","Nupr1")
plot_genes_branched_heatmap(mycds[ferroptosis_all,],
                                 branch_point = 1,
                                 num_clusters = 1, #这些基因被分成几个group
                                 cores = 1,
                                 hmcols = colorRampPalette(c("#6666FF","white","#FF9999"))(60),
                                 use_gene_short_name = T,
                                 show_rownames = T)#是否返回一些重要信息
plot_genes_branched_heatmap(mycds[ferroptosis_all_fe,],
                            branch_point = 1,
                            num_clusters = 1, #这些基因被分成几个group
                            cores = 1,
                            hmcols = colorRampPalette(c("#6666FF","white","#FF9999"))(60),
                            use_gene_short_name = T,
                            show_rownames = T)#是否返回一些重要信息
plot_genes_branched_heatmap(mycds[ferroptosis_all_gugt,],
                            branch_point = 1,
                            num_clusters = 1, #这些基因被分成几个group
                            cores = 1,
                            hmcols = colorRampPalette(c("#6666FF","white","#FF9999"))(60),
                            use_gene_short_name = T,
                            show_rownames = T)#是否返回一些重要信息
plot_genes_branched_heatmap(mycds[ferroptosis_all_gyh,],
                            branch_point = 1,
                            num_clusters = 1, #这些基因被分成几个group
                            cores = 1,
                            hmcols = colorRampPalette(c("#6666FF","white","#FF9999"))(60),
                            use_gene_short_name = T,
                            show_rownames = T)#是否返回一些重要信息
plot_genes_branched_heatmap(mycds[ferroptosis_all_tf,],
                            branch_point = 1,
                            num_clusters = 1, #这些基因被分成几个group
                            cores = 1,
                            hmcols = colorRampPalette(c("#6666FF","white","#FF9999"))(60),
                            use_gene_short_name = T,
                            show_rownames = T)#是否返回一些重要信息
pdf("branched_heatmap.pdf",width = 5,height = 6)
tmp1$ph_res
dev.off()

plot_genes_branched_pseudotime(mycds[tsw,],
                               branch_point = 1,
                               color_by = "seurat_clusters",
                               cell_size=2,
                               ncol = 2)
####拟时序分析小胶####
install.packages("BiocManager")
BiocManager::install("monocle")
library(monocle)
#抽取部分细胞进行拟时序分析
allCells=names(Idents(sub_microglia1))
allType = levels(Idents(sub_microglia1))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sub_microglia1)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,1000,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sub_microglia1[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(cg_sce))
DimPlot(cg_sce, label = T, pt.size = 0.5)
cg_sce <- subset(cg_sce, idents=c(0,1,5))
scRNAsub <- cg_sce
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
# meta表转成特定格式
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# 基因名表转成特定格式
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
expressed_genes <- row.names(subset(fData(mycds)))
top10 <- sub_microglia.markers[sub_microglia.markers$cluster %in% c(0,1,5),] %>% 
  group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

diff.wilcox = FindAllMarkers(scRNAsub)
all.markers = diff.wilcox %>%dplyr::select(gene, everything()) %>% subset(p_val<0.05)
diff.genes <- subset(all.markers,p_val_adj<0.01)$gene

ordering_genes<-as.matrix(top10)
diff.genes <- c("Fth1","Slc40a1","Tfrc","Vdac2","Vdac3")
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
p1
mycds <- setOrderingFilter(mycds, ordering_genes)
p1 <- plot_ordering_genes(mycds)
p1
mycds <- reduceDimension(mycds, max_components = 4, method = 'DDRTree')
#排序
mycds <- orderCells(mycds, root_state = 3)
mycds <- orderCells(mycds)
mycds <- orderCells(mycds,reverse = T)
p1 = plot_cell_trajectory(mycds, color_by = "mgcelltype")
p1 + facet_wrap(~sample, nrow = 1)

plot_cell_trajectory(mycds,color_by = "Pseudotime")
p1 = plot_cell_trajectory(mycds, color_by = "Pseudotime") 
p1 + scale_color_viridis_c()

cg=as.character(head(top10$gene)) 
plot_genes_in_pseudotime(mycds[diff.genes,],color_by = "seurat_clusters")
ferroptosis_all <- c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                     "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                     "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                     "Cdo1","Me1","Coq2","Gch1","Slc25a39","Sqstm1","Hif1a",
                     'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                     "Phgdh","Hmgb1","Pparg")
ferroptosis_all <- c("Fth1","Atf3","Vdac2","Keap1","Cisd1","Nfe2l2","Fdx1","Hif1a","Hmgb1","Sqstm1","Gpx4","Pebp1","Phgdh","Acsl3","Coq2","Hmgcr")
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[ferroptosis_all,], num_clusters = c(3,5,8), show_rownames = TRUE,
                                                 return_heatmap = TRUE)
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[ferroptosis_all,], num_clusters = c(3), show_rownames = TRUE,
                                                 return_heatmap = TRUE,
                                                 hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
my_pseudotime_cluster <- plot_pseudotime_heatmap(mycds[ferroptosis_all,], num_clusters = c(3), show_rownames = TRUE,
                                                 return_heatmap = TRUE,
                                                 hmcols = colorRampPalette(c("#6666FF","white","#FF9999"))(100))



BEAM_res <- BEAM(mycds, cores = 2)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(cds[row.names(BEAM_res)[1:50]], branch_point = 1, num_clusters = 3, 
                            cores=4, use_gene_short_name=TRUE, show_rownames=TRUE)
####相关性分析基因####
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(5)),])
sub_astrocyte <- subset(sub_astrocyte, cells = cell.used)
allCells=names(Idents(sub_astrocyte))
allType = levels(Idents(sub_astrocyte))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sub_astrocyte)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,1000,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sub_astrocyte[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(sub_astrocyte))
exprSet <- sub_astrocyte@assays[["RNA"]]@data
exprSet <- as.data.frame(t(exprSet))
y <- as.numeric(exprSet[,"Tfrc"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)
for(i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type = "pearson")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
View(cor_data_df)
install.packages("ggstatsplot")
install.packages("dplyr")
library(ggstatsplot)
ggscatterstats(data = exprSet,
               y = Tfrc,
               x = Hif1a,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram",
               title = "relationship")
FeatureScatter(cg_sce,feature1 = "Tfrc", feature2 = "Raft1")
FeatureScatter(sub_astrocyte,feature1 = "Tfrc", feature2 = "Hif1a")
####yq-bc1####
load("bs_code/bs_totalcell.Rdata")
#clustree找最合适的分辨率点
seq <- seq(0.1, 1, by = 0.1)
for(res in seq){
  sce.all <- FindClusters(sce.all, resolution = res)
}
p1 <- clustree(sce.all, prefix = 'RNA_snn_res.') + coord_flip()
p2 <- DimPlot(sce.all, group.by = 'RNA_snn_res.0.3', label = T)
p1 + p2 + plot_layout(widths = c(3, 1))
#看一下三维及标准化后的数据
DimHeatmap(sce.all, dims = 1:20,cell = 500)
ElbowPlot(sce.all, ndims = 30)
View(sce.all[["RNA"]]@scale.data[1:30,1:30])
#cell cycle看一下是否有细胞周期效应
sce.all<- CellCycleScoring(object = sce.all, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = sce.all@meta.data)
DimPlot(sce.all,reduction = "umap",label = F,group.by="Phase",pt.size = 0.2)
ggsave('figure/cellcycle.tiff',height = 7,width = 11)
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "mito.percent","ribo.percent"), ncol = 4)
####yq-fig1####
###fig1-j-k
library(rio)
library(tidyr)
library(paletteer)
library(corrplot)
library(psych)
##读取数据及标准化
dat_PZ <- import("yq_fig/Fig1/cor_PZ1.xlsx")
dat_HM <- import("yq_fig/Fig1/cor _ HM1.xlsx")
dat_HY <- import("yq_fig/Fig1/cor _XQ1.xlsx")
dat_ST <- import("yq_fig/Fig1/cor_WZT.xlsx")
dat_HH <- import("yq_fig/Fig1/xjx.xlsx")
#
#pz
xwx_PZ <- dat_PZ[,c("FST","SPT","CTCS")]
xwx_PZ <- na.omit(xwx_PZ)
FEP_PZ <- dat_PZ[,c("Fe","GSH","MDA","Gpx4","Ptgs2")]
FEP_PZ <- na.omit(FEP_PZ)
#HM
xwx_HM <- dat_HM[,c("FST","SPT","CTCS")]
xwx_HM <- na.omit(xwx_HM)
FEP_HM <- dat_HM[,c("Fe","GSH","MDA","Gpx4","Ptgs2")]
FEP_HM <- na.omit(FEP_HM)
#HY
xwx_HY <- dat_HY[,c("FST","SPT","CTCS")]
xwx_HY <- na.omit(xwx_HY)
FEP_HY <- dat_HY[,c("Fe","GSH","MDA","Gpx4","Ptgs2")]
FEP_HY <- na.omit(FEP_HY)


#pz
xwx_PZ <- dat_PZ[,c("FST","TST","GJB","KCJL","KCB","SPT")]
xwx_PZ <- na.omit(xwx_PZ)
FEP_PZ <- dat_PZ[,c("Fe","GSH","MDA","Gpx4","Ptgs2")]
FEP_PZ <- na.omit(FEP_PZ)
#HM
xwx_HM <- dat_HM[,c("FST","TST","GJB","KCJL","KCB","SPT")]
xwx_HM <- na.omit(xwx_HM)
FEP_HM <- dat_HM[,c("Fe","GSH","MDA","Gpx4","Ptgs2")]
FEP_HM <- na.omit(FEP_HM)
#HY
xwx_HY <- dat_HY[,c("FST","TST","GJB","KCJL","KCB","SPT")]
xwx_HY <- na.omit(xwx_HY)
FEP_HY <- dat_HY[,c("Fe","GSH","MDA","Gpx4","Ptgs2")]
FEP_HY <- na.omit(FEP_HY)
#ST
xwx_ST <- dat_ST[,c("FST","TST","GJB","KCJL","KCB","SPT")]
xwx_ST <- na.omit(xwx_ST)
FEP_ST <- dat_ST[,c("Fe","GSH","MDA","Gpx4","Ptgs2")]
FEP_ST <- na.omit(FEP_ST)
#hh
mir_HH <- dat_HH[,c(1)]

FEP_HH <- dat_HH[,c(2:7)]

#计算相关系数
auto_PZ <- cor(x = FEP_PZ,y = xwx_PZ,method = "pearson")
auto_HM <- cor(x = FEP_HM,y = xwx_HM,method = "pearson")
auto_HY <- cor(x = FEP_HY,y = xwx_HY,method = "pearson")
auto_ST <- cor(x = FEP_ST,y = xwx_ST,method = "pearson")
#psych计算相关性
testRes = corr.test(x = FEP_PZ,y = xwx_PZ,use = "complete")
auto_PZ <- testRes[["r"]]
p_PZ <- testRes[["p"]]
testRes = corr.test(x = FEP_HM,y = xwx_HM,use = "complete")
auto_HM <- testRes[["r"]]
p_HM <- testRes[["p"]]
testRes = corr.test(x = FEP_HY,y = xwx_HY,use = "complete")
auto_HY <- testRes[["r"]]
p_HY <- testRes[["p"]]
testRes = corr.test(x = FEP_ST,y = xwx_ST,use = "complete")
auto_ST <- testRes[["r"]]
p_ST <- testRes[["p"]]
testRes = corr.test(x = mir_HH,y = FEP_HH,use = "complete")
auto_HH <- testRes[["r"]]
p_HH <- testRes[["p"]]
#绘图
my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10)
corrplot(auto_PZ,
         method = "number",            # 图案形状 "square"方框,"circle"圆, "ellipse"椭圆, "number"数字, "shade"阴影花纹, "color"颜色方框, "pie饼图"
         type = "lower",               # "full", "lower", "upper"
         col= my_color, # 主体颜色
         bg = "white",                 # 背景颜色
         is.corr = T,                  # 输入的矩阵是否是相关性矩阵，如果是的话，数据范围会限制到-1到1
         add = T,                      # 是否在原来的图层上添加图形
         diag = F,                     # 是否显示主对角
         addCoefasPercent = F,         # 是否把相关性数值改为百分数
         order = "original",           # 排序方式 c("original", "AOE", "FPC", "hclust", "alphabet"), original：原始状态，alphabet：字母顺序 hclust，分层聚类顺序
         tl.pos = "n",                 # 坐标轴标签的位置'lt', 'ld', 'td', 'd' or 'n'   # 左边 d中间
         cl.pos = "n",                 # 图例位置：r：右边 b：下边 n：不显示
         na.label = "",                # 当为NA时，显示的内容
         p.mat = p_PZ,            # P值矩阵
         sig.level = 0.05,             # 当p大于sig.level时触发动作
)
#PZ
corrplot(auto_PZ,
         tl.col="black",
         method="number",
         tl.srt=45,
         addCoef.col="black",
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100), # 主体颜色
         cl.length = NULL,
         number.cex = 3,              # 相关性数字标签的字体大小
         number.font = 5,             # 相关性数字标签的字体
         insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
         pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
         pch.col = "black",           # 图案颜色
         pch.cex = 10,                 # 图案大小
         tl.cex = 2,                  # 坐标轴标签字体的大小
         cl.cex = 2,                # 图例的字体大小
         p.mat = p_PZ,            # P值矩阵
         sig.level = 0.05  )           # 当p大于sig.level时触发动作
#HM
corrplot(auto_HM,
         tl.col="black",
         method="number",
         tl.srt=45,
         addCoef.col="black",
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100), # 主体颜色
         cl.length = NULL,
         number.cex = 3,              # 相关性数字标签的字体大小
         number.font = 5,             # 相关性数字标签的字体
         insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
         pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
         pch.col = "black",           # 图案颜色
         pch.cex = 10,                 # 图案大小
         tl.cex = 2,                  # 坐标轴标签字体的大小
         cl.cex = 2,                # 图例的字体大小
         p.mat = p_HM,            # P值矩阵
         sig.level = 0.05  )           # 当p大于sig.level时触发动作
#HY
corrplot(auto_HY,
         tl.col="black",
         method="number",
         tl.srt=45,
         addCoef.col="black",
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100), # 主体颜色
         cl.length = NULL,
         number.cex = 3,              # 相关性数字标签的字体大小
         number.font = 5,             # 相关性数字标签的字体
         insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
         pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
         pch.col = "black",           # 图案颜色
         pch.cex = 10,                 # 图案大小
         tl.cex = 2,                  # 坐标轴标签字体的大小
         cl.cex = 2,                # 图例的字体大小
         p.mat = p_HY,            # P值矩阵
         sig.level = 0.05  )           # 当p大于sig.level时触发动作
#ST
corrplot(auto_ST,
         tl.col="black",
         method="number",
         tl.srt=45,
         addCoef.col="black",
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100), # 主体颜色
         cl.length = NULL,
         number.cex = 3,              # 相关性数字标签的字体大小
         number.font = 5,             # 相关性数字标签的字体
         insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
         pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
         pch.col = "black",           # 图案颜色
         pch.cex = 10,                 # 图案大小
         tl.cex = 2,                  # 坐标轴标签字体的大小
         cl.cex = 2,                # 图例的字体大小
         p.mat = p_ST,            # P值矩阵
         sig.level = 0.05  )           # 当p大于sig.level时触发动作
#HH-HH
corrplot(auto_HH,
         tl.col="black",
         method="number",
         tl.srt=45,
         addCoef.col="black",
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100), # 主体颜色
         cl.length = NULL,
         number.cex = 2,              # 相关性数字标签的字体大小
         number.font = 2,             # 相关性数字标签的字体
         insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
         pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
         pch.col = "black",           # 图案颜色
         pch.cex = 10,                 # 图案大小
         tl.cex = 2,                  # 坐标轴标签字体的大小
         cl.cex = 2,                # 图例的字体大小
         p.mat = p_HH,            # P值矩阵
         sig.level = 0.05  )           # 当p大于sig.level时触发动作
corrplot(auto_PZ,
         method = "number",           # 图案形状 "square"方框,"circle"圆, "ellipse"椭圆, "number"数字, "shade"阴影花纹, "color"颜色方框, "pie饼图"
         type = "full",               # 绘制范围"full"全部, "lower"下半部分, "upper"半部分
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100), # 主体颜色
         bg = "white",                # 背景颜色
         # col.lim = c(-1,1),         # 数据颜色的范围，是相关性数据的话，直接is.corr = T就好
         title = "",                  # 标题
         is.corr = T,                 # 输入的矩阵是否是相关性矩阵，如果是的话，数据范围会限制到-1到1
         add = F,                     # 是否在原来的图层上添加图形
         diag = T,                    # 是否显示主对角
         outline = F,                 # 图案的轮廓，True或False或某一颜色值
         mar = c(0, 0, 0, 0),         # 下 左 上 右 边距
         addgrid.col = NULL,          # 网格线的颜色，NA为不绘制，NULl为默认的灰色
         addCoef.col = NULL,          # 当method!="number"时，是否显示相关性数值，显示的颜色
         addCoefasPercent = F,        # 是否把相关性数值改为百分数
         order = "original",          # 排序方式 c("original", "AOE", "FPC", "hclust", "alphabet"), original：原始状态，alphabet：字母顺序 hclust，分层聚类顺序
         hclust.method = c("complete", "ward", "ward.D", "ward.D2", "single", "average","mcquitty", "median", "centroid"), # 当order = "hclust"时，分层聚类的算法
         tl.pos = "lt",               # 坐标轴标签的位置'lt', 'ld', 'td', 'd' or 'n'   # 左边 d中间
         tl.cex = 1,                  # 坐标轴标签字体的大小
         tl.col = "black",            # 坐标轴标签字体的颜色
         tl.offset = 0.4,             # 坐标轴标签离图案的距离
         tl.srt = 90,                 # 坐标轴标签旋转角度
         cl.pos = "r",                # 图例位置：r：右边 b：下边 n：不显示
         cl.length = NULL,            # 数字越大，图例的分隔越稠
         cl.cex = 1,                # 图例的字体大小
         cl.ratio = 0.15,             # 图例的宽度
         cl.align.text = "c",         # 图例文字的对齐方式 l左对齐 c居中 r右对齐
         cl.offset = 0.5,             # 图例文字距离图例颜色条的距离 居中时无效
         number.cex = 3,              # 相关性数字标签的字体大小
         number.font = 5,             # 相关性数字标签的字体
         number.digits = 2,           # 相关性数字标签，保留的小数点位数
         na.label = "",               # 当为NA时，显示的内容
         p.mat = p_PZ,           # P值矩阵
         sig.level = 0.05,            # 当p大于sig.level时触发动作
         insig = "pch",               # p值大于sig.level时的方案"pch"图案, "p-value"P值数字, "blank"空白, "n"无操作, "label_sig"星号,
         pch = 4,                     # 当insig = "pch"时的图案形状 4为叉
         pch.col = "black",           # 图案颜色
         pch.cex = 3,                 # 图案大小
         plotCI = "n",                # c("n", "square", "circle", "rect"),  # p值置信区间的方案
         lowCI.mat = testRes$lowCI,   # p值置信区间下边界数据
         uppCI.mat = testRes$uppCI,   # p值置信区间上边界数据
)
corrplot(auto_PZ, 
         add = TRUE, 
         method = "number", 
         col = "black", 
         diag = T, 
         tl.pos = "n", 
         cl.pos = "n")
#HM
corrplot(auto_HM,
         tl.col="black",
         method="pie",
         col=my_color,
         tl.srt=45,
         addCoef.col="black")

#HY
corrplot(auto_HY,
         tl.col="black",
         method="pie",
         col=my_color,
         tl.srt=45,
         addCoef.col="black")

#ST
corrplot(auto_ST,
         tl.col="black",
         method="pie",
         col=my_color,
         tl.srt=45,
         addCoef.col="black")
#hh-fe
corrplot(auto_HH,
         tl.col="black",
         method="pie",
         col=my_color,
         tl.srt=45,
         addCoef.col="black")
####yq-fig2####
#差异基因图
head(sub_astrocyte.markers)
diff_cell <- sub_astrocyte.markers
#显著性
diff_cell$label <- ifelse(diff_cell$p_val_adj<0.01,"adjust P-val<0.01", "adjust P-val>=0.01")
top10_cell = diff_cell %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 
diff_cell$size <- case_when(!(diff_cell$gene %in% top10_cell$gene)~1,
                            diff_cell$gene %in% top10_cell$gene ~2)
dt <- filter(diff_cell,size ==1)

dfbar <- data.frame(x=c("0","1","2","3","4","5","8"),
                    y=c(3.5,3,3,3,3,3,3))
p <- ggplot() + 
  geom_col(data = dfbar,
           mapping = aes(x = x, y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 0.85,
              width = 0.4)+
  geom_jitter(data = top10_cell,
              aes(x = cluster, y = avg_log2FC, color = label),
              size = 1,
              width = 0.4)

#添加显著基因
dfcol <- data.frame(x=c("0","1","2","3","4","5","8"),
                    y=0,
                    label=c("0","1","2","3","4","5","8"))
mycol <- c("#E64B357F","#00A0877F","#3C54887F","#F39B7F7F","#1E90FF","#BC8F8F7F","#00CD00")
library(ggrepel) 
p1 <-  p + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.5,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F) +
  geom_text_repel(
    data=top10_cell,
    aes(x = cluster, y = avg_log2FC, label = gene),
    size = 3,
    arrow = arrow(length = unit(0.008,"npc"),
                  type = "open",ends = "last"))
#对图进行修饰
p2 <- p1+
  scale_color_manual(name = NULL,
                     values = c("red","grey"))+
  labs(x="cluster",y = "avg_log2FC")+
  geom_text(data = dfcol,
            aes(x = x, y = y, label = label),
            size = 3,
            color = "black")+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 12,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 0.5),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 10)
  )

####yq-fig6####
library(ggunchained)
library(reshape2)
library(plyr)
library(reshape2)  
library(ggsci)
library(ggpubr)
library(scales)
library(tidyverse)
library(rstatix)
library(tidyverse)
library(ggpubr)
#差异表达小提琴图
tsw <- c("Tfrc","Slc40a1","Fth1","Ftl1","Aco1","Slc11a2","Ireb2")
VlnPlot(subset(sub_astrocyte, downsample = 100),group.by = "sample",
        features = unique(tsw), ncol = 5,assay = "RNA",cols = c("#6495ED","#FF6A6A"))
#差异表达箱线图
gene=c("Tfrc","Slc40a1","Fth1","Ftl1","Aco1","Slc11a2","Ireb2")
sce.sub.input <- GetAssayData(sub_astrocyte, assay = "RNA", slot = "data") # normalized data matrix
sce.sub.input=sce.sub.input[rownames(sce.sub.input)%in%gene,]
sce.sub.input=as.data.frame(sce.sub.input)
sce.sub.input2 <- data.frame(t(sce.sub.input))
sce.sub.input2$sample = row.names(sce.sub.input2)
data_new = sce.sub.input2
data_new$group = NA
data_new$group[which(str_detect(data_new$sample, "^C"))] <- "Control"
data_new$group[which(str_detect(data_new$sample, "^P"))] <- "Pb"
data_new = melt(data_new)
colnames(data_new) = c("group","sample","gene", "expression")
mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                 axis.title = element_text(size = 12,color ="black"), 
                 axis.text = element_text(size= 12,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.text.x = element_text(angle = 45, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 12),
                 legend.title= element_text(size= 12)
)
p = ggboxplot(data_new,x="gene", y="expression",color = "group",
          palette =c("#E7B800", "#00AFBB"),#分组着色,
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          size=1, #箱型图边线的粗细
          #outlier.shape=NA, #不显示outlier
          legend = "right")
my_comparisons <- list(c("Control", "Pb"))
p + stat_compare_means(aes(group = group),method = "t.test")+mytheme

####yq-fig7####
library(rio)
library(tidyr)
TF <- import("yq_fig/Fig7/Mus_musculus_TF.txt")
TF_symbol <- TF$Symbol
TF_symbol <- TF_symbol[!duplicated(TF_symbol)]
TF <- import("yq_fig/Fig7/AnimalTFDB_tfbs_predict_result.txt")
####yq-fig8####
#a
install.packages("ggupset")
R.Version()
devtools::install_github("YuLab-SMU/DOSE")
devtools::install_github("YuLab-SMU/HDO.db")
devtools::install_github("YuLab-SMU/clusterProfiler")
packageVersion("clusterProfiler")
library(ggupset)
load(file = "bs_code/AS/bs_AS_cx.Rdata")
cell_type_cols <- c(brewer.pal(9, "Set2"), "#FF34B3", "#BC8F8F", "#20B2AA", "#00F5FF", 
                    "#FFA500", "#ADFF2F", "#FF6A6A", "#7FFFD4", "#AB82FF", "#90EE90", 
                    "#00CD00", "#008B8B", "#6495ED", "#FFC1C1", "#CD5C5C", "#8B008B",
                    "#FF3030", "#7CFC00", "#000000", "#708090")
DimPlot(sub_astrocyte, label = T, cols = cell_type_cols,  pt.size = 0.6, repel = T, split.by = "sample")
cell.used <- rownames(sub_astrocyte@meta.data[which(sub_astrocyte@meta.data$seurat_clusters%in% c(5,6,8)),])
sub_astrocyte <- subset(sub_astrocyte, cells = cell.used)
allCells=names(Idents(sub_astrocyte))
allType = levels(Idents(sub_astrocyte))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sub_astrocyte)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,1000,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sub_astrocyte[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(sub_astrocyte))
exprSet <- sub_astrocyte@assays[["RNA"]]@data
exprSet <- as.data.frame(t(exprSet))
y <- as.numeric(exprSet[,"Fth1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)
for(i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type = "pearson")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
View(cor_data_df)
library(ggstatsplot)
ggscatterstats(data = exprSet,
               y = Fth1,
               x = Tubb2b,
               centrality.para = "mean",
               margins = "both",
               xfill = "#CC79A7",
               yfill = "#009E73",
               marginal.type = "densigram",
               title = "relationship")
FeatureScatter(cg_sce,feature1 = "Tfrc", feature2 = "Tubb2b")
FeatureScatter(sub_astrocyte,feature1 = "Tfrc", feature2 = "Hif1a")
#b
diffgene <- import(file = "bs_table/AS/sub_astrocyte.diffgene.csv")
diffgene <- diffgene[,c(1:3)]
names(diffgene) <- c("symbol","pvalue","logfc")
all_marker = diffgene$symbol
all_marker <- bitr(all_marker,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)
diffgene <- merge(all_marker,diffgene,by = 1,all = F)
geneList = diffgene[,4]
names(geneList) = as.character(diffgene[,2])
head(geneList)
geneList <- geneList[!duplicated(names(geneList))]
geneList = sort(geneList, decreasing = TRUE)
diffgene_gseago_BP <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                                            OrgDb = org.Mm.eg.db,
                                            ont = "BP",#可选择bp.MF,CC,ALL
                                            nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                                            minGSSize = 1,#最小基因集的基因数
                                            maxGSSize = 500,#最大基因集的基因数
                                            pvalueCutoff = 10,
                                            pAdjustMethod = "BH",#p值的阈???
                                            verbose      = FALSE)#是否输出提示信息，默认为false
write.table(diffgene_gseago_BP,file="diffgene_gseago_BP.csv",sep=",", quote=F,row.names = F)
diffgene_gseago_CC <- gseGO(geneList = geneList,#排序后的基因列表,100基因以上的要求
                         OrgDb = org.Mm.eg.db,
                         ont = "CC",#可选择bp.MF,CC,ALL
                         nPerm = 1000,#置换检验的次数，默???1000，保持默认即???
                         minGSSize = 1,#最小基因集的基因数
                         maxGSSize = 500,#最大基因集的基因数
                         pvalueCutoff = 10,
                         pAdjustMethod = "BH",#p值的阈???
                         verbose      = TRUE)#是否输出提示信息，默认为false
save(geneList, diffgene_gseago_BP, diffgene_gseago_CC, file = "gsea_result.Rda")
load(file = "gsea_result.Rda")
write.table(diffgene_gseago_CC,file="diffgene_gseago_CC.csv",sep=",", quote=F,row.names = F)
diffgene_GSEAKEGG <-  gseKEGG(geneList = geneList, organism = 'mmu',nPerm = 1000,minGSSize = 3, pvalueCutoff = 1, verbose  = FALSE)

gseaplot2(diffgene_gseago_BP, 10,pvalue_table = TRUE, color = "#0AFF99")
####yq-fig9####
library(enrichplot)
load(file = "gsea_result.Rda")
mito_gene <- import("yq_fig/Fig8/mito_gene.xlsx")
diff_gene <- import("yq_fig/Fig8/sub_astrocyte.diffgene.csv")
diff_gene <- diff_gene$V1
a <- intersect(diff_gene$V1,mito_gene$Symbol)
diff_gene <- merge(a,diff_gene,by = 1,all = F)
write.csv(diff_gene,file = "yq_fig/Fig8/diff_mitogene.csv")
gseaplot2(diffgene_gseago_CC, 24,pvalue_table = TRUE, color = "#0AFF99")#线粒体mitochondrion
gseaplot2(diffgene_gseago_CC, 20,pvalue_table = TRUE, color = "#0AFF99")#lysosome
gseaplot2(diffgene_gseago_CC, 38,pvalue_table = TRUE, color = "#0AFF99")#endoplasmic 
gseaplot2(diffgene_gseago_CC, 330,pvalue_table = TRUE, color = "#0AFF99")#reticulumGolgi apparatus
gseaplot2(diffgene_gseago_CC, geneSetID = c(24,20,38,330), subplots = 1:3,
          color = c("#E495A5", "#86B875", "#7DB0DD","#DEAAFF"))
gseaplot2(diffgene_gseago_CC, geneSetID = c(24,20,38,330), subplots = 1:3,
          color = c("#E495A5", "#86B875", "#7DB0DD","#DEAAFF"))
gseaplot2(diffgene_gseago_BP, geneSetID = c(5,21,412), subplots = 1:3,
          color = c("#E495A5", "#86B875", "#7DB0DD"))

####rq####
install.packages("doBy")
library(doBy)
library(patchwork)
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
data <- import("rqsj.xlsx")
colnames(data)<- c("bh","nl","xl","xq","rz","moc","xy","yj","jq")

data <- data[,-1]
moresatas <- function(x){c(mean=mean(x),sd=sd(x))}
summaryBy(nl+xl+xq+rz+moc~jq,data = data,FUN = moresatas)
##control
data_control <- data[which(data$jq%in%c(0)),]
table(data_control$xl)
##Pb
data_Pb <- data[which(data$jq%in%c(1)),]
table(data_Pb$xl)

####rdfx####
load("bs_code/bs_totalcell.Rdata")
cell.used <- rownames(sce.all@meta.data[-which(sce.all$celltype %in% c('Mic_Ast',"Neural stem cell","Ependymal cell", "Neutrophils")),])
cell.used <- rownames(sce.all@meta.data[-which(sce.all$celltype %in% c('Mic_Ast')),])
table(Idents(sce.all))
length(cell.used)
sce.all <- subset(sce.all, cells = cell.used)
DimPlot(sce.all, label = F, pt.size = 0.5)

data("sc_sim")
head(sc_sim@meta.data)
augur = calculate_auc(sc_sim)
head(augur$AUC, 11)
head(sce.all@meta.data)
cg_sce <- subset(cg_sce, idents=c(0,1,5))

allCells=names(Idents(sce.all))
allType = levels(Idents(sce.all))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,5000,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sce.all[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(cg_sce))
DimPlot(cg_sce, label = T, pt.size = 0.5)
meta.data1 <- data.frame(lable = NA, cell_type = NA)

augur = calculate_auc(cg_sce, cg_sce@meta.data, cell_type_col = "celltype", label_col = "sample",n_threads = 8)
a <- augur$feature_importance
b <- augur$AUC
write.csv(b,file="b.csv")

augur = calculate_auc(sce.all, sce.all@meta.data, cell_type_col = "celltype", label_col = "sample",n_threads = 8)
a <- augur$feature_importance
b <- augur$AUC
write.csv(b,file="b1.csv")
astro_mimportance <- a[a$cell_type %in% c("Astrocyte"),]
write.csv(astro_mimportance,file="astro_mimportance.csv")
Mic_mimportance <- a[a$cell_type %in% c("Microglia"),]
write.csv(Mic_mimportance,file="Mic_mimportance.csv")

####衰老打分####
library(patchwork)
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(Augur)
library(ggplot2) 
library(cowplot) 
library(paletteer)  
library(gplots)
library(ggpubr)    
library(ggsci) 
library(stringr)
library(rio)
load("bs_code/bs_totalcell.Rdata")
HH <- import("CellAge_Senescence Genes.csv")
GENE <- HH$`Gene Symbol`

capitalize_first <- function(x) {
  # 将第一个字母大写，其余字母小写
  paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
}

LIST <- capitalize_first(GENE) 
OLD <- list(c(LIST))

allCells=names(Idents(sce.all))
allType = levels(Idents(sce.all))
choose_Cells = unlist(lapply(allType, function(x){
  cgCells = allCells[Idents(sce.all)== x ] #这里是将属于每个亚群的细胞取出来
  cg=sample(cgCells,1000,replace = T) #从取出的亚群中进行随机抽样，这里使用的是sample函数
  cg #返回随机取到的样本，储存在choose_Cells中
}))
cg_sce = sce.all[, allCells %in% choose_Cells]
class(cg_sce)
table(Idents(cg_sce))
DimPlot(cg_sce, label = T, pt.size = 0.5)


ferroptosis_all <- list(c("Fth1","Slc40a1","Tfrc","Slc11a2","Slc25a28","Slc25a37","Vdac2","Vdac3","Slc39a14",
                          "Steap3","Ncoa4","Nfs1","Cisd1","Cisd2","Mfsd7b","Abcb7","Abcb8","Fdx1","Aco1","Ireb2",
                          "Slc7a11","Gpx4","Hmox1","Nfe2l2","Dhodh","Slc3a2","Aifm2","Gclc","Keap1","Abcc1","Atf3","Atf4",
                          "Cdo1","Me1","Coq2","Gch1","Slc25a39","Slc3a2","Sqstm1","Ncoa4","Hif1a",
                          'Acsl4',"Trp53","Ptgs2","Pebp1","Acsl1","Acsl3","Cs","Hmgcr","Sqle","Alox5","Alox15","Alox12",
                          "Phgdh","Hmgb1","Pparg"))
WNT_features <- OLD
cg_sce <- AddModuleScore(cg_sce,features = WNT_features,ctrl = 100,name = "WNT_features")
p <- ggboxplot(cg_sce@meta.data, x="celltype", y="WNT_features1", width = 0.6, 
               color = "sample",#轮廓颜色
               palette =c("#1E90FF", "#FF6347"),#分组着色
               xlab = F, #不显示x轴的标签
               x.text.angle = 90,
               bxp.errorbar=T,#显示误差条
               bxp.errorbar.width=0.5, #误差条大小
               size=0.5, #箱型图边线的粗细
               outlier.shape=NA, #不显示outlier
               legend = "right")
compare_means(WNT_features1 ~ sample, data = cg_sce@meta.data,group.by = "celltype")
p + stat_compare_means(aes(sample = sample),label = "p.signif")