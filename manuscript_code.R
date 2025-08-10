####加载包####
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
library(ggplot2) 
library(cowplot) 
library(paletteer)  
library(gplots)
library(ggpubr)    
library(ggsci) 
library(stringr)
library(rio)
library(irGSEA)
library(doBy)
library(patchwork)
library(dplyr)
library(stringr)
library(ggplot2)
library(Seurat)
library(augur)
####Figure2####

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
#铁死亡打分
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

####Figure3####
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
load(file = "bs_code/AS/bs_AS_cx.Rdata")
diff_as3 <- import("bs_table/AS/sub_astrocyte.diffgene3.csv")
diff_as5 <- import("bs_table/AS/sub_astrocyte.diff5.csv")
diff_as8 <- import("bs_table/AS/sub_astrocyte.diffgene8.csv")
a <- intersect(diff_as3$V1,intersect(diff_as5$gene,diff_as8$gene))
#拟时序分析
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
####Figure5####
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