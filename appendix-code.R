library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)
library(tidyverse)
library(reshape2)
library(scales)
library(gdata)
library(ggraph)

library(scater)
library(cowplot)
library(SingleR)
library(data.table)
library(paletteer)
library(R.utils)
library(ComplexHeatmap)
library(clustree)
library(viridis)
library(pheatmap)
library(devtools)
library(Scillus)
library(ggcharts)
library(Cairo)
library(cartography)
library(knitr)
library(patchwork)
library(BiocParallel)
library(nichenetr)
library(purrr)
library(scFunctions)
library(kableExtra)
library(gsfisher)
library(easyalluvial)
library(MySeuratWrappers)

####HUMAN####
####input####
GM136 <- CreateSeuratObject(Read10X("./original date/GM136/"),project = 'GM136')
GM143 <- CreateSeuratObject(Read10X("./original date/GM143/"),project = 'GM143')
GM144 <- CreateSeuratObject(Read10X("./original date/GM144/"),project = 'GM144')
GM147 <- CreateSeuratObject(Read10X("./original date/GM147/"),project = 'GM147')
GM148 <- CreateSeuratObject(Read10X("./original date/GM148/"),project = 'GM148')
GM169 <- CreateSeuratObject(Read10X("./original date/GM169/"),project = 'GM169')
GM183 <- CreateSeuratObject(Read10X("./original date/GM183/"),project = 'GM183')
GM184a <- CreateSeuratObject(Read10X("./original date/GM184a/"),project = 'GM184a')
GM238 <- CreateSeuratObject(Read10X("./original date/GM238/"),project = 'GM238')
GM241 <- CreateSeuratObject(Read10X("./original date/GM241/"),project = 'GM241')
GM242 <- CreateSeuratObject(Read10X("./original date/GM242/"),project = 'GM242')
GM283 <- CreateSeuratObject(Read10X("./original date/GM283/"),project = 'GM283')
GM289 <- CreateSeuratObject(Read10X("./original date/GM289/"),project = 'GM289')

PD134 <- CreateSeuratObject(Read10X("./original date/PD134/"),project = 'PD134')
PD153 <- CreateSeuratObject(Read10X("./original date/PD153/"),project = 'PD153')
PD161 <- CreateSeuratObject(Read10X("./original date/PD161/"),project = 'PD161')
PD161b <- CreateSeuratObject(Read10X("./original date/PD161b/"),project = 'PD161b')
PD164 <- CreateSeuratObject(Read10X("./original date/PD164/"),project = 'PD164')
PD164b <- CreateSeuratObject(Read10X("./original date/PD164b/"),project = 'PD164b')
PD164c <- CreateSeuratObject(Read10X("./original date/PD164c/"),project = 'PD164c')
PD170 <- CreateSeuratObject(Read10X("./original date/PD170/"),project = 'PD170')

GM <- merge(GM136, y = c(GM143, GM144, GM147, GM148, GM183, GM184a, GM238, GM241, GM242, GM169, GM283, GM289), 
            add.cell.ids = c("GM136", "GM143", "GM144", "GM147", "GM148", "GM183", "GM184a", 
                             "GM238", "GM241", "GM242", "GM169", "GM283", "GM289"), 
            project = "GM")
PD=merge(PD134,y=c(PD153,PD161,PD161b,PD164,PD164b,PD164c,PD170),
         add.cell.ids = c('PD134','PD153','PD161','PD161b','PD164','PD164b','PD164c',
                          'PD170'),
         project = "PD")

GM[["percent.mt"]] <- PercentageFeatureSet(GM, 
                                           pattern = "^MT-")
PD[["percent.mt"]] <- PercentageFeatureSet(PD, 
                                           pattern = "^MT-")

GM <- subset(GM, 
             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
PD <- subset(PD, 
             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)GM <- NormalizeData(GM)

GM <- NormalizeData(GM)
PD <- NormalizeData(PD)

GM <- FindVariableFeatures(GM, nfeatures = 4000)
PD <- FindVariableFeatures(PD, nfeatures = 4000)

GM@meta.data$group='H'
PD@meta.data$group='P'

####cluster####
sampleList <- list(GM, PD)
options(future.globals.maxSize= 5000*1024^2)
oral.anchors <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
oralIntegrated <- IntegrateData(anchorset = oral.anchors, dims = 1:50)

#Dimensionality reduction and cell clustering
# cell-cycle scoring and regression
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
oralIntegrated <- CellCycleScoring(oralIntegrated, 
                                   s.features = s.genes, 
                                   g2m.features = g2m.genes,
                                   set.ident = TRUE)
# standard workflow for clustering
oralIntegrated <- ScaleData(oralIntegrated, 
                            vars.to.regress = c("S.Score", "G2M.Score"), 
                            verbose = FALSE)
oralIntegrated <- RunPCA(oralIntegrated, 
                         npcs = 50, 
                         verbose = FALSE)
oralIntegrated <- FindNeighbors(oralIntegrated, 
                                reduction = "pca", 
                                dims = 1:50)
oralIntegrated <- FindClusters(oralIntegrated, 
                               resolution = seq(from = 0.1, 
                                                to = 1.0, 
                                                by = 0.1))
oralIntegrated <- RunUMAP(oralIntegrated, 
                          reduction = "pca", 
                          dims = 1:50)


oralIntegrated_HvP=oralIntegrated
DefaultAssay(oralIntegrated_HvP) <- "RNA"
DefaultAssay(oralIntegrated_HvP) <- "integrated"
Idents(oralIntegrated_HvP) <- "integrated_snn_res.1"
all.markers_HvP <- FindAllMarkers(oralIntegrated_HvP, 
                                  only.pos = TRUE)

significant.markers_HvP <- all.markers_HvP[all.markers_HvP$p_val_adj < 0.2, ]

top20.markers_HvP <- significant.markers_HvP %>% 
  group_by(cluster) %>% 
  top_n(n=20, wt=avg_log2FC)
write.csv(top20.markers_HvP, "integrated_HvP_top20Markers.csv")
write.table(significant.markers_HvP, "./Export/HvP marker.txt", sep="\t",quote = F,row.names = F)

#rename#
oralIntegrated_HvP <- subset(oralIntegrated_HvP,idents = c('23'), invert = TRUE)
oralIntegrated_HvP <- RenameIdents(oralIntegrated_HvP, 
                                   "0"="Endothelial", 
                                   "3"="Endothelial", 
                                   "8"="Endothelial", 
                                   "11"="Endothelial", 
                                   "32"="Endothelial", 
                                   "21"="Endothelial", 
                                   "4"="Endothelial",
                                   "28"="Endothelial",
                                   "12"="Endothelial", 
                                   
                                   "1"="Fibroblast", 
                                   "13"="Fibroblast", 
                                   "14"="Fibroblast",
                                   "16"="Fibroblast",
                                   "20"="Fibroblast", 
                                   "22"="Fibroblast",
                                   
                                   
                                   "17"="Immune", 
                                   "12"="Immune",
                                   "19"="Immune",
                                   "29"="Immune",
                                   "5"="Immune",
                                   "10"="Immune",
                                   "15"="Immune",
                                   "2"="Immune",
                                   "6"="Immune",
                                   "9"="Immune",
                                   "18"="Immune",
                                   "26"="Immune",
                                   "31"="Immune",
                                   
                                   "7"="Epithelial", 
                                   "24"="Epithelial", 
                                   "27"="Epithelial", 
                                   "30"="Epithelial", 
                                   
                                   "25"="Other")
levels(oralIntegrated_HvP)=c("Endothelial","Fibroblast","Immune","Epithelial","Other")
oralIntegrated_HvP$generalCellTypes <- oralIntegrated_HvP@active.ident

#figure 1A
pdf(file ='./Export-NIKI/Dim.pdf',width = 8,height = 8 )
oralIntegrated_HvP2=subset(oralIntegrated_HvP,idents=c("Endothelial","Fibroblast","Immune","Epithelial"))
Seurat::DimPlot(object = oralIntegrated_HvP2, 
                label = F, pt.size = 1.5,
                reduction = "umap", raster = T,raster.dpi = c(1200, 1200),
                cols = pal) +
  theme(legend.text = element_text(size = 6),legend.position='right')
dev.off()

#figure S1A
# selected genes to display
GenesEndo <- c("ACKR1", 
               "RAMP2", 
               "SELE", 
               "VWF", 
               "PECAM1")
GenesFib <- c("LUM", 
              "COL3A1", 
              "DCN", 
              "COL1A1", 
              "CFD")
GenesIm <- c("CD69", 
             "CD52", 
             "CXCR4", 
             "PTPRC", 
             "HCST")
GenesEpi <- c("KRT14", 
              "KRT5", 
              "S100A2", 
              "CSTA", 
              "SPRR1B")
clustergeen=c(GenesEpi,GenesFib,GenesEndo,GenesIm)

pdf(file = './Export/clustergeen.pdf' ,width =4,height = 4)
DotPlot(subset(oralIntegrated_HvP,idents = c("Endothelial",'Epithelial','Fibroblast','Immune')),
        features = clustergeen,dot.scale = 4)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
dev.off()

#figure 1B
IL36s=c('IL36G','IL36B','IL36A','IL1F10','IL36RN')
pdf(file = './Export/IL36s.pdf' ,width = 8,height = 4)
DotPlot(subset(oralIntegrated_HvP,idents = c("Endothelial",'Epithelial','Fibroblast','Immune')),
        features = IL36s,dot.scale = 15)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
dev.off()

####epithelial####

epi_HvP <- subset(oralIntegrated_HvP, 
                  idents = c("Epithelial"))
Idents(epi_HvP) <- "integrated_snn_res.1"
epi2 <- epi_HvP
DefaultAssay(epi2) <- "integrated"
epi2 <- FindNeighbors(epi2, 
                      dims = 1:50)
epi2 <- FindClusters(epi2, 
                     resolution = seq(from = 0.1, 
                                      to = 1.0, 
                                      by = 0.1))
Idents(epi2) <- "integrated_snn_res.0.1"
epi2Markers <- FindAllMarkers(epi2, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.75)
epi2Sig <- epi2Markers[epi2Markers$p_val_adj < 0.2, ]
epi2Top <- epi2Sig %>%
  group_by(cluster) %>%
  top_n(n=5, wt=avg_log2FC)
write.table(epi2Sig, "./Export/epi_HvPsig.txt", sep="\t",quote = F,row.names = F)
print(epi2Top, n=40)
epi_HvP <- epi2
epi2 <- NULL
DimPlot(epi_HvP)
epi_HvP <- RunUMAP(epi_HvP, dims = 1:50)
#renmae#
GenesBasal=c('COL17A1','DST','KRT15')
GenesJunc= c('FDCSP','ODAM','SLPI','IL36G')
GenesGranu=c('CRNN','SPRR3','CNFN')
GenesSpin=c('KRT1','KRT16','KRT5') 
GenesMel=c('PMEL','MLANA','TYRP1')
clusterepi=c(GenesJunc,GenesGranu,GenesSpin,GenesBasal,GenesMel)

epi_HvP <- RenameIdents(epi_HvP, '0'='Spinosum', 
                        '1'='Basale',
                        "2"="Junctional", 
                        '3'='Granulosum',
                        "4"="Melanocytes")
#figure 1D
pdf(file = './Export/Epi-IL36s.pdf' ,width = 4,height = 4)
Seurat::VlnPlot(subset(epi_HvP,idents=c("Junctional", "Granulosum", 'Spinosum','Basale')),
                'IL36G',split.by = 'group',pt.size = 0,adjust = 1,cols = pal[c(2,1)])
dev.off()
#figure S1B
pdf(file = './Export/clustergeen-2.pdf' ,width =4,height = 4)
DotPlot(epi_HvP,
        features = clusterepi,dot.scale = 4)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
dev.off()
####cell chat####
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
data.input=oralIntegrated_HvP[["RNA"]]@counts 
meta = oralIntegrated_HvP@meta.data 

meta=meta[,-c(2:3)]
meta$condition=meta$group
meta$labels=meta$generalCellTypes
colnames(meta)[1]='patient.id'
table(meta$condition) 
meta2222=meta
cell.use = rownames(meta)[meta$condition == c("P")] 

data.input = data.input[, cell.use]
meta = meta[cell.use, ]

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "truncatedMean",trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F)

cellchat@netP$pathways

pathways.show <- c("CXCL") 
# [1] "CXCL"      "MIF"       "CCL"       "VISFATIN"  "MK"        "ncWNT"     "PERIOSTIN" "IGF"      
#[9] "PTN"       "ANGPTL"    "FGF"       "CALCR"     "VEGF"      "TGFb"      "GALECTIN"  "PDGF"     
#[17] "GAS"       "IL6"       "EGF"       "TWEAK"     "PARs"      "ANGPT"     "WNT"       "PROS"     
#[25] "BMP"       "TRAIL"     "SEMA3"     "EDN"       "IL16"      "CD70"      "NRG"       "BAFF"     
#[33] "HGF"  

#figure 4A#
pdf(file = './Export-NIKI/cellchat.pdf' ,width = 6,height = 4)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 

#figure 4B#
pdf(file = './Export-NIKI/cellchat2.pdf' ,width = 6,height = 4)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,cluster.cols = pal,
                                  width = 8, height = 2.5, font.size = 10)
dev.off()

#figure 4C#
pdf(file = './Export-NIKI/IL36R.pdf' ,width = 8,height = 4)
DotPlot(oralIntegrated_HvP, features = 'IL1RL2',dot.scale = 15,cols = c('green',"red"))+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))
dev.off()

####MOUSE####
####input####
C <- CreateSeuratObject(Read10X("./original mice/filtered_feature_bc_matrix-c/"),project = 'C')
P <- CreateSeuratObject(Read10X("./original mice/filtered_feature_bc_matrix-12/"),project = 'P')
GM=C
BM=P

GM[["percent.mt"]] <- PercentageFeatureSet(GM, 
                                           pattern = "^mt-")
BM[["percent.mt"]] <- PercentageFeatureSet(BM, 
                                           pattern = "^mt-")

BM@meta.data$group='P'
GM@meta.data$group='C'

GM <- subset(GM, 
             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
BM <- subset(BM, 
             subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

BM <- NormalizeData(BM)
BM <- FindVariableFeatures(BM, nfeatures = 4000)
GM <- NormalizeData(GM)
GM <- FindVariableFeatures(GM, nfeatures = 4000)
sampleList <- list(GM, BM)
options(future.globals.maxSize= 5000*1024^2)
oral.anchors <- FindIntegrationAnchors(object.list = sampleList, dims = 1:50)
oralIntegrated <- IntegrateData(anchorset = oral.anchors, dims = 1:50)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes=tolower(s.genes)
g2m.genes=tolower(g2m.genes)

library(stringr)
s.genes=str_to_title(s.genes)
g2m.genes=str_to_title(g2m.genes)


oralIntegrated <- CellCycleScoring(oralIntegrated, 
                                   s.features = s.genes, 
                                   g2m.features = g2m.genes,
                                   set.ident = TRUE)
# standard workflow for clustering
oralIntegrated <- ScaleData(oralIntegrated, 
                            vars.to.regress = c("S.Score", "G2M.Score"), 
                            verbose = FALSE)

oralIntegrated <- RunPCA(oralIntegrated, 
                         npcs = 50, 
                         verbose = FALSE)
oralIntegrated <- FindNeighbors(oralIntegrated, 
                                reduction = "pca", 
                                dims = 1:50)
oralIntegrated <- FindClusters(oralIntegrated, 
                               resolution = seq(from = 0.1, 
                                                to = 1.0, 
                                                by = 0.1))
oralIntegrated <- RunUMAP(oralIntegrated, 
                          reduction = "pca", 
                          dims = 1:50)
# assess cluster tree
clustree(oralIntegrated)

DefaultAssay(oralIntegrated) <- "RNA"
Idents(oralIntegrated) <- "integrated_snn_res.1"
all.markers.1 <- FindAllMarkers(oralIntegrated, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.75)
significant.markers.1 <- all.markers.1[all.markers.1$p_val_adj < 0.2, ]
write.table(significant.markers.1, "./Export/mar.txt", sep="\t",quote = F,row.names = F)
MICE_oralintegrated=oralIntegrated

MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('0','5','8','3','12','16','17')] = 'Epithelial'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('20','28','18','39','28','35')] = 'Epithelial'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('2','9','13','14','21')] = 'Fibroblast'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('7','36')] = 'Endothelial'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('4','23','27')] = 'T cell'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('6','22')] = 'B cell'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('1','11','26','29','34','15')] = 'Neutrophil'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('10','19','24','32','37')] = 'Mono/Mac'
MICE_oralintegrated$Label[MICE_oralintegrated$integrated_snn_res.1 %in% c('33','25','30','31','38')] = 'Other'

Idents(MICE_oralintegrated)='Label'
levels(MICE_oralintegrated)=c("Epithelial","Fibroblast","Endothelial","Neutrophil",
                              "Mono/Mac", "T cell" ,"B cell" ,"Other" )
pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,10,4,5,9,6,7)]


#figure S2A
GenesEndo <- c("Egfl7", "Ramp2","Pecam1",'Cdh5','Cd93')
GenesFib <- c("Lum","Col3a1","Dcn","Col1a1",'Col1a2')

GenesEpi <- c("Krt14", "Krt5","Epcam",'Krt10','Krt16')

GenesTcell <- c("Cd3d","Cd3e","Il7r","Cd3g", "Trbc2")
GenesBcell <- c("Cd79a", "Cd37", "Cd19", "Cd79b","Ms4a1")
GenesNeutro=c('Csf3r','Ly6g','S100a8','S100a9','Il1b')
Genesmono=c("Cd68", "C1qa", "C1qb", "Ctss",'Cfp')
micemar=c(GenesEpi,GenesFib,GenesEndo,GenesNeutro,Genesmono,GenesTcell,GenesBcell)

pdf(file = './Export mice0613/clustergeen.pdf' ,width =4,height = 6)
DotPlot(MICE_oralintegrated,
        features = micemar,dot.scale = 4)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
dev.off()

#figure 2A
pdf(file ='./Export/Dim.pdf',width = 8,height = 8 )
MICE_oralintegrated2=subset(MICE_oralintegrated,idents=c("Epithelial","Fibroblast","Endothelial","Neutrophil",
                                                         "Mono/Mac", "T cell" ,"B cell"))

Seurat::DimPlot(object = MICE_oralintegrated2, 
                label = F, pt.size = 2,
                reduction = "umap", raster = T,raster.dpi = c(1200, 1200),
                cols = pal) +
  theme(legend.text = element_text(size = 6),legend.position='right')
dev.off()

#figure 2B
IL36s=c('Il1f9','Il1f8','Il1f6','Il1f10','Il1f5')
pdf(file = './Export/IL36s.pdf' ,width = 5,height = 6)
DotPlot(MICE_oralintegrated2,
        features = IL36s,dot.scale = 15)+coord_flip()+
  theme_bw()+
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('grey','#336699','#66CC66','#FFCC33'))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))
dev.off()

#fugure 2C
pdf(file = './Export/IL36g.pdf' ,width = 5,height = 6)
Seurat::VlnPlot(MICE_oralintegrated2,'Il1f9',split.by = 'group',pt.size = 0,
                adjust = 2)
dev.off()

####cellchat####
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

data.input=MICE_oralintegrated[["RNA"]]@counts 
meta = MICE_oralintegrated@meta.data 

meta$condition=meta$group
meta$labels=meta$Label
colnames(meta)[1]='patient.id'
table(meta$condition) 

cell.use = rownames(meta)[meta$condition == c("P")] 

data.input = data.input[, cell.use]
meta = meta[cell.use, ]

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents)=c( "Epithelial", "Fibroblast","Endothelial",
                           "Mono/Mac","Neutrophil","B cell","T cell","Other")


groupSize <- as.numeric(table(cellchat@idents)) 

CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "truncatedMean",trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "./Export mice0613/cellchat.rds")


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F)


cellchat@netP$pathways

pathways.show <- c("CXCL") 

#figure 5A
pdf(file = './Export mice0613/cellchat-contribution.pdf' ,width = 6,height = 4)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)

#figure 5B
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
pdf(file = './Export mice0613/cellchat-sender.pdf' ,width = 6,height = 4)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show,cluster.cols = pal,
                                  width = 8, height = 2.5, font.size = 10)
dev.off()
