##integrate data with seurat v3
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

samples = read.delim2("Sample_annotation_file.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")

nCoV.list = list()
x=0
for(sample_s in samples$sample){
  print(sample_s)
  sample_i = samples %>% dplyr::filter(.,sample == sample_s)
  #sample_i=sample_s
  datadir = paste(sample_s,"_filtered_feature_bc_matrix.h5",sep="")
  sample.tmp = Read10X_h5(datadir,use.names = TRUE, unique.features = TRUE)
  sample.tmp.seurat <- CreateSeuratObject(counts = sample.tmp, min.cells = 3, min.features = 200,project = sample_s)
  sample.tmp.seurat[['percent.mito']] <- PercentageFeatureSet(sample.tmp.seurat, pattern = "^MT-")
  sample_i$nFeature_RNA_low = as.numeric(sample_i$nFeature_RNA_low)
  sample_i$nFeature_RNA_high = as.numeric(sample_i$nFeature_RNA_high)
  sample_i$nCount_RNA = as.numeric(sample_i$nCount_RNA)
  sample_i$percent.mito = as.numeric(sample_i$percent.mito)
  sample.tmp.seurat <- subset(x = sample.tmp.seurat, subset = nFeature_RNA > sample_i$nFeature_RNA_low & nFeature_RNA < sample_i$nFeature_RNA_high 
                              & nCount_RNA > sample_i$nCount_RNA & percent.mito < sample_i$percent.mito)
  sample.tmp.seurat <- NormalizeData(sample.tmp.seurat, verbose = FALSE)
  sample.tmp.seurat <- FindVariableFeatures(sample.tmp.seurat, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
  nCoV.list[sample_s] = sample.tmp.seurat
}
nCoV <- FindIntegrationAnchors(object.list = nCoV.list, dims = 1:50)
nCoV.integrated <- IntegrateData(anchorset = nCoV, dims = 1:50,features.to.integrate = rownames(nCoV))

####add  sample info
sample_info = as.data.frame(colnames(nCoV.integrated))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = nCoV.integrated@meta.data$orig.ident
sample_info = dplyr::left_join(sample_info,samples)
rownames(sample_info) = sample_info$ID
nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)

###first generate data and scale data in RNA assay
DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))

##change to integrated assay
DefaultAssay(nCoV.integrated) <- "integrated"
dpi = 300
png(file="qc.png", width = dpi*16, height = dpi*8, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

png(file="umi-gene.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
dev.off()

# Run the standard workflow for visualization and clustering
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
png(file="pca.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = nCoV.integrated,ndims = 100)
dev.off()

###cluster
nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2) 

###tsne and umap
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- RunUMAP(nCoV.integrated, reduction = "pca", dims = 1:50)
png(file="tsne.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, label.size = 8)
dev.off()

DefaultAssay(nCoV.integrated) <- "RNA"
# find markers for every cluster compared to all remaining cells, report only the positive ones
nCoV.integrated@misc$markers <- FindAllMarkers(object = nCoV.integrated, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')
write.table(nCoV.integrated@misc$markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')

dpi = 300
png(file="feature.png", width = dpi*24, height = dpi*5, units = "px",res = dpi,type='cairo')
VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"))
dev.off()
saveRDS(nCoV.integrated, file = "nCoV.rds")

hc.markers = read.delim2("marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
tt1 = DoHeatmap(object = subset(nCoV.integrated, downsample = 500), features = top10$gene) + NoLegend()
ggplot2::ggsave(file="marker_heatmap_MAST.pdf",plot = tt1,device = 'pdf',width = 20, height = 16, units = "in",dpi = dpi,limitsize = FALSE)

#draw heatmap
markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF','DCN',
            'FCGR3A','TREM2','KRT18')
hc.markers = read.delim2("marker_MAST.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
hc.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) -> top30
var.genes = c(nCoV.integrated@assays$RNA@var.features,top30$gene,markers)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"),features = var.genes)
saveRDS(nCoV.integrated, file = "nCoV.rds")

dpi = 300
png(file="tsne_m.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap_m.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE)
dev.off()

dpi = 300
png(file="nCoV-umap-group-sample.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'sample_new')
dev.off()

png(file="nCoV-umap-split-sample.png", width = dpi*16, height = dpi*16, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'sample_new', ncol = 4)
dev.off()

png(file="nCoV-umap-group-group.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'group')
dev.off()

png(file="nCoV-umap-split-group.png", width = dpi*12, height = dpi*4, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'group', ncol = 3)
dev.off()

png(file="nCoV-umap-group-disease.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = FALSE, group.by = 'disease')
dev.off()

png(file="nCoV-umap-split-disease.png", width = dpi*10, height = dpi*4.5, units = "px",res = dpi,type='cairo')
DimPlot(object = nCoV.integrated, reduction = 'umap',label = TRUE, split.by = 'disease', ncol = 2)
dev.off()

####marker expression
dpi = 300
markers = c('AGER','SFTPC','SCGB3A2','TPPP3','KRT5',
            'CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
            'CLEC9A','IL3RA',
            'CD3D','CD8A','KLRF1',
            'CD79A','IGHG4','MS4A1',
            'VWF','DCN',
            'FCGR3A','TREM2','KRT18','HBB')
#markers = c('HBB')
markers = c('CLEC9A', 'CD1C', 'FCER1A', 'LILRA4', 'CD68', 'CD14')

markers= c('TPPP3', 'KRT18','IGHG4','CD3D')
png(file="marker/violin_marker.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
print(VlnPlot(object = nCoV.integrated, features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
dev.off()
png(file="marker/umap_marker.png", width = dpi*30, height = dpi*36, units = "px",res = dpi,type='cairo')
print(FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000")))
dev.off()
for(marker in markers){
  png(file=paste("marker/violin_",marker,".png",sep=''), width = dpi*8, height = dpi*3, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = nCoV.integrated, features = marker,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
  dev.off()
  png(file=paste("marker/umap_",marker,".png",sep=''), width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
  print(FeaturePlot(object = nCoV.integrated, features = marker,cols = c("lightgrey","#ff0000")))
  dev.off()
}

library(ggplot2)
pdf(file="marker_heatmap.pdf", width = 10, height = 8)
pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12)) + labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6))
print(pp)
dev.off()

#####################Sub clusters MNP
SelectedClusters <- subset(x = nCoV.integrated, idents = c(22,25,27,28,0,1,2,3,4,5,6,9,10,11,12,13,21,23,29))
 
# Run the standard workflow for visualization and clustering
SelectedClusters <- ScaleData(SelectedClusters, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
SelectedClusters <- RunPCA(SelectedClusters, verbose = FALSE,npcs = 100)
SelectedClusters <- ProjectDim(object = SelectedClusters)
png(file="pca_Selected.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = SelectedClusters,ndims = 100)
dev.off()

###cluster
SelectedClusters <- FindNeighbors(object = SelectedClusters, dims = 1:50)
SelectedClusters <- FindClusters(object = SelectedClusters, resolution = 1.2) 

###tsne and umap
SelectedClusters <- RunTSNE(object = SelectedClusters, dims = 1:50)
SelectedClusters <- RunUMAP(SelectedClusters, reduction = "pca", dims = 1:50)
png(file="tsne_selected.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = SelectedClusters, reduction = 'tsne',label = TRUE)
dev.off()
png(file="umap_selected.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = SelectedClusters, reduction = 'umap',label = TRUE)
dev.off()

###############################

#markers = c('CLEC9A', 'CD1C', 'FCER1A', 'CD68', 'CD14')
markers = c('CD123','AXL','BSG','CCR2','CD38','CD16','CD86','CD200R','CD15','CD44','FceR1','CD163','CD141','CD34','CD3','CD7','CD19','HLA-DR','CD206','CD1c','CD116','CD126','CD5','CD117','C5AR1','CLEC9A','CD45RA','CD14') 
markers = c('CLEC9A', 'CD14', 'CD1C', 'IL3RA', 'LILRA4', 'FCGR3A', 'FCER1A', 'THBD', 'CD68', 'IGHG4', 'MS4A1', 'KLRD1', 'FCGR3B', 'CD3D', 'TPPP3', 'KRT18', 'CSF2RA', 'CD101', 'MME', 'CXCR4', 'S100A8', 'HLA-DRA', 'HLA-DRB1' )
#markers = c('HLA-DRA')
png(file="marker/violin_marker_selected.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
print(VlnPlot(object = SelectedClusters, features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
dev.off()

png(file="MNP/MNP_umap_marker.png", width = dpi*30, height = dpi*36, units = "px",res = dpi,type='cairo')
print(FeaturePlot(object = SelectedClusters, features = markers,cols = c("lightgrey","#ff0000")))
dev.off()

dpi=300
png(file="MNP/MNP-umap-group-group.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = SelectedClusters, reduction = 'umap',label = FALSE, group.by = 'group')
dev.off()

png(file="MNP/MNP-umap-group-disease.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = SelectedClusters, reduction = 'umap',label = FALSE, group.by = 'disease')
dev.off()

png(file="MNP/MNP-umap-group-sample.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = SelectedClusters, reduction = 'umap',label = FALSE, group.by = 'sample_new')
dev.off()


MNP.Marker.Disease <- FindMarkers(object = SelectedClusters,ident.1 = "Y",group.by = SelectedClusters$disease, assay = 'RNA', test.use = 'bimod', logfc.threshold = 0.25)
MNP.Marker.M_HC <- FindMarkers(object = SelectedClusters,ident.1 = "M",ident.2 = "HC",group.by = SelectedClusters$group, assay = 'RNA', test.use = 'bimod')
MNP.Marker.S_HC <- FindMarkers(object = SelectedClusters,ident.1 = "S",ident.2 = "HC",group.by = SelectedClusters$group, assay = 'RNA', test.use = 'bimod')

write.table(MNP.Marker.Disease,file='MNP/MNP_Marker_Disease_Y_N.txt', row.names = TRUE, quote = FALSE, sep = '\t')
write.table(MNP.Marker.M_HC,file='MNP/MNP_Marker_M_HC.txt', row.names = TRUE, quote = FALSE, sep = '\t')
write.table(MNP.Marker.S_HC,file='MNP/MNP_Marker_S_HC.txt', row.names = TRUE, quote = FALSE, sep = '\t')

#####################Sub clusters MDC
MDC.Clusters <- subset(x = SelectedClusters, idents = c(11,27,29))
# Run the standard workflow for visualization and clustering
MDC.Clusters <- ScaleData(MDC.Clusters, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
MDC.Clusters <- RunPCA(MDC.Clusters, verbose = FALSE,npcs = 100)
MDC.Clusters <- ProjectDim(object = MDC.Clusters)
png(file="MDC/pca_MDC.png", width = dpi*10, height = dpi*6, units = "px",res = dpi,type='cairo')
ElbowPlot(object = MDC.Clusters,ndims = 100)
dev.off()

###cluster
MDC.Clusters <- FindNeighbors(object = MDC.Clusters, dims = 1:50)
MDC.Clusters <- FindClusters(object = MDC.Clusters, resolution = 1.2) 

###tsne and umap
MDC.Clusters <- RunTSNE(object = MDC.Clusters, dims = 1:50)
MDC.Clusters <- RunUMAP(MDC.Clusters, reduction = "pca", dims = 1:50)
png(file="MDC/tsne_MDC.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = MDC.Clusters, reduction = 'tsne',label = TRUE)
dev.off()
png(file="MDC/umap_MDC.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = MDC.Clusters, reduction = 'umap',label = TRUE)
dev.off()

################Draw selected markers on TSNE/UMAP
markers = c('CLEC9A', 'CD14', 'CD1C', 'IL3RA', 'LILRA4', 'FCGR3A', 'FCER1A', 'THBD', 'CD68', 'IGHG4', 'MS4A1', 'KLRD1', 'FCGR3B', 'CD3D', 'TPPP3', 'KRT18', 'CSF2RA', 'CD101', 'MME', 'CXCR4', 'S100A8', 'HLA-DRA', 'HLA-DRB1' )
png(file="MDC/MDC_violin_marker.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
print(VlnPlot(object = MDC.Clusters, features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)))
dev.off()

png(file="MDC/MDC_umap_marker.png", width = dpi*30, height = dpi*36, units = "px",res = dpi,type='cairo')
print(FeaturePlot(object = MDC.Clusters, features = markers,cols = c("lightgrey","#ff0000")))
dev.off()

dpi=300
png(file="MDC/MDC-umap-group-group.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = MDC.Clusters, reduction = 'umap',label = FALSE, group.by = 'group')
dev.off()

png(file="MDC/MDC-umap-group-disease.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = MDC.Clusters, reduction = 'umap',label = FALSE, group.by = 'disease')
dev.off()

png(file="MDC/MDC-umap-group-sample.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = MDC.Clusters, reduction = 'umap',label = FALSE, group.by = 'sample_new')
dev.off()
########List differencial expressed genes
MDC.Marker.Disease <- FindMarkers(object = MDC.Clusters,ident.1 = "Y",group.by = MDC.Clusters$disease, assay = 'RNA', test.use = 'bimod', logfc.threshold = 0.25)
MDC.Marker.M_HC <- FindMarkers(object = MDC.Clusters,ident.1 = "M",ident.2 = "HC",group.by = MDC.Clusters$group, assay = 'RNA', test.use = 'bimod')
MDC.Marker.S_HC <- FindMarkers(object = MDC.Clusters,ident.1 = "S",ident.2 = "HC",group.by = MDC.Clusters$group, assay = 'RNA', test.use = 'bimod')

write.table(MDC.Marker.Disease,file='MDC/MDC_Marker_Disease_Y_N.txt', row.names = TRUE, quote = FALSE, sep = '\t')
write.table(MDC.Marker.M_HC,file='MDC/MDC_Marker_M_HC.txt', row.names = TRUE, quote = FALSE, sep = '\t')
write.table(MDC.Marker.S_HC,file='MDC/MDC_Marker_S_HC.txt', row.names = TRUE, quote = FALSE, sep = '\t')

markers = c('C5AR1')
png(file="figures/violin_C5AR1.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
#print(VlnPlot(object = SelectedClusters,group.by ='group', features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_fill_manual(values=c("#FFFFFF", "#A5CBA1", "#A5CBA1")))
print(VlnPlot(object = SelectedClusters,group.by ='group', features = markers,pt.size = 0)+scale_fill_manual(values=c("#FFFFFF", "#A5CBA1", "#FF6347")))
dev.off()

markers =c('CCL2', 'CCL7', 'CCL8') #4B
#markers=c('S100A9', 'LYZ', 'CD163', 'NAMPT' ,'SRGN', 'LITAF', 'IFI6','MAFB','ISG15', 'LST1', 'TCF7L2', 'MS4A7') #4C 
#markers=c('IFITM1', 'IFITM2','IFITM3','GRINA','HLA-DRA','HLA-DQB1','HLA-DPA1','CD74') #4D
#markers =c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F') #4E

####Not in use
#markers = c('C5AR1','FCGR3A', 'LST1', 'CTSS', 'TCF7L2', 'MS4A7') #4E
#markers = c('IFITM1','IL1B','SOCS3','GRINA','IFITM2','CD84','CXCR4','HLA-DRA','HLA-DQB1','HLA-DPA1','CD74','CAMP')#4F
#markers=c('MPO', 'MKI67', 'MPO', 'PLAC8', 'IL1R2', 'SELL')
#markers=c('ISG15', 'IFI6', 'IFITM3', 'APOBEC3A', 'TYMP', 'ISG15', 'CST3', 'LGALS2', 'IFI6')
#markers=c('CD226', 'CD69', 'CXCR3','MPO', 'ELANE', 'PRTN3', 'CD24', 'LCN2', 'PADI4','ZFP36L2', 'VCAN', 'LYZ', 'S100A9',
#          'IL6R', 'CSF2RB', 'CD38', 'CD44', 'BSG', 'THBD','IL10', 'AREG', 'TGFBI', 'TGFB1', 'VEGFA', 'KITLG')
#markers=c('IL1B', 'CXCR4')
#markers=c('PLBD1', 'MAFB', 'LGALS1', 'CD14', 'FCN1','ISG15', 'LGALS2','COTL1', 'AIF1', 'PSAP', 'CSF', 'FLT', 'CTSS', 'TNFAIP2', 'FGL2')
#markers=c('CXCR2','FCGR3B')
for(marker in markers){
  png(file=paste("figures/Fig4/4B_violin_",marker,".png",sep=''), width = dpi*3, height = dpi*3, units = "px",res = dpi,type='cairo')
  print(VlnPlot(object = MNP.Clusters, group.by='group'
                , features = marker,pt.size =0)+
          theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_fill_manual(values=c("#FFFFFF", "#A5CBA1", "#007475")))
  dev.off()
}  

#markers = c('C5AR1','FCGR3A', 'LST1', 'CTSS', 'TCF7L2', 'MS4A7')
#png(file="figures/Fig4E_violin_markers.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
#print(VlnPlot(object = SelectedClusters, group.by = 'group', features = markers,pt.size = 0)+ theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + scale_fill_manual(values=c("#FFFFFF", "#A5CBA1", "#007475")))
#dev.off()
VlnPlot(object = SelectedClusters, group.by='group'
        , features = marker, pt.size =0)+ 
  theme(axis.text.y = element_text(size = 20),axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_fill_manual(values=c("#FFFFFF", "#A5CBA1", "#007475"))

##cluster 4 and 7 and others
MDC.Cluster.4 <- subset(x = MDC.Clusters, idents = c(4))
MDC.Cluster.7 <- subset(x = MDC.Clusters, idents = c(7))
MDC.Cluster.other <- subset(x = MDC.Clusters, idents = c(0,1,2,3,5,6,8 ))
MDC.Cluster.1358 <- subset(x = MDC.Clusters, idents = c(1,3,5,8 ))

MDC.Cluster.4.S_M <- FindMarkers(object = MDC.Cluster.4,ident.1 = "S", ident.2 = "M", group.by = MDC.Cluster.4$group, assay = 'RNA', test.use = 'bimod')
MDC.Cluster.4.S_Rest <- FindMarkers(object = MDC.Cluster.4, ident.1 = "S", group.by = MDC.Cluster.4$group, assay = 'RNA', test.use = 'bimod')
MDC.Cluster.7.S_HC <- FindMarkers(object = MDC.Cluster.7,ident.1 = "S",ident.2 = "HC",group.by = MDC.Cluster.7$group, assay = 'RNA', test.use = 'bimod')
MDC.Cluster.7.S_Rest <- FindMarkers(object = MDC.Cluster.7,ident.1 = "S",group.by = MDC.Cluster.7$group, assay = 'RNA', test.use = 'bimod')

MDC.Cluster.other.D <- FindMarkers(object = MDC.Cluster.other,ident.1 = "Y",ident.2 = "N",group.by = MDC.Cluster.other$disease, assay = 'RNA', test.use = 'bimod')
MDC.Cluster.1358.S_M <- FindMarkers(object = MDC.Cluster.1358,ident.1 = "S", ident.2 = "M", group.by = MDC.Cluster.1358$group, assay = 'RNA', test.use = 'bimod')

write.table(MDC.Cluster.4.S_M,file='MDC/MDC_Cluster_4_S_M.txt', row.names = TRUE, quote = FALSE, sep = '\t')
write.table(MDC.Cluster.4.S_Rest, file = 'MDC/MDC_Cluster_4_S_Rest.txt', row.names = TRUE, quote = FALSE, sep='\t')
write.table(MDC.Cluster.7.S_HC,file='MDC/MDC_Cluster_7_S_HC.txt', row.names = TRUE, quote = FALSE, sep = '\t')
write.table(MDC.Cluster.7.S_Rest,file='MDC/MDC_Cluster_7_S_Rest.txt', row.names = TRUE, quote = FALSE, sep = '\t')

write.table(MDC.Cluster.other.D,file='MDC/MDC_Cluster_other_Disease.txt', row.names = TRUE, quote = FALSE, sep = '\t')
write.table(MDC.Cluster.1358.S_M,file='MDC/MDC_Cluster_1358_S_M.txt', row.names = TRUE, quote = FALSE, sep = '\t')

# find markers for every cluster compared to all remaining cells, report only the positive ones
MDC.Clusters@misc$markers <- FindAllMarkers(object = MDC.Clusters, assay = 'RNA', only.pos = TRUE,test.use = 'MAST')
write.table(MDC.Clusters@misc$markers,file='marker_MDC_clusters.txt', row.names = FALSE, quote = FALSE, sep = '\t')
MDC.markers = read.delim2("marker_MDC_clusters.txt",header = TRUE, stringsAsFactors = FALSE,check.names = FALSE, sep = "\t")
MDC.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC) -> MDC.top30
tt1 = DoHeatmap(object = subset(MDC.Clusters, downsample = 500), features = MDC.top30$gene) + NoLegend()
ggplot2::ggsave(file="MDC/marker_heatmap_MAST.pdf",plot = tt1,device = 'pdf',width = 20, height = 16, units = "in",dpi = dpi,limitsize = FALSE)


#saveRDS(nCoV.integrated, file = "SaveRDS/ncov_integrated.rds")
#saveRDS(SelectedClusters, file = "SaveRDS/SelectedClusters.rds")
#saveRDS(MDC.Clusters, file = "SaveRDS/MDC.Clusters.rds")

write.table(MDC.top30,file='MDC/MDC_Heatmap_genes.txt', row.names = TRUE, quote = FALSE, sep = '\t')

#clustering wtih 3 colors
png(file="MDC/umap_MDC_3color.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = MDC.Clusters, reduction = 'umap',label = FALSE,
        cols = c('0' = '#C1BDBD','1' = '#C1BDBD','2' = '#C1BDBD',
                 '3' = '#C1BDBD','4' = '#27BA9C', '5' = '#C1BDBD',
                 '6' = '#C1BDBD','7' = '#F3756D', '8' = '#C1BDBD'))
dev.off()

markers=c('CLEC9A', 'CADM1')
png(file="MDC/MDC_umap_CLEC9A_CADM1.png", width = dpi*30, height = dpi*36, units = "px",res = dpi,type='cairo')
print(FeaturePlot(object = SelectedClusters, features = markers,cols = c("lightgrey","#ff0000")))
dev.off()
#Fig 4F
markers=c('CLEC9A', 'CADM1','LILRA4', 'IL3RA')
for(marker in markers){
  png(file=paste("MDC/4F_",marker,".png",sep=''), width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
  #print(FeaturePlot(object = MDC.Clusters, features = marker,cols = c("lightgrey","red")))
  print(VlnPlot(object = MDC.Clusters, features = marker,pt.size =0)+
          theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
        +scale_fill_manual(values=c('0' = '#C1BDBD','1' = '#C1BDBD','2' = '#C1BDBD',
                                    '3' = '#C1BDBD','4' = '#27BA9C', '5' = '#C1BDBD',
                                    '6' = '#C1BDBD','7' = '#F3756D', '8' = '#C1BDBD')))
  
  dev.off()
}

##Supplimentory figure 5B
markers=c('CD14', 'CD68','CLEC9A', 'CADM1', 'CD1C', 'FCER1A','LILRA4'
          , 'IL3RA','FCGR3B', 'CXCR2', 'IGHG4','MS4A1','TPPP3', 'KRT18','CD3D','KLRD1')
markers=c('CD14', 'CD68', 'IGHG4','MS4A1','CLEC9A','LILRA4', 'CD1C', 'FCER1A'
          ,'TPPP3', 'KRT18','CD3D','KLRD1')

#markers=c('CXCR2','FCGR3B')
for(marker in markers){
  png(file=paste("figures/supp_fig/2B_umap_",marker,".png",sep=''), width = dpi*6, height = dpi*6, units = "px",res = dpi,type='cairo')
  print(FeaturePlot(object = nCoV.integrated, features = marker,label.size = 6, cols = c("lightgrey","#ff0000"))+
          theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20)
                , axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)
                , legend.text = element_text(size = 20))
  )
  dev.off()
}

markers=c('CD14', 'CD68','CLEC9A', 'CADM1', 'CD1C', 'FCER1A','LILRA4', 'IL3RA','CXCR2','FCGR3B')
for(marker in markers){
  png(file=paste("figures/supp_fig/2G_umap_",marker,".png",sep=''), width = dpi*6, height = dpi*6, units = "px",res = dpi,type='cairo')
  print(FeaturePlot(object = SelectedClusters, features = marker,label.size = 6,cols = c("lightgrey","#ff0000"))+
          theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20)
                , axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)
                , legend.text = element_text(size = 20))
        )
  dev.off()
}

#Data for distribution stack bar plot
stack <- table(Idents(nCoV.integrated),nCoV.integrated$sample)
stack <- table(Idents(nCoV.integrated),nCoV.integrated$disease)
stack <- table(Idents(nCoV.integrated),nCoV.integrated$group)

stack <- table(Idents(SelectedClusters),SelectedClusters$sampleID)
stack <- table(Idents(SelectedClusters),SelectedClusters$disease)
stack <- table(Idents(SelectedClusters),SelectedClusters$group)

stack <- table(Idents(MDC.Clusters),MDC.Clusters$sampleID)
stack <- table(Idents(MDC.Clusters),MDC.Clusters$disease)
stack <- table(Idents(MDC.Clusters),MDC.Clusters$group)

stack
#Stack bar plots
#stackplot <- read.table("figures/stack bar graphs/cell population per sample.txt", sep="\t", header=T, quote = '')
#stackplot
# Stacked + percent
#ggplot(stackplot, aes(fill=C100, y=C100, x=ID)) + 
#  geom_bar(position="fill", stat="identity")
#marker="CCL2"
markers=c('AREG','CCL2','NEAT1','FOS','CD83','NR4A3','BCL2A1','ID2','KLF6','IFI27','TSC22D3','CXCR4','JUN','CD55','TGFB1','RHOG','ITGA4','LAMP5','TXNIP','CLEC4C','HLA-DQA2')
markers=c('IFITM3','IFITM1','TAOK1','IFI27','IFI30','MT-ATP8','MT-ND4L','FOS','IFITM2','MX1','TNFSF10','NR4A3','IFI44L','IRF7','NME2','CD83','DUSP4','ISG15','DUSP1','COX5A','CCL2','SAMD9L','SAMD9','ALDOA','CXCL10','GNG5','VRK2','CD74','LYZ')
for(marker in markers){
  png(file=paste("MDC/IPA results/DC1/violin_",marker,".png",sep=''), width = dpi*3, height = dpi*3, units = "px",res = dpi,type='cairo')
  print(VlnPlot(MDC.Cluster.7, group.by='group'
                , features = marker,pt.size =0)+
          theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_fill_manual(values=c("#FFFFFF", "#A5CBA1", "#007475")))
  dev.off()
}  


ggplot(MDC.Cluster.4.S_M, aes(x=marker, y=avg_logFC, fill=avg_logFC)) + 
  geom_violin()

dpi=300
png(file="figures/umap_MNP_New.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
DimPlot(object = MNP.Clusters, reduction = 'umap',label = TRUE, label.size = 6)+
theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(size = 20)
      , axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)
      , legend.text = element_text(size = 20))
dev.off()

MNP.Clusters <- subset(x = SelectedClusters, idents = c(1:14,16:31))
###############23-06-2020

MDC.Cluster.4.S_Rest <- FindMarkers(object = MDC.Cluster.4, ident.1 = "S", group.by = MDC.Cluster.4$group, assay = 'RNA', test.use = 'bimod')
write.table(MDC.Cluster.4.S_Rest, file = 'MDC/MDC_Cluster_4_S_Rest.txt', row.names = TRUE, quote = FALSE, sep='\t')
MDC.Cluster.7.S_Rest <- FindMarkers(object = MDC.Cluster.7,ident.1 = "S",group.by = MDC.Cluster.7$group, assay = 'RNA', test.use = 'bimod')
write.table(MDC.Cluster.7.S_Rest,file='MDC/MDC_Cluster_7_S_Rest.txt', row.names = TRUE, quote = FALSE, sep = '\t')
MNP.Clusters_M_HC <- FindMarkers(object = MNP.Clusters,ident.1 = "M", ident.2 = "HC", group.by = MNP.Clusters$group, assay = 'RNA', test.use = 'bimod')
write.table(MNP.Clusters_M_HC,file = 'MNP/MNP.Clusters_M_HC.txt',row.names = TRUE, quote = FALSE, sep = '\t')
MNP.Clusters_S_HC <- FindMarkers(object = MNP.Clusters,ident.1 = "S", ident.2 = "HC", group.by = MNP.Clusters$group, assay = 'RNA', test.use = 'bimod')
write.table(MNP.Clusters_S_HC,file = 'MNP/MNP.Clusters_S_HC.txt',row.names = TRUE, quote = FALSE, sep = '\t')
MNP.Clusters_S_M <- FindMarkers(object = MNP.Clusters,ident.1 = "S", ident.2 = "M", group.by = MNP.Clusters$group, assay = 'RNA', test.use = 'bimod')
write.table(MNP.Clusters_S_M,file = 'MNP/MNP.Clusters_S_M.txt',row.names = TRUE, quote = FALSE, sep = '\t')

##############26-06-2020### gene names with HLA
markers=c('HLA-A','HLA-B','HLA-C','HLA-DMA','HLA-DMB','HLA-DOA','HLA-DOB','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DQB2','HLA-DRA','HLA-DRB1','HLA-DRB5','HLA-E','HLA-F')
for(marker in markers){
  png(file=paste("MDC/HLA/other.Cluster_",marker,".png",sep=''), width = dpi*3, height = dpi*3, units = "px",res = dpi,type='cairo')
  print(VlnPlot(MDC.Cluster.other, group.by='group'
                , features = marker,pt.size =0)+
          theme(axis.text.y = element_text(size = 20), axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+scale_fill_manual(values=c("#FFFFFF", "#A5CBA1", "#007475")))
  dev.off()
}  

