Args <- commandArgs()
print(Args)

library('Seurat')
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)

varGenes <- as.integer(Args[6])
CCAdim <- as.integer(Args[7])
PCAdim <- as.integer(Args[8])
#TSNE/UMAP
DimRed <- Args[9]



dataDir <- '/Path/to/scRNAseq/rawData'
PATH <- '/Path/to/Outdir'
SCN1_dir <- paste0(dataDir, '/SCN1/outs/filtered_gene_bc_matrices/hg19')
SCP1_dir <- paste0(dataDir, '/SCP1/outs/filtered_gene_bc_matrices/hg19')
SCP2_dir <- paste0(dataDir, '/SCP2/outs/filtered_gene_bc_matrices/hg19')


outDir <- paste0(PATH,'/varGenes',Args[6],'_CCA',Args[7],'_PCA',Args[8],'_DimRed',Args[9])
if (!dir.exists(outDir)){
  dir.create(outDir)}
FilePrefix=paste0(outDir,'/Aggr3')

SCN1.data <- Read10X(SCN1_dir)
SCP1.data <- Read10X(SCP1_dir)
SCP2.data <- Read10X(SCP2_dir)

SCN1 <- CreateSeuratObject(counts = SCN1.data, project = "SCN1", min.cells = 10, min.features = 400)
SCP1 <- CreateSeuratObject(counts = SCP1.data, project = "SCP1", min.cells = 10, min.features = 400)
SCP2 <- CreateSeuratObject(counts = SCP2.data, project = "SCP2", min.cells = 10, min.features = 400)


SLE.big <- merge(SCN1, y = c(SCP1, SCP2), add.cell.ids = c("SCN1", "SCP1", "SCP2"), project = "SLE")

mito.genes=grep(pattern='^MT-',x=rownames(x=GetAssayData(object = SLE.big)),value=TRUE)
percent.mito=Matrix::colSums(GetAssayData(object = SLE.big, slot = "counts")[mito.genes,])/Matrix::colSums(GetAssayData(object = SLE.big, slot = "counts"))
SLE.big=AddMetaData(object=SLE.big,metadata=percent.mito,col.name='percent.mt')
SLE.big <- subset(SLE.big, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 5)





SLE.list <- SplitObject(SLE.big, split.by = "orig.ident")
for (i in 1:length(SLE.list)) {
  SLE.list[[i]] <- NormalizeData(SLE.list[[i]], verbose = FALSE)
  SLE.list[[i]] <- FindVariableFeatures(SLE.list[[i]], selection.method = "vst", nfeatures = varGenes, 
                                             verbose = FALSE)
#Remove chrY Gene
  YGene=c('DDX3Y','TTTY15','EIF1AY','RPS4Y1','XGY2','SRY','RPS4Y1','ZFY','TGIF2LY','PCDH11Y','TTTY23B','TTTY23','TSPY2','FAM197Y9','LINC00280','TTTY1B','TTTY1','TTTY2B','TTTY2','TTTY21B','TTTY21','TTTY7B','TTTY7','TTTY8B','TTTY8','AMELY','TBL1Y','TTTY16','TTTY12','LINC00279','TTTY18','TTTY19','TTTY11','RBMY1A3P','TTTY20','TSPY10','TSPY3','TSPY4','TSPY8','FAM197Y2','FAM197Y4','FAM197Y7','FAM197Y8','FAM197Y6','TSPY1','FAM197Y3','RBMY3AP','TTTY22','GYG2P1','TTTY15','USP9Y','UTY','TMSB4Y','VCY1B','VCY','NLGN4Y','NLGN4Y-AS1','FAM41AY1','FAM41AY2','FAM224B','FAM224A','XKRY2','XKRY','CDY2B','CDY2A','HSFY1','HSFY2','TTTY9A','TTTY9B','TTTY14','CD24P4','LOC107987345','TXLNGY','KDM5D','TTTY10','EIF1AY','RPS4Y2','ERVH-6','PRORY','RBMY2EP','RBMY1A1','RBMY1B','RBMY1D','RBMY1E','TTTY13','PRY2','PRY','LOC101929148','TTTY6','TTTY6B','RBMY1F','RBMY1J','TTTY5','RBMY2FP','LOC100652931','TTTY17A','TTTY17B','TTTY17C','TTTY4B','TTTY4C','TTTY4','BPY2B','BPY2C','BPY2','DAZ1','DAZ2','DAZ3','DAZ4','TTTY3B','TTTY3','CDY1','CDY1B','CSPG4P1Y','GOLGA2P2Y','GOLGA2P3Y','PRYP4')
  ordering_genes=VariableFeatures(object = SLE.list[[i]])
  print('Length of ordering_genes:')
  print(length(ordering_genes))
  print('Remove ChrY gene!')
  for (g in VariableFeatures(object = SLE.list[[i]])){
    if (g %in% YGene){
      print(g)
      ordering_genes=ordering_genes[ordering_genes!=g]
    }
  }
  VariableFeatures(object = SLE.list[[i]])=ordering_genes
}

reference.list <- SLE.list[c("SCN1", "SCP1", "SCP2")]
SLE.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:CCAdim)
print('FindIntegrationAnchors')
SLE.integrated <- IntegrateData(anchorset = SLE.anchors, dims = 1:CCAdim)
DefaultAssay(SLE.integrated) <- "integrated"
print('IntegrateData')
#top50 <- head(VariableFeatures(SLE.integrated), 50)
#pdf(file=paste0(FilePrefix,'VarGene.pdf'),width=20, height=7)
#plot1 <- VariableFeaturePlot(SLE.integrated)
#plot2 <- LabelPoints(plot = plot1, points = top50, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))
#dev.off()

SLE.integrated <- ScaleData(SLE.integrated, verbose = FALSE)
SLE.integrated <- RunPCA(SLE.integrated, npcs = 40, verbose = FALSE)
print("Run PCA Done!")

pdf(file=paste0(outDir,'/PCA.pdf'))
ElbowPlot(object = SLE.integrated, ndims = 45)
dev.off()

SLE.integrated <- FindNeighbors(SLE.integrated, dims = 1:PCAdim)
SLE.integrated <- FindClusters(SLE.integrated, resolution = 1)
print("FindClusters Done!")

PCA=as.data.frame(Embeddings(object = SLE.integrated, reduction = "pca"))
write.table(PCA,file=paste0(FilePrefix,'.pca.txt'),sep='\t',quote=F)

if (DimRed=='TSNE'){
print('TSNE')
SLE.integrated <- RunTSNE(SLE.integrated, dims = 1:PCAdim, reduction = "pca")
pdf(file=paste0(outDir,'/Cluster.pdf'))
DimPlot(SLE.integrated, reduction = "tsne")
dev.off()

pdf(file=paste0(outDir,'/Batch.pdf'))
DimPlot(SLE.integrated, reduction = "tsne", group.by = "orig.ident")
dev.off()
TSNE=as.data.frame(Embeddings(object = SLE.integrated, reduction = "tsne"))
write.table(TSNE,file=paste0(FilePrefix,'.tsne.txt'),sep='\t',quote=F)
print("Save TSNE Done!")

}else if (DimRed=='UMAP'){
  print("UMAP")
  SLE.integrated <- RunUMAP(SLE.integrated, dims = 1:PCAdim, reduction = "pca")
  pdf(file=paste0(outDir,'/Cluster.pdf'))
  DimPlot(SLE.integrated, reduction = "umap")
  dev.off()
  
  pdf(file=paste0(outDir,'/Batch.pdf'))
  DimPlot(SLE.integrated, reduction = "umap", group.by = "orig.ident")
  dev.off()  
  TSNE=as.data.frame(Embeddings(object = SLE.integrated, reduction = "umap"))
  write.table(TSNE,file=paste0(FilePrefix,'.umap.txt'),sep='\t',quote=F)
  print("Save UMAP Done!")
}

#Save MetaData
MetaData=as.data.frame(SLE.integrated@meta.data)
write.table(MetaData,file=paste0(FilePrefix,'.MetaData.txt'),sep='\t',quote=F)


pdf(paste0(outDir,'/FeaturePlot1.pdf'),width=12, height=10)
FeaturePlot(
  SLE.integrated,
  features = c('FOXP3','TIGIT','RTKN2','TBX21','CCR5','CCL5','CXCR3','IFNG')
)
dev.off()

pdf(paste0(outDir,'/FeaturePlot2.pdf'),width=12, height=10)
FeaturePlot(
  SLE.integrated,
  features = c('IL17A','CCR6','RORA','RORC','SELL','CCR7')
)
dev.off()

pdf(paste0(outDir,'/FeaturePlot3.pdf'),width=12, height=10)
FeaturePlot(
  SLE.integrated,
  features = c('NKG7','FGFBP2', 'GZMA','GZMB','GZMK','GNLY')
)
dev.off()

#SaveExpData
Data=as.data.frame(as.matrix(GetAssayData(object = SLE.integrated)))
write.table(Data,file=paste0(FilePrefix,'.dataNorm.txt'),sep='\t',quote=F)

VarGenes=VariableFeatures(object = SLE.integrated)
VarGeneData=Data[VarGenes,]
write.table(VarGeneData,file=paste0(FilePrefix,'.VarGeneData.txt'),sep='\t',quote=F)


#Marker Gene
#pbmc.markers=FindAllMarkers(object=SLE.integrated,test.use = "wilcox",only.pos=TRUE,min.pct=0.1,return.thresh=0.01,logfc.threshold=0.2)

#MarkerGene=as.data.frame(pbmc.markers %>% group_by(cluster) %>% top_n(30,avg_logFC))
#write.table(MarkerGene,file=paste(FilePrefix,'_30MarkerGene.wilcox.txt'),quote =F,sep='\t')
#print("Save MarkerGene File Done!")


#top10=pbmc.markers %>% group_by(cluster) %>% top_n(15,avg_logFC)
#pdf(file=paste0(FilePrefix,'_10MarkerGene.png'))
#DoHeatmap(object=SLE.integrated,features=top10$gene)+ NoLegend()
#dev.off()
#print("Save MarkerGene Heatmap Done!")


#nUMI 
pdf(file=paste0(FilePrefix,'.nUMI.Count.pdf'))
MetaData$UMAP1=TSNE$UMAP_1
MetaData$UMAP2=TSNE$UMAP_2
ggplot(MetaData,aes(x=UMAP1,y=UMAP2,colour=nUMI))+ geom_point()
dev.off()




