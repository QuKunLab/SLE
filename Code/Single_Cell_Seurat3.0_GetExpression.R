Args <- commandArgs()
print(Args)

library('Seurat')
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(SAVER)



dataDir <- '/Path/to/scRNAseq/rawData'
PATH <- '/Path/to/Outdir'
SCN1_dir <- paste0(dataDir, '/SCN1/outs/filtered_gene_bc_matrices/hg19')
SCP1_dir <- paste0(dataDir, '/SCP1/outs/filtered_gene_bc_matrices/hg19')
SCP2_dir <- paste0(dataDir, '/SCP2/outs/filtered_gene_bc_matrices/hg19')


FilePrefix=paste0(PATH,'/Aggr3')

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
SLE.big <- NormalizeData(SLE.big, normalization.method = "LogNormalize", scale.factor = 10000)




SLE.big=NormalizeData(SLE.big, verbose = FALSE)



#SaveExpData
Data=as.data.frame(as.matrix(GetAssayData(object = SLE.big)))
write.table(Data,file=paste0(FilePrefix,'.dataNorm.txt'),sep='\t',quote=F)

#Do SAVER Analysis
print("Do SAVER Analysis!")
save=saver(Data,size.factor = 1)
print("Do SAVER Analysis Done!")
DataSaver=as.data.frame(save$estimate)
write.table(DataSaver,file=paste0(FilePrefix,'.dataSAVER.txt'),quote=F,sep='\t')















