# ##################Single cell analysis
library(data.table)
library("R.utils")
library(dyno)
library(irlba)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(RColorBrewer)
library(ggplot2)

###Build the file
pca1 <- fread("Input5.GSE86469.unique.txt",
              data.table = F)
pca1[1:4,1:4]
d1=pca1

d1 <- aggregate(.~SYMBOL,d1,max)
rownames(d1)=d1[,1]

save(d1,file="GSE86469.rds")



pbmc1 <- CreateSeuratObject(counts = d1,
                            min.cells = 3, 
                            min.features = 200,
                            project = "GSE86469")


pbmc =pbmc1

# pbmc = pbmc1
as.data.frame(pbmc@assays$RNA@counts[1:10, 1:2])
head(pbmc@meta.data)


save(pbmc,file="pbmc.GSE86469.rds")  #save

load("pbmc.GSE86469.rds")

###################################01.filtration###################################
logFCfilter=1               #logFC
adjPvalFilter=0.05          #adjPval


#使用PercentageFeatureSet Percentage of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc=subset(x = pbmc, subset = nFeature_RNA > 100 & percent.mt < 5 & nCount_RNA > 3 )   
#Plotting a violin of genetic characteristics
pdf(file="01.featureViolin.pdf", width=10, height=6)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()



#Correlation plot of sequencing depth
pdf(file="01.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()



#Standardize data
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#Genes with a large coefficient of variation between cells were extracted
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
#Output a characteristic variance plot
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()



###################################02.PCA###################################
##PCA
pbmc=ScaleData(pbmc)          #Standard preprocessing steps prior to PCA dimensionality reduction
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     

#Map the characteristic genes for each PCA component
pdf(file="02.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#Graph principal component analysis
pdf(file="02.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#heatmap
pdf(file="02.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#P-value distribution
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()

pdf(file="02.ElbowPlot.pdf",width=8,height=6)
ElbowPlot(pbmc)
dev.off()



###################################03.TSNE###################################
##TSNE
pcSelect=8
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)       #Calculate the adjacency distance
pbmc <- FindClusters(object = pbmc, resolution = 0.5)         #Group cells and standardize cell modularity
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             #TSNE
pdf(file="03.TSNE.pdf",width=6.5,height=6)
# TSNEPlot(object = pbmc, pt.size = 2, label = TRUE,group.by = "groups")    #visualization TSNE
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE) 
dev.off()
write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)




pdf(file="03_Cells_cluster.pdf",width=10,height=8)
DimPlot(pbmc, reduction = "umap",split.by = "groups",label = TRUE, repel = FALSE)
print(p)
dev.off()


##Find the differential gene for each cluster
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.2,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#Draw a heatmap of the marker in each cluster
pdf(file="03.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

#Draw marker
pdf(file="03.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = row.names(sig.markers)[1:2])
dev.off()


showGenes=c("TM4SF4","ADCYAP1","CDKN1C","AKAP12","SOD2","PNLIP","COL4A1","SPARC") 


#Draw a scatter plot of the marker in each cluster
pdf(file="03.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = showGenes, cols = c("green", "red"))
dev.off()

#Draw a bubble chart of marker in each cluster
pdf(file="03.markerBubble.pdf",width=12,height=6)
cluster10Marker=showGenes
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()




###################################04.SingleR###################################
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
#ref=get(load("ref_Human_all.RData"))
ref=celldex::HumanPrimaryCellAtlasData()
singler=SingleR(test=counts, ref =ref,
                labels=ref$label.main, clusters = clusters)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
singler2=SingleR(test=counts, ref =ref, 
                 labels=ref$label.main)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#Visualization after cluster annotation
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    
dev.off()

##Variance analysis after cluster annotation
pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)





###################################05.monocle###################################
#Prepare the documents required for cell trajectory analysis
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

save(monocle.matrix,file="monocle.matrix.rds")
save(monocle.sample,file="monocle.sample.rds")
save(monocle.geneAnn,file="monocle.geneAnn.rds")
save(monocle.clusterAnn,file="monocle.clusterAnn.rds")
save(monocle.markers,file="monocle.markers.rds")

load("monocle.matrix.rds")
load("monocle.sample.rds")
load("monocle.geneAnn.rds")
load("monocle.clusterAnn.rds")
load("monocle.markers.rds")



monocle.sample <- read.table('monocle.sample11.txt',
                             header = T,sep = '\t',check.names = F,row.names = 1)

ID <- rownames(monocle.sample)
monocle.matrix <- monocle.matrix[,ID]


#Convert the Seurat results into the cell matrix required by the monocle
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])



#Add cell clustering data
clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
clusterAnn <- clusterAnn[2:6]
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

# pData(cds) <- subset(pData(cds), cell_type2 == "Neurons"||cell_type2 == "Epithelial_cells"||cell_type2 == "Smooth_muscle_cells")
save(cds,file="cds.last.rds")

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)


save(cds,file="cds.CT.rds")

load("cds.CT.rds")

##Cell dimensionality reduction

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')



# Sort the cells
cds <- orderCells(cds, reverse = F) 
pdf(file="05.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
pdf(file="05.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()

# View cell annotation information
head(cds@phenoData@data)
pdf(file="05.trajectory.Cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Cluster")
dev.off()



pdf(file="05.trajectory.3333groups.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "group")
dev.off()
pdf(file="05.trajectory.222cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()

HSMM_expressed_genes <-  row.names(subset(fData(cds)))
HSMM_filtered <- cds[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("TTR", "SST", "CDKN1C","DLK1","SDC4")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "State")
plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")

pdf(file="06..3333cellType.pdf",width=6.5,height=6)
plot_genes_in_pseudotime(cds_subset, color_by =  "group")
dev.off()



############cellphonedb
/data2/project-YangHM/anaconda/bin/cellphonedb method statistical_analysis --counts-data ensembl --output-path test_output meta.txt GSM4006644_BC1_gene_cell_exprs_table.txt


#Bubble chart
cellphonedb plot dot_plot --pvalues-path out/pvalues.txt --output-path test_output  --output-name out.dotplot.pdf

#Heatmap
cellphonedb plot heatmap_plot --pvalues-path out/pvalues.txt --output-path test_output --pvalue 0.05 --count-name test.heatmap_count.pdf --log-name test.heatmap_log_count.pdf --count-network-name test.count_network.txt --interaction-count-name test.interaction_count.txt meta.txt




############Analysis of single-cell group-to-group differences
###keys/control

exp =read.table("Input5.GSE81608.unique.txt", header=T, sep="\t", check.names=F)

namess <- read.table("name.txt", header=T, sep="\t", check.names=F)

colnames(exp) <- namess$Acc


control =read.table("control.GSE81608.test.txt", header=T, sep="\t", check.names=F)

keys =read.table("keys.GSE81608.test.txt", header=T, sep="\t", check.names=F)

control.exp <- exp[,c("SYMBOL",control$Acc)]
keys.exp <- exp[,c("SYMBOL",keys$Acc)]



write.table(control.exp,file="control.GSE81608.test.txt",sep="\t",row.names=F,quote=F)

write.table(keys.exp,file="keys.GSE81608.test.txt",sep="\t",row.names=F,quote=F)




library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(monocle)


#Data from the control group were read and collated
rt=read.table("control.GSE81608.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
conData=avereps(data)
colnames(conData)=paste0("C.", colnames(conData))

#Read the data from the experimental group and collate the data
rt=read.table("keys.GSE81608.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
treatData=avereps(data)
colnames(treatData)=paste0("T.", colnames(treatData))


#Data merging
sameGene=intersect(row.names(conData), row.names(treatData))
data=cbind(conData[sameGene,], treatData[sameGene,])

write.table(data,file="data.test.txt",sep="\t",row.names=T,quote=F)


pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_")

pbmc[["percent.mt"]]=PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmcCon=subset(x = pbmc, subset = nFeature_RNA > 100 & percent.mt < 5 & nCount_RNA > 3 )    



pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)


#Analysis of differences between groups
logFCfilter=0.5
adjPvalFilter=0.05
groups=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(groups)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
pbmc.markers=FindMarkers(pbmc, ident.1 = "T", ident.2 = "C", group.by = 'group', min.pct = 0.25)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val))<adjPvalFilter),]
sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
write.table(sig.markers,file="diffGene.GSE81608.test.txt",sep="\t",row.names=F,quote=F)


################################################################Transcriptome difference analysis


####Keys/Control

gene_exp <- read.table('Input5.GSE86468.unique.txt',
                       sep = '\t',check.names = F,header = T)

exprSet <-gene_exp

gene_exp1 <- aggregate(.~SYMBOL,gene_exp,max)
rownames(gene_exp1) <- gene_exp1$SYMBOL
gene_exp2 <- as.data.frame(gene_exp1[,-1])

nNormal <- 15
nTumor <- 9


ID <- colnames(gene_exp2)
dim(gene_exp2)
Type <- c(rep('Control',nNormal),rep('RPL',nTumor))
ID <- cbind(ID,Type)
colnames(ID) <- c("id","Type")
ID <- as.data.frame(ID)

gene_exp3 <- t(gene_exp2)
id <- rownames(gene_exp3)
gene_exp3 <- cbind(id,gene_exp3)



riskTime <- merge(ID,gene_exp3,by = 'id')
rownames(riskTime) <- riskTime$id
riskTime <- riskTime[ID$id,]


group_list= riskTime$Type %>% factor(.,levels = c("RPL","Control"),ordered = F)
pvalue <-0.05
logFoldChange <- 0.5 #Adjust the difference value


design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <-rownames(riskTime)


gene_exp1[is.na(gene_exp1)] <- 0
dat1 <- apply(gene_exp1[,-1],2,function(x){log2(x+1)})
str(dat1)



dat1 <- dat1[,rownames(design)]

contrast.matrix <- makeContrasts(RPL-Control, levels = design)#Case比Control
fit <- lmFit(dat1, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)


allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.table(cbind(Symbol=rownames(allDiff),allDiff),file="Tumor-Normal.edgeR.Out.txt",sep="\t",row.names = F,quote = F)
diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffSig),diffSig),file="Tumor-Normal.diffSig.0.5.txt",sep="\t",row.names = F,quote = F)



#Volcano map
library(ggplot2)
pvalue <-0.05
logFC <- 0.5
allDiff$Significant <- ifelse(allDiff$P.Value<pvalue &
                                abs(allDiff$logFC)>= logFC,
                              ifelse(allDiff$logFC> logFC,'up','down'),'no')
mycol <- c("#03A9F4","#212121","#E91E63")###
pdf(file="Tumor-Normal.pdf",width=8,height=8)
p <- ggplot(allDiff, aes(logFC, -log10(P.Value), colour= Significant))+
  geom_point(size=1.2,alpha=0.4)+theme_bw()+
  scale_color_manual(values = mycol,name='Significant')+
  labs(title="BLCA vs. Control",x="Log2FC",y="-log10 (P.value)")+
  geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 1)+
  geom_vline(xintercept = c(-logFC, logFC), linetype=3,lwd = 1)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
plot(p)
dev.off()



#Heatmap
library(pheatmap)
diff <- rownames(diffSig)
diffexp <- dat1[diff,]


annotation_col <- data.frame(Type = factor(group_list, levels = c("Control","RPL")))
rownames(annotation_col) <- rownames(riskTime)

color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")

# dev.new()
pdf(file = 'Tumor-Normal.pdf',width = 10,height = 15) 
pheatmap(diffexp1,cellwidth = 8,cellheight = 1,
         method="spearman", #"pearson" (default), "kendall", or "spearman"
         scale="row", #
         cluster_rows=T,#
         cluster_cols=F,#
         color = colorRampPalette(color.key)(200),###
         show_colnames=F,show_rownames =F,
         annotation_col = annotation_col,
         treeheight_row = "1",treeheight_col = "1",#
         border_color = "NA")

dev.off()





#GO/KEGG
library(org.Hs.eg.db) 
library(clusterProfiler)
library(dplyr) 
library(ggplot2)

gene <- gene[,1]
# gene

hg<-bitr(gene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")
id <- hg[,2]
id


###GO
go <- enrichGO(gene,OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'SYMBOL')
dim(go)

write.csv(go,file="go.csv")



#Visualize the results
barplot(go,showCategory=20,drop=T)

dotplot(go,showCategory=20)

pdf(file = "GO.pdf", width = 10, height = 10)
barplot(go, split="ONTOLOGY")+ facet_grid(ONTOLOGY~.,scale="free")
dev.off()

##KEGG
ego <- enrichKEGG(
  gene = id,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 1
)

write.csv(ego,file="KEGG.csv")

pdf(file = "KEGG.pdf", width = 10, height = 10)
dotplot(ego,font.size=8)	
dev.off()





########################ROC curve

library("pROC")
library(ggpubr)
library(ggplot2)
library(pROC)

load('exprSet.all.rda')


gene <- read.table('result.txt',
                   sep = '\t',check.names = F,header = F)

gene_exp <- gene_exp1[gene$V1,]


gene_exp_EDG <- t(gene_exp[,-1])
# ene_exp_EDG <- scale(gene_exp_EDG)

trait <- read.table('GSE86468.sample.txt',
                    sep = '\t',check.names = F,header = T)

colnames(gene_exp)
gene_exp_EDG <- as.data.frame(gene_exp_EDG)
gene_exp_EDG$ID <- rownames(gene_exp_EDG)

trait2 <- merge(trait,gene_exp_EDG,by.x="Acc",by.y = "ID")
# trait2 <- cbind(trait,gene_exp_EDG)
rownames(trait2) <- trait2$Acc
Key <- as.factor(trait2$Type)
trait2 <- cbind(Key,trait2)
str(trait2)

trait2 <- as.data.frame(trait2[,-(2:3)])

write.table(trait2,file="GSE86468.corr.txt",sep="\t",row.names = T,quote = F)



####################Validation set
gene_exp1 <- read.table('Input5.GSE81608.unique.txt',
                        sep = '\t',check.names = 1,header = T)

rownames(gene_exp1) <- gene_exp1$SYMBOL

gene <- read.table('result.txt',
                   sep = '\t',check.names = F,header = F)


gene_exp = gene_exp1 %>% filter(gene_exp1[,1] %in% gene$V1)

gene_exp <- aggregate(.~SYMBOL,gene_exp,max)
rownames(gene_exp) <- gene_exp$SYMBOL

# gene_exp <- gene_exp1[gene$V1,]


gene_exp_EDG <- t(gene_exp[,-1])
gene_exp_EDG <- log2(gene_exp_EDG+1)

# gene_exp_EDG <- scale(gene_exp_EDG)


trait <- read.table('GSE81608.sample.txt',
                    sep = '\t',check.names = F,header = T)

colnames(gene_exp)

trait2 <- cbind(trait,gene_exp_EDG)
Key <- as.factor(trait2$Type)
trait2 <- cbind(Key,trait2)
str(trait2)

trait2 <- as.data.frame(trait2[,-(2:3)])

name <- read.table('GSE81608.names.txt',
                   sep = '\t',check.names = F,header = T)
trait2 <- trait2[name$Acc,]
tezheng.gene <- trait2




rt_roc <- tezheng.gene[,c("Key","CDKN1C")]
colnames(rt_roc) <- c("group","gene")
roc_pt_CDKN1C <- roc(rt_roc$group, rt_roc$gene,levels=c("RPL", "Control"))
AUC_val_CDKN1C <- round(as.numeric(roc_pt_CDKN1C$auc), digits = 3)#AUC

rt_roc <- tezheng.gene[,c("Key","DLK1")]
colnames(rt_roc) <- c("group","gene")
roc_pt_DLK1 <- roc(rt_roc$group, rt_roc$gene,levels=c("RPL", "Control"))
AUC_val_DLK1 <- round(as.numeric(roc_pt_DLK1$auc), digits = 3)#AUC


rt_roc <- tezheng.gene[,c("Key","SDC4")]
colnames(rt_roc) <- c("group","gene")
roc_pt_SDC4 <- roc(rt_roc$group, rt_roc$gene,levels=c("RPL", "Control"))
AUC_val_SDC4 <- round(as.numeric(roc_pt_SDC4$auc), digits = 3)#AUC



pdf("ROC.pdf",width = 8,height = 8)
gl <- ggroc(list(DLK1=roc_pt_DLK1,CDKN1C=roc_pt_CDKN1C),legacy.axes = TRUE, size = 1)


gl+ 
  xlab("1- specificities") + ylab("sensitivities") + 
  ggtitle("GSE86468 ROC")+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="red", linetype="dashed")+
  annotate(geom="text", x = 0.75, y = 0.4, label = paste ('DLK1    AUC = ' , AUC_val_DLK1, sep = "" )) +
  annotate(geom="text", x = 0.75, y = 0.3, label = paste ('CDKN1C    AUC = ' , AUC_val_CDKN1C, sep = "" )) +
  
  
  theme_bw()
dev.off()




#####################################################Correlation analysis


library(ggstatsplot)

corr <- read.table('GSE86468.corr.txt',
                   sep = '\t',check.names = F,header = T)

rownames(corr) <- corr[,1]
corr <- corr[,-1]


pdf('CDKN1C.age.cor.pdf',width = 10,height = 10)
ggscatterstats(corr, 
               y ="CDKN1C", 
               x ="age",
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "densigram",  
               title = "Relationship between CDKN1C and age")

dev.off()




####################Single-gene GSEA analysis
# Samples are grouped by the median of the gene expression value
library ("DESeq2")
library(dplyr)
library(GSEABase)
library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)

library(tidyr)
library(dplyr)
library(stringi)
library(stringr)
library(limma)
library(edgeR)



load('exprSet.all.rda')

gene_exp <- gene_exp1[,-1]

min(gene_exp)
max(gene_exp)

gene_exp <- apply(gene_exp[,-1],2,function(x){log2(x+1)})

gene <- "CDKN1C"
gene.exp <- gene_exp[gene,]


label <- if_else(gene.exp < median(as.numeric(gene.exp)), 0, 1)

group.low <- gene_exp[,label == 0]
group.high <- gene_exp[,label == 1]

group <- cbind(group.high,group.low)
group1 <- t(group)


Type <- c(rep('High',12),rep('Low',12))

riskTime <-as.data.frame(cbind(Type,group1))



group_list= riskTime$Type %>% factor(.,levels = c("High","Low"),ordered = F)


design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <-rownames(riskTime)
design <- as.data.frame(design)

dat1 <-gene_exp
dat1 <- as.data.frame(dat1)


max(dat1)
min(dat1)
dat1 <- log2(dat1+1)
dat1 <- dat1[,rownames(design)]

contrast.matrix <- makeContrasts(High-Low, levels = design)#Case比Control
fit <- lmFit(dat1, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)


allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)

allDiff <- allDiff[order(allDiff$logFC,decreasing=T),]

gene <- allDiff$logFC
names(gene) <-rownames(allDiff)


#2、GO####
#gmt

c5 <- read.gmt("F:/workspace/GSEA_gmt/c5.go.v7.4.symbols.gmt")
gsea <- GSEA(gene, TERM2GENE=c5, verbose=FALSE, pvalueCutoff = 0.05); head(gsea)
write.table(gsea,'gsea.CDKN1C.go.txt',sep = '\t',quote = F,row.names = F)
pdf('GSEA.CDKN1C.GO.pdf',width = 10,height = 10)
gseaplot2(gsea, geneSetID = 1:10,pvalue_table = FALSE,
          rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
dev.off()


#KEGG
c2 <- read.gmt("F:/workspace/GSEA_gmt/c2.cp.kegg.v7.4.symbols.gmt")
gsea <- GSEA(gene, TERM2GENE=c2, verbose=FALSE, pvalueCutoff = 0.05); head(gsea)
write.table(gsea,'gsea.CDKN1C.kegg.txt',sep = '\t',quote = F,row.names = F)
pdf('GSEA.CDKN1C.KEGG.pdf',width = 10,height = 10)
gseaplot2(gsea, geneSetID = 1:3,pvalue_table = FALSE,
          rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
dev.off()

