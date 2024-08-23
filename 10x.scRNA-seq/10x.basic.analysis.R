#!/usr/bin/R
#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :

argv = commandArgs(T)
if( length(argv) != 2 ) {
	print( 'usage: R --slave --args <filtered.h5> <output.prefix> < basic.analysis.R' )
	q()
}

sid=as.character( argv[2] )

options( stringsAsFactors=F )
suppressPackageStartupMessages({
	library(Seurat);library(plyr);library(dplyr);library(cowplot);
	library(ggplot2);library(gridExtra);library(ggthemes);
	library(patchwork);library(forcats);library(RColorBrewer);
	library(DoubletFinder);library(SingleR);
})

outfileName=paste0(sid, '.basic.pdf')
pdf( outfileName, width=10, height=8 )

## load data and build Seurat object
cr = Read10X_h5( argv[1] )
sc.raw = CreateSeuratObject(counts=cr, project=sid, min.cells=3, min.features=100)

## Quality control
sc.raw[["percent.mito"]] = PercentageFeatureSet(sc.raw, pattern="^MT-")
sc.raw[["percent.ribo"]] = PercentageFeatureSet(sc.raw, pattern="^RP[SL]")
sc.raw[["percent.hbab"]] = PercentageFeatureSet(sc.raw, pattern="^HB[AB][12]?$")

p1=VlnPlot(sc.raw, features='nCount_RNA', y.max=50000,
		   pt.size=0.01)+NoLegend()
p2=VlnPlot(sc.raw, features='nFeature_RNA', y.max=10000,
		   pt.size=0.01)+NoLegend()
p3=VlnPlot(sc.raw, features='percent.mito', y.max=100,
		   pt.size=0.01)+NoLegend()
p4=VlnPlot(sc.raw, features='percent.ribo', y.max=100,
		   pt.size=0.01)+NoLegend()
px=VlnPlot(sc.raw, features='percent.hbab', y.max=100,
		   pt.size=0.01)+NoLegend()
grid.arrange(p1, p2, p3, p4, px, nrow=1,
			 top=paste0('Before QC, cells=', dim(sc.raw)[2]))

p5=FeatureScatter(sc.raw, feature1="nCount_RNA", feature2="percent.mito")
p6=FeatureScatter(sc.raw, feature1="nCount_RNA", feature2="nFeature_RNA")
p7=FeatureScatter(sc.raw, feature1="nCount_RNA", feature2="percent.ribo")
p8=FeatureScatter(sc.raw, feature1="percent.mito", feature2="percent.ribo")
grid.arrange(p5, p6, p7, p8, nrow=2)

## filter outliners and cells with abnormally high chrM (i.e., dead cells)
sce = subset(sc.raw, percent.mito < 25 & percent.hbab < 5)

## Perform cell cycle scoring, will plot this afterwards
sce = CellCycleScoring(sce, g2m.features=cc.genes.updated.2019$g2m.genes,
					   s.features=cc.genes.updated.2019$s.genes)

## regress out unwanted genes
all.gene = rownames( sce )
cc.genes = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
RP=grep(rownames(sce), pattern="^RP[SL]",
		value=T, invert=F, ignore.case=T)
MT=grep(rownames(sce), pattern="^MT-",
		value=T, invert=F, ignore.case=T)
keep.genes = all.gene[!(all.gene %in% c(cc.genes, RP, MT))]
sce = subset(sce, features=keep.genes)

## normalize & scale data, run PCA
sce = NormalizeData(sce, normalization.method='LogNormalize', scale.factor=10000, verbose=F)
sce = FindVariableFeatures(sce, selection.method='vst', nfeatures=2000, verbose=F)
#LabelPoints(plot=VariableFeaturePlot(sce), points=head(VariableFeatures(sce), 10), repel=T);
sce = ScaleData(sce)
sce = RunPCA(sce, features=VariableFeatures(sce), verbose=F)

## run t-SNE
sce = FindNeighbors(sce, dims=1:30)
sce = FindClusters(sce, resolution=c(0.4,0.5,0.6), verbose=F)
sce = RunTSNE(sce, dims=1:30)
sce = RunUMAP(sce, dims=1:30, verbose=F)
## I will fix the resolution to 0.5 here
## because this is a basic analysis and I will pool all the cells
## and re-analyze the data later
sce = SetIdent(sce,value='RNA_snn_res.0.5')
sce$seurat_clusters=sce$RNA_snn_res.0.5

#DimPlot(sce,reduction = "tsne", label=T,split.by ='orig.ident');
#LabelClusters(DimPlot(sce, reduction = "tsne"),id = 'ident')
TSNEPlot(sce, pt.size=1, label=T, group.by='RNA_snn_res.0.5')+
ggtitle('res=0.5, before doublet removal')
UMAPPlot(sce, pt.size=1, label=T, group.by='RNA_snn_res.0.5')+
ggtitle('res=0.5, before doublet removal')

## detect doublets using DoubletFinder
## https://www.jianshu.com/p/6770c6a05287
## find the optimal parameter in DoubletFinder
## if log-normalization is used, set sct=F (default)
## otherwise set it to T (is SCT-normalization is used)
sweep.res.list = paramSweep_v3(sce, PCs=1:30, sct=F)
sweep.stats = summarizeSweep(sweep.res.list, GT=F)
bcmvn = find.pK(sweep.stats)
pK_good = bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

## in 10x, doublet rate is related to the No. of recovered cells,
## increases ~0.8% when recovering additional 1000 cells
## Note that here I use the sc.raw as the recovered cells
DoubletRate = dim(sc.raw)[2] * 8e-6
homotypic.prop = modelHomotypic(sce$seurat_clusters)
nExp_poi = round(DoubletRate*ncol(sce))
nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))
sce = doubletFinder_v3(sce, PCs=1:30, pN=0.25, pK=pK_good, nExp=nExp_poi.adj, reuse.pANN=F, sct=F)

## Plot results
## note that DoubletFinder will add 1 column named "DF.classifications_XXXX"
## which is parameter-specific and need to be renamed
index = grep("DF.classifications", colnames(sce@meta.data))
colnames(sce@meta.data)[index] = "DFpred"
DimPlot(sce, reduction="umap", group.by="DFpred",
		order=c("Singlet", "Doublet"), cols=c("red", "grey"))+
ggtitle( paste0("DoubletFinder result, pK=", pK_good) )
write.table(sce$DFpred,	file=paste0(sid, ".DoubletFinder.result"),
			quote=F, sep='\t', row.names=T, col.names=F)

## show No. of cells per cluster before and after doublet removal
count.before.DF = as.data.frame(table(sce@meta.data$RNA_snn_res.0.5))

## filter out doublets
sce = subset( sce, DFpred=="Singlet" );
count.after.DF = as.data.frame(table(sce@meta.data$RNA_snn_res.0.5))

cluster.num = length( count.before.DF$Var1 )
count.matrix = data.frame( cbind(
	rep(count.before.DF$Var1, 2),
	c(rep("0.Before", cluster.num), rep("1.After", cluster.num)),
	c(count.before.DF$Freq, count.after.DF$Freq)
))
colnames(count.matrix)=c("ClusterID", "Condition", "CellNum")
write.table( count.matrix, file=paste0(sid, ".cell.counts"), sep="\t", quote=F)
count.matrix = read.table( paste0(sid, ".cell.counts") )
ggplot(count.matrix, aes(x=ClusterID, y=CellNum))+
geom_bar(stat='identity', aes(fill=Condition),position=position_dodge(0.9))+
labs(title="Cell counts before/after DoubletFinder")+
xlab("Cluster ID")+ylab("No. of cells")+
theme_bw()+theme(axis.text.x=element_text(size=6))+
scale_fill_manual(values=c("orange","steelblue"))+
guides(fill=guide_legend(title=NULL))

## Plot the cluster results in various resolutions
#TSNEPlot(sce, pt.size=.5, label=T, group.by='RNA_snn_res.0.3')+
#ggtitle('res=0.3, after doublet removal');
#UMAPPlot(sce, pt.size=.5, label=T, group.by='RNA_snn_res.0.3')+
#ggtitle('res=0.3, after doublet removal');

TSNEPlot(sce, pt.size=.5, label=T, group.by='RNA_snn_res.0.5')+
ggtitle('res=0.5, after doublet removal');
UMAPPlot(sce, pt.size=.5, label=T, group.by='RNA_snn_res.0.5')+
ggtitle('res=0.5, after doublet removal');

## plot statistics in each cluster
p1=VlnPlot(sce, features='nCount_RNA',   y.max=30000, pt.size=0)+NoLegend()
p2=VlnPlot(sce, features='nFeature_RNA', y.max=10000, pt.size=0)+NoLegend()
p3=VlnPlot(sce, features='percent.mito', y.max=25, pt.size=0)+NoLegend()
p4=VlnPlot(sce, features='percent.ribo', y.max=100,pt.size=0)+NoLegend()
grid.arrange(p1, p2, p3, p4, nrow=2, top=paste0('After QC (MT%<25) & doublet removal, cells=', dim(sce)[2]))

## show the UMI/feature/MT numbers and phase per cluster
pa = FeaturePlot(sce, reduction="umap", features='nCount_RNA')
pb = FeaturePlot(sce, reduction="umap", features='nFeature_RNA')
pc = FeaturePlot(sce, reduction="umap", features='percent.mito')
pd = DimPlot(sce, reduction="umap", group.by="Phase")
grid.arrange(pa, pb, pc, pd, nrow=2)

## find cluster markers
cm = FindAllMarkers(sce, logfc.threshold=0.2, min.pct=0.25, only.pos=T)
cm = cm[with(cm,order(cluster, p_val_adj, -avg_log2FC)),]

## plot top markers
source( "bubblePlot.R" )
top10 = cm %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
ngene = length(top10$gene)
part1 = round(ngene * 0.5)
part2 = ngene - round(ngene * 0.5)
# wt could be avg_logFC or -p_val_adj
DoHeatmap(sce, features=head(top10$gene, n=part1)) + NoLegend()
qqq(sce, features=unique(top10$gene)[1 : part1], cols=c('#08519c','#b30000'),
	group.by='seurat_clusters', dot.scale=4., scale=F,
	col.min=0, col.max=3, scale.min=0,scale.max=100)+
xlab('Cluster ID') + ggtitle('Top10 marker genes per cluster, part1')+
theme(axis.text.y=element_text(size=8))

DoHeatmap(sce, features=tail(top10$gene, n=part2)) + NoLegend()
qqq(sce, features=unique(top10$gene)[part1+1 : ngene],
	cols=c('#08519c','#b30000'), group.by='seurat_clusters', dot.scale=4., scale=F,
	col.min=0, col.max=3, scale.min=0,scale.max=100)+
xlab('Cluster ID') + ggtitle('Top10 marker genes per cluster, part2')+
theme(axis.text.y=element_text(size = 8))

write.table(cm, file=paste0(sid, ".cluster.markers"), quote=F, sep='\t', row.names=F, col.names=T)

## plot known markers
## by default, UMAP results are shown
FeaturePlot(sce, reduction="umap", features=c("PTPRC", "EPCAM", "PECAM1", "COL3A1") ) +
scale_colour_gradientn( colours=rev(brewer.pal(n=11, name="Spectral")) )

## automated annotation using SingleR
load("ref_Human_all.RData")
sce_for_SingleR = GetAssayData(sce, slot="data");
pred = SingleR(test=sce_for_SingleR, ref=ref_Human_all,
			   labels=ref_Human_all$label.main, method="cluster",
			   clusters=sce@meta.data$seurat_clusters,
			   assay.type.test="logcounts", assay.type.ref="logcounts");
celltype = data.frame(ClusterID=rownames(pred), celltype=pred$labels, stringsAsFactors=F);
sce@meta.data$SingleR = celltype[match(sce@meta.data$seurat_clusters, celltype$ClusterID), 'celltype'];
TSNEPlot(sce, group.by="SingleR", label=T, label.size=5);
UMAPPlot(sce, group.by="SingleR", label=T, label.size=5);
write.table(sce$SingleR, file=paste0(sid, ".SingleR.result"), quote=F, sep='\t', row.names=T, col.names=F);

## example manual annotation
#new.cluster.ids = c("Plasma", "T", "Goblet", "NK", "B", "Fibroblast", "Monocyte", "Erythroblast", "Myofibroblast", "Endothelial" )
#names(new.cluster.ids) = levels(sce);
#sce = RenameIdents(sce, new.cluster.ids);
#DimPlot(sce, reduction="umap", label=T, pt.size=0.8, label.size=8) + NoLegend();

## save project
save(sce, file=paste0(sid, ".rds"))

