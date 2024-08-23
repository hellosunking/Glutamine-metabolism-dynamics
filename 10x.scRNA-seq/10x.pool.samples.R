#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

## usage: R --slave [--args sid sample.info] < basic.analysis.R

argv = commandArgs(T)
if( length(argv) != 2 ) {
	print( 'usage: R --slave --args <sample.info> <output.prefix> < pool.samples.R' )
	print( 'sample.info should contains 3 columns: sample.id /path/to/h5 /path/to/DoubletFinder.result')
	q()
}

suppressPackageStartupMessages({
	library(Seurat);library(plyr);library(dplyr);library(cowplot);
	library(ggplot2);library(gridExtra);library(ggthemes);library(patchwork);
	library(gridExtra);library(forcats);library(SingleR);
});
options( stringsAsFactors=F );

sid = as.character( argv[2] )
pdf( paste0(sid, '.pdf'), width=10, height=8 );

info = read.table( argv[1] )
## sid h5.path DoubletFinder.path

samples = info$V1
h5.path = info$V2
df.path = info$V3

valid   = list()
sc.list = list()
to.merge= c()

## load files and add sampleID to each barcode
for(i in 1:length(samples)) {
	print( samples[i] )
	counts = Read10X_h5( h5.path[i] )
	DoubletFinder = read.table( df.path[i] )
	pass.cells = subset( DoubletFinder, V2=="Singlet" )

	sce = CreateSeuratObject( counts, project=samples[i], min.features=1 )
	sce = subset( sce, cells=pass.cells$V1 )
	new.cell.id = paste0( samples[i], "-", gsub("-1", "", colnames(sce)) ) ## incase the cell id is like ACCGTGA-1
	sce = RenameCells(sce, new.names=new.cell.id )

	sc.list[[i]] = sce
	valid[[i]] = colnames( sce )

	if( i != 1 ) {
		to.merge = c( to.merge, sc.list[[i]] )
	}
}

## this is Seurat's merge
sce = merge(sc.list[[1]], y=to.merge);

## add a label of sclib
sce$sclib = sce$orig.ident;

## Perform cell cycle scoring, will plot this afterwards
#sce = CellCycleScoring(sce, g2m.features=cc.genes.updated.2019$g2m.genes,
#					   s.features=cc.genes.updated.2019$s.genes);

## regress out unwanted genes
all.gene = rownames( sce );
cc.genes = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes);
RP.genes = grep(rownames(sce), pattern="^RP[SL]", value=T, invert=F, ignore.case=T);
MT.genes = grep(rownames(sce), pattern="^MT-", value=T, invert=F, ignore.case=T);
keep.genes = all.gene[!(all.gene %in% c(cc.genes, RP.genes, MT.genes))];
sce = subset(sce, features=keep.genes);

## use Harmony package
library( harmony );
sce = NormalizeData(sce, normalization.method='LogNormalize', scale.factor=10000, verbose=F);
sce = FindVariableFeatures(sce, selection.method="vst", verbose=F);
sce = ScaleData( sce, verbose=F );
sce = RunPCA(sce, npcs=50, verbose=F);

sce = RunHarmony(sce, "sclib", plot_convergence=F, max.iter.harmony=20);
sce = FindNeighbors(sce, reduction="harmony", dims=1:30, verbose=F);
sce = FindClusters(sce, resolution=0.5, verbose=F);
#identity( sce );

sce = RunTSNE(sce, reduction="harmony", dims=1:30, verbose=F);
sce = RunUMAP(sce, reduction="harmony", dims=1:30, verbose=F);
sce = SetIdent(sce,value='RNA_snn_res.0.5');
sce$seurat_clusters = sce$RNA_snn_res.0.5;

DimPlot(sce, reduction="umap", group.by="sclib", pt.size=0.1, shuffle=T);
DimPlot(sce, reduction="umap", group.by="ident", label=T, pt.size=0.1);

## plot per sample
#for( i in 1:length(samples) ) {
	## note that DimPlot does not do ploting in loops
#	per.cell = DimPlot(sce, reduction="umap", cells=valid[[i]])+ggtitle(samples[i]);
#	grid.arrange( per.cell );
#}

## plot phase for cells, which could be used to evaluate the performance of integration
#DimPlot(sce, reduction="umap", group.by="Phase")

## plot known markers
## by default, UMAP results are shown
FeaturePlot(sce, reduction="umap", features=c("PTPRC", "EPCAM", "PECAM1", "COL3A1") )
#+scale_colour_gradientn( colours=rev(brewer.pal(n=11, name="Spectral")) )

## save project
save( sce, file=paste0(sid, ".rds") );

q()

########################################### OPTIONAL: add annotations ########################
cell.anno = read.table( "cell.anno", head=T )
#Donor	Cell	ClusterID	ClusterName
#BM1	BM1-2-TACAGTGCACCAACCG	1	CD34+.pre-B

anno = cell.anno$ClusterName
names(anno) = cell.anno$Cell
sce  = AddMetaData( sce, anno, col.name="CellType" )

donor = cell.anno$Donor
names(donor) = cell.anno$Cell
sce  = AddMetaData( sce, donor, col.name="Donor" )

#xx = data.frame(sce$seurat_clusters, sce$CellType)
#yy = paste0(xx$sce.CellType, "-", xx$sce.seurat_clusters)
#names(yy) = rownames(xx)
#sce = AddMetaData(sce, yy, col.name="NewCluster")

## GLS expression
#GLS.expr = sce[["RNA"]]
#GLS.mean.per.type=AverageExpression(sce, features="GLS", group.by = "CellType")

