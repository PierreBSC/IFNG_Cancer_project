
#Library Loading
library("BiocParallel")
library(DESeq2)
library(FactoMineR)
library(fifer) 
library(pheatmap) 

#Multicore registration :
N_core =6  # Be careful of not overloading your computer....
register(MulticoreParam(N_core))


##I) Loading/annotation/QC of the  data 
count_table=read.table("In_vitro_experiment/Data/featureCounts_count_matrix.txt",header = T,row.names = 1)
head(count_table)
dim(count_table)

##Creating annotation table

annotation_table=data.frame(Condition=c("Stimulated_1h","Stimulated_1h","Stimulated_1h",
                                        "Stimulated_24h","Stimulated_24h","Stimulated_24h",
                                        "Stimulated_6h","Stimulated_6h","Stimulated_6h",
                                        "Control","Control","Control"))

barplot(colSums(count_table)) ## -> deep sequencing + homogenous data

##II)Analysis using DESeq

dds <- DESeqDataSetFromMatrix(countData = count_table,
                              colData = annotation_table,
                              design = ~ Condition)
hist(log10(1+rowSums(count_table)),n=100)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds$condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)
plotDispEsts(dds)
plotMA(dds, ylim=c(-3,3),alpha = 0.01)
res <- results(dds)
resultsNames(dds)
resLFC_1h <- lfcShrink(dds, coef="Condition_Stimulated_1h_vs_Control", type="apeglm")
resLFC_6h <- lfcShrink(dds, coef="Condition_Stimulated_6h_vs_Control", type="apeglm")
resLFC_24h <- lfcShrink(dds, coef="Condition_Stimulated_24h_vs_Control", type="apeglm")

plotMA(resLFC_1h, ylim=c(-6,6),alpha = 0.01)
plotMA(resLFC_6h, ylim=c(-6,6),alpha = 0.01)
plotMA(resLFC_24h, ylim=c(-6,6),alpha = 0.01)

table_1h=as.matrix(resLFC_1h)
table_6h=as.matrix(resLFC_6h)
table_24h=as.matrix(resLFC_24h)

gene_info_bis=gene_info
rownames(gene_info_bis)=gene_info_bis$Ensembl_Gene_ID
gene_info_bis=gene_info_bis[rownames(resLFC_1h),]

rownames(table_1h)=gene_info_bis$Gene_name
rownames(table_6h)=gene_info_bis$Gene_name
rownames(table_24h)=gene_info_bis$Gene_name

plot(table_1h[,"log2FoldChange"],table_6h[,"log2FoldChange"])

transformed_expression <- rlog(dds, blind=FALSE)
transformed_expression=assay(transformed_expression)
rownames(transformed_expression)=gene_info_bis$Gene_name

gene_plot_function=function(gene) {
  par(las=1)
  u=data.frame(Gene=transformed_expression[gene,],
               Condition=factor(annotation_table$Condition,levels =c("Control","Stimulated_1h","Stimulated_6h","Stimulated_24h")))
  prism.plots(Gene~Condition,u,pch=21,bg="orange",cex=2,main=gene,col="black",cex.main=1.5,ylab="Gene expression (Log2)")
}
gene_plot_function("H2-K1")
gene_plot_function("Irf8")

###III)Export of the various results for later analysis

###A)Extract all genes differentially expressed at one time point at least once

gene_1h=rownames(table_1h)[abs(table_1h[,"log2FoldChange"])>1 & (table_1h[,"padj"])<0.01]
gene_1h=gene_1h[!is.na(gene_1h)]
gene_6h=rownames(table_6h)[abs(table_6h[,"log2FoldChange"])>1 & (table_6h[,"padj"])<0.01]
gene_6h=gene_6h[!is.na(gene_6h)]
gene_24h=rownames(table_24h)[(abs(table_24h[,"log2FoldChange"])>1 & (table_24h[,"padj"])<0.01)]
gene_24h=gene_24h[!is.na(gene_24h)]


par(las=1,mar=c(6,6,4,4))
barplot(c(length(gene_1h),length(gene_6h),length(gene_24h)),
        xlab="Stimulation time",ylab="Number of DE genes",cex.lab=1.4,
        names=c("1h","6h","24h"),cex.axis=1.4,main="Number of DE genes across conditions")

gene_total=unique(c(gene_1h,gene_6h,gene_24h)) ## 268 genes in total

MHC_genes=gene_total[grepl("H2",gene_total)]
pheatmap(transformed_expression[MHC_genes,],scale = "row")

write.table(gene_1h,file = "In_vitro_experiment/DE_Gene_list/1h_stim.txt",sep="\t",quote=F,row.names = F)
write.table(gene_6h,file = "In_vitro_experiment/DE_Gene_list/6h_stim.txt",sep="\t",quote=F,row.names = F)
write.table(gene_24h,file = "In_vitro_experiment/DE_Gene_list/24h_stim.txt",sep="\t",quote=F,row.names = F)

#B)Export of the plots



pdf("In_vitro_experiment/Plot/Gene_plot.pdf")
par(las=1)
for (k in gene_total) {
  gene_plot_function(k)
  cat(paste(k,"\n"))
}
dev.off()


pdf("In_vitro_experiment/Plot/Heatmap_in_vitro_1.pdf",width = 12,height = 30)
pheatmap(transformed_expression[gene_total,c("unstimulated.1","unstimulated.2","unstimulated.3",
                                             "X1h.1","X1h.2","X1h.3",
                                             "X6h.1","X6h.2","X6h.3",
                                             "X24h.1","X24h.2","X24h.3")],gaps_col = c(3,6,9),scale = "row",cluster_cols = F)
dev.off()



pdf("In_vitro_experiment/Plot/Heatmap_in_vitro_2.pdf",width = 8,height = 8)
pheatmap(transformed_expression[gene_total,c("unstimulated.1","unstimulated.2","unstimulated.3",
                                             "X1h.1","X1h.2","X1h.3",
                                             "X6h.1","X6h.2","X6h.3",
                                             "X24h.1","X24h.2","X24h.3")],
         gaps_col = c(3,6,9),scale = "row",cluster_cols = F,show_colnames = F,show_rownames = F)
dev.off()

#### Optional analysis : nice Venn plot
require(gplots) 

v=get.venn.partitions(list(gene_1h,gene_6h,gene_24h))
pdf("In_vitro_experiment/Plot/Venn_diagram.pdf")
venn(list("1h Stimulation"=gene_1h,
          "6h Stimulation"=gene_6h,
          "24h Stimulation"=gene_24h),small = 1)
dev.off()

####
####Barplot Motif enrichment
Enrichment_6h = read.delim("In_vitro_experiment/Analysis_output/Enrichment_6h.tsv",skip=11)
pdf("In_vitro_experiment/Plot/Enrichment_barplot_6h.pdf")
par(las=1)
barplot(Enrichment_6h$NES[1:10],ylim=c(0,20),ylab="NES score",cex.lab=1.4,cex.axis=1.4,
        xlab="Ranked motif",cex.lab=1.3,col=c(rep("orange",3),rep("grey",7)))
dev.off()

Enrichment_24h = read.delim("In_vitro_experiment/Analysis_output/Enrichment_24h.tsv",skip=11)

pdf("In_vitro_experiment/Plot/Enrichment_barplot_24h.pdf")
par(las=1)
barplot(Enrichment_24h$NES[1:10],ylim=c(0,20),ylab="NES score",cex.lab=1.4,cex.axis=1.4,
        xlab="Ranked motif",cex.lab=1.3,col=c(rep("orange",3),rep("grey",7)))
dev.off()
