library(fifer)
library(uwot)
library(pagoda2)
library(Matrix)
library(igraph)
library(GEOquery)
library(RColorBrewer)
library(Desktop/HNSCC_Ronan/Figures)

#Analysis of IFNg signature in infiltrating immune cells of human melanoma 
#Original paper : "Dysfunctional CD8+ T cells form a proliferative, dynamically regulated compartment within human melanoma" by Li et al. 
#Single-cell technology used : MARS-seq 



##Two auxiliary functions
get_otsus_threshold = function(x) {
  y = x[order(x)]
  intra_variance = vapply(X = 1:length(y),FUN = function(i) {
    s = i * var(y[1:i]) + (length(y)-i) * var(y[(i+1):length(y)])
    return(s)},FUN.VALUE = 1)
  identified_value = which.min(intra_variance)
  
  identified_value = y[identified_value]
  return(identified_value)
}


color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.999,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}

#I)Loading and aggregation of the data

#A)Expression data


##All UMI tabs were manually downloaded from the GEO website (GSE123139)
##Link used : https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE123139&format=file

l_files =list.files("IFNG_project/Melanoma_MARSseq/UMI_tab/",full.names = T)

data_raw = read.table(l_files[1],header = T,row.names = 1,sep = "\t")

for (k in l_files[59:length(l_files)]) {
  u = read.table(k,header = T,row.names = 1,sep = "\t")
  u = u[rownames(data_raw),]
  data_raw = cbind(data_raw,u)
  print(k)
}
data_raw = as(as.matrix(data_raw),"dgCMatrix")

#B)Annotation data

annotation_plates = getGEO("GSE123139")
annotation_plates= annotation_plates$GSE123139_series_matrix.txt.gz
annotation_plates = pData(annotation_plates)

Gating_strategy = as.character(annotation_plates$`facs gate:ch1`)
Gating_strategy = rep(Gating_strategy,each=384)
names(Gating_strategy) = colnames(data_raw)

Patient = as.character(annotation_plates$`patient id:ch1`)
Patient = rep(Patient,each=384) #Each plate : 384 cells...
names(Patient) = colnames(data_raw)

Tissue_type = as.character(annotation_plates$`sample source:ch1`)
Tissue_type = rep(Tissue_type,each=384)
names(Tissue_type) = colnames(data_raw)

#C)Filtering genes and cell

lib_size = Matrix::colSums(data_raw)
hist(log10(lib_size),100)
abline(v=log10(350))

gene_size = Matrix::rowSums(data_raw)
hist(log10(gene_size),100)
abline(v=log10(300))

n_cell_patients = table(Patient)
n_cell_patients = n_cell_patients[order(n_cell_patients,decreasing = T)]
selected_patients = names(n_cell_patients[1:9])
sum(Patient%in%selected_patients)

data_count = data_raw[gene_size>300,lib_size>350 & Patient%in%selected_patients ]

Patient_count = Patient[colnames(data_count)]
Gating_count = Gating_strategy[colnames(data_count)]

#II)Data analysis per se

r <- Pagoda2$new(data_count,log.scale=F)
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine',n.cores = 10)
r$getKnnClusters(method=multilevel.community,type='PCA')

UMAP_plot = umap(r$reductions$PCA,verbose=T,metric= "cosine",n_neighbors = 30,spread = 3)

#Distribution of clusters across UMAP space
plot(UMAP_plot,pch=21,bg=string.to.colors(r$clusters$PCA$community),xaxt="n",yaxt="n",xlab="UMAP 1",ylab="UMAP 2",bty='n')
#Gene expression of one gene (Granzyme A gene) across UMAP space
plot(UMAP_plot,pch=16,col=alpha(color_convertion(r$counts[,"GZMA"]),alpha = 0.1))

r$getDifferentialGenes(type = "PCA",clusterType = "community",upregulated.only = T,verbose = T)

#III)Studying IFNg production

#A)Which cells produce IFNG ?

boxplot(r$counts[,"CD8A"]~factor(r$clusters$PCA$community,levels = order_cluster))

##Computing production of IFNg by each cell type from the CD3+ compartment
Total_IFNG_prod = aggregate(as.matrix(r$misc$rawCounts[Gating_count=="CD3","IFNG"]),by=list(r$clusters$PCA$community[Gating_count=="CD3"]),FUN=sum)
Total_IFNG_prod= Total_IFNG_prod$V1
names(Total_IFNG_prod) = 1:length(Total_IFNG_prod)
Total_IFNG_prod = Total_IFNG_prod/sum(Total_IFNG_prod)
Total_IFNG_prod = Total_IFNG_prod[order(Total_IFNG_prod,decreasing = F)]

col_list =brewer.pal(7,"Paired")

pdf("IFNG_project/Melanoma_MARSseq/Figures/IFNG_distribution_CD3.pdf")
par(las=1,mar=c(5,7,4,22))
barplot(matrix(Total_IFNG_prod,ncol=1)*100,horiz = F,col=col_list,
        ylab="Total IFNG UMI production",cex.lab=1.3,cex.axis = 1.3)
dev.off()

##3 clusters : classic CD8+ T , super-activated CD8+ T cells and dividing CD8+ produces most of the IFNg

#B)Production of IFNg across patients 

mean_IFNg_CD3_per_patient = aggregate(log2(1+r$counts[CD3_cells,"IFNG"]),
                                      by=list(Patient_count[Gating_count=="CD3"]),FUN=mean)
mean_IFNg_CD3_per_patient=mean_IFNg_CD3_per_patient$x


par(las=1)
barplot(mean_IFNg_CD3_per_patient,xlab="Patient",cex.lab=1.3,
        ylab="Mean IFNG expression (Log(1+TPM))",main="Mean IFNG gene expression by CD3+ cells")

#C)Proportion of IFNG positive CD8+ T cells

CD8_cells =r$clusters$PCA$community%in%c(14,17,19)
hist(log10(1+r$counts[CD8_cells,"IFNG"]),100)
CD8_cells_IFNg_pos = CD8_cells & r$counts[,"IFNG"]>0

CD3_proportion_CD8_IFNG_pos = table(Patient_count[Gating_count=="CD3"],CD8_cells_IFNg_pos[Gating_count=="CD3"])
CD3_proportion_CD8_IFNG_pos = CD3_proportion_CD8_IFNG_pos/rowSums(CD3_proportion_CD8_IFNG_pos)
barplot(t(CD3_proportion_CD8_IFNG_pos),col=c("grey","orange"))


#IV)Studying IFNG effect 

#A)Three cell populations as potential IFNg target : Macrophages, Neutrophils and B cells

#Which cells express IFNg receptor (IFNGR1 gene)

boxplot(log10(1+r$counts[,"IFNGR1"])~r$clusters$PCA$community,main="IFNGR1 expression across clusters")

#B)Studying Neutro heterogeneity 

#The cluster number of each cell type need to be adjusted... 
#Here neutrophils correspond to cluster 2 (check for S100A8/9 positive cluster)

r_neutro = Pagoda2$new(data_count[,r$clusters$PCA$community==2],log.scale=T)
r_neutro$adjustVariance(plot=T,gam.k=10)
r_neutro$calculatePcaReduction(nPcs=50,n.odgenes=5e2)
r_neutro$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine',n.cores = 10)
r_neutro$getKnnClusters(method=multilevel.community,type='PCA')

UMAP_neutro = umap(r_neutro$reductions$PCA,verbose=T,metric= "cosine",n_neighbors = 30,spread = 4)
plot(UMAP_neutro,pch=21,bg=string.to.colors(Patient_count[r$clusters$PCA$community==2]),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
plot(UMAP_neutro,pch=21,bg=string.to.colors(r_neutro$clusters$PCA$community))
plot(UMAP_neutro,pch=21,bg=color_convertion(r_neutro$counts[,"GBP1"]))

#C)Studying Macro heterogeneity 

r_macro = Pagoda2$new(data_count[,r$clusters$PCA$community==27],log.scale=T)
r_macro$adjustVariance(plot=T,gam.k=10)
r_macro$calculatePcaReduction(nPcs=50,n.odgenes=5e2)
r_macro$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine',n.cores = 10)
r_macro$getKnnClusters(method=multilevel.community,type='PCA')

UMAP_macro = umap(r_macro$reductions$PCA,verbose=T,metric= "cosine",n_neighbors = 30,spread = 4)
plot(UMAP_macro,pch=21,bg=string.to.colors(Patient_count[r$clusters$PCA$community==27]),
     xaxt="n",yaxt="n",xlab= "UMAP 1",ylab="UMAP 2",bty="n")
plot(UMAP_macro,pch=21,bg=string.to.colors(r_macro$clusters$PCA$community))
plot(UMAP_macro,pch=21,bg=color_convertion(r_macro$counts[,"GBP1"]))

r_macro$getDifferentialGenes(type = "PCA",clusterType = "community",upregulated.only = T,verbose = T)

#D)Studying B cells  diversity

r_B_cells = Pagoda2$new(data_count[,r$clusters$PCA$community==9],log.scale=T)
r_B_cells$adjustVariance(plot=T,gam.k=2)
r_B_cells$calculatePcaReduction(nPcs=50,n.odgenes=5e2)
r_B_cells$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine',n.cores = 10)
r_B_cells$getKnnClusters(method=multilevel.community,type='PCA')

UMAP_B_cells = umap(r_B_cells$reductions$PCA,verbose=T,metric= "cosine",n_neighbors = 30,spread = 4)
plot(UMAP_B_cells,pch=21,bg=string.to.colors(Patient_count[r$clusters$PCA$community==9]),
     xaxt="n",yaxt="n",xlab= "UMAP 1",ylab="UMAP 2",bty="n")
plot(UMAP_macro,pch=21,bg=string.to.colors(r_macro$clusters$PCA$community))
plot(UMAP_B_cells,pch=21,bg=color_convertion(r_B_cells$counts[,"GBP1"]))

#V)Computing IFNg score for each cell

#A different IFNG score need to be compute for each cell type 

#A)For macrophages 

#Performing de-novo over-dispersion test  


#Clustering genes based on correlation

cor_macro_cluster = cor(as.matrix(r_macro$counts[,r_macro$getOdGenes(200)]),method = "pearson")
Gene_clustering_macro=hclust(dist(cor_macro_cluster),method = "ward")
Gene_clustering_macro=cutree(Gene_clustering_macro,k = 15) #15 clusters studied

gene_env=c()
for (k in 1:length(unique(Gene_clustering_macro))) {
  gene_env[[k]]=names(which(Gene_clustering_macro==k))
}
names(gene_env)=1:15
gene_env_2 <- list2env(gene_env) # convert to an environment

#Over-dispersion test 
r_macro$testPathwayOverdispersion(setenv = gene_env_2,type = "counts",verbose = T,
                                   plot = F,max.pathway.size = 100,min.pathway.size = 1,recalculate.pca=T)
pathway_info=r_macro$misc$pathwayODInfo
pathway_info=pathway_info[order(pathway_info$cz,decreasing=T),]
View(pathway_info)
dispersed_genes=c()

for (i in (substr(rownames(pathway_info),7,10))) {
  u=r_macro$misc$pwpca[[i]]$xp$rotation
  u=u[unique(rownames(u)),]
  u=(u[order(u,decreasing = F)])
  dispersed_genes[[i]]=u
  print(i)
}

#Which gene set corresponds to IFNg signature ?
#To identify the good one we used the references signatures provided in the Hallmark gene sets from the MSigDB 
#IFN alpha signature : software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=HALLMARK_INTERFERON_ALPHA_RESPONSE&fileType=txt
#IFN gamma signature : software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=HALLMARK_INTERFERON_GAMMA_RESPONSE&fileType=txt

IFNG_gene_set = read.table("IFNG_project/Melanoma_MARSseq/Geneset_IFNG.txt",skip = 2)
IFNG_gene_set = intersect(IFNG_gene_set$V1,colnames(r_macro$counts))

IFNB_gene_set = read.table("IFNG_project/Melanoma_MARSseq/Geneset_IFNB.txt",skip = 2)
IFNB_gene_set = intersect(IFNB_gene_set$V1,colnames(r_macro$counts))

#Looking for enrichment in the de-novo gene sets of the two different hallmark signatures 

IFN_GS_enrichment = c()
for (k in names(dispersed_genes)) {
  genes = names(dispersed_genes[[k]])
  e_IFNG = sum(genes%in%c(IFNG_gene_set))/length(genes)
  e_IFNB = sum(genes%in%c(IFNB_gene_set))/length(genes)
  p_IFNG = binom.test(x = sum(genes%in%c(IFNG_gene_set)),
                      n = length(genes),p = sum(r_macro$getOdGenes(200)%in%IFNG_gene_set)/length(r_macro$getOdGenes(200)) )
  p_IFNB = binom.test(x = sum(genes%in%c(IFNB_gene_set)),
                      n = length(genes),p = sum(r_macro$getOdGenes(200)%in%IFNB_gene_set)/length(r_macro$getOdGenes(200)) )
  IFN_GS_enrichment = rbind(IFN_GS_enrichment,c(e_IFNG,e_IFNB,p_IFNG$p.value,p_IFNB$p.value))
}
rownames(IFN_GS_enrichment) = names(dispersed_genes)
colnames(IFN_GS_enrichment) = c("Enrichment_IFNG","Enrichment_IFNB","P_IFNG","P_IFNB")
IFN_GS_enrichment = as.data.frame(IFN_GS_enrichment)
IFN_GS_enrichment$P_IFNG = p.adjust(IFN_GS_enrichment$P_IFNG,method = "BH")
IFN_GS_enrichment$P_IFNB = p.adjust(IFN_GS_enrichment$P_IFNB,method = "BH")


pdf("IFNG_project/Melanoma_MARSseq/Figures/Enrichment_signatures.pdf",useDingbats = F)
par(las=1)
barplot(rbind(-log10(IFN_GS_enrichment$P_IFNG),-log10(IFN_GS_enrichment$P_IFNB)),cex.lab=1.3,cex.axis=1.3,
        beside = T,names.arg = 1:15,xlab="Pathway",ylab="-Log10(P-value)",col=c("orange","grey"))
legend("topright",legend = c("IFNG signature","IFNB signature"),fill=c("orange","grey"),bty='n')
abline(h=2,lty=2,lwd=2)
dev.off()

#One signature is only enriched in IFNg genes but not in IFN-alpha -> the good one

pdf("IFNG_project/Melanoma_MARSseq/Figures/Contribution_gene_macrophages_score.pdf",useDingbats = F)
par(las=1,mar=c(6,6,4,10))
barplot(dispersed_genes$`6`,horiz = T,xlab="Contribution to Macrophage IFNG signature",cex.lab=1.3,col="dodgerblue3")
dev.off()

par(las=1,mar = c(4,6,4,2))
IFN_score_macro = as.numeric(r_macro$misc$pwpca$`6`$xp$scores)
IFN_score_macro = IFN_score_macro-min(IFN_score_macro)
IFN_score_macro = data.frame(Score = IFN_score_macro , 
                             Patient = Patient_count[r$clusters$PCA$community==27])

Mean_IFN_score_patient_macro = aggregate(IFN_score_macro$Score,by=list(IFN_score_macro$Patient),mean)

n_patient = table(IFN_score_macro$Patient)
n_patient = n_patient[order(n_patient,decreasing = T)]
top_patients = names(n_patient[1:5]) #Taking only the 5 patients with the highest number of Macrophages

pdf("IFNG_project/Melanoma_MARSseq/Figures/Distribution_IFNG_score_macrophages.pdf",width = 7,height = 5,useDingbats = F)
p= ggplot(IFN_score_macro[IFN_score_macro$Patient%in%top_patients,], aes(x = Score, fill = Patient)) + geom_density(alpha = 0.4) + labs(x = "New x label")
p = p + labs(x = "IFNG score", y = "Density") +theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot(p)
dev.off()


#B)For neutrophils

cor_neutro_cluster = cor(as.matrix(r_neutro$counts[,r_neutro$getOdGenes()]),method = "pearson")
Gene_clustering_neutro=hclust(dist(cor_neutro_cluster),method = "ward")
Gene_clustering_neutro=cutree(Gene_clustering_neutro,k = 15)


gene_env=c()
for (k in 1:length(unique(Gene_clustering_neutro))) {
  gene_env[[k]]=names(which(Gene_clustering_neutro==k))
}
names(gene_env)=1:15
gene_env_2 <- list2env(gene_env) # convert to an environment

r_neutro$testPathwayOverdispersion(setenv = gene_env_2,type = "counts",verbose = T,
                                   plot = F,max.pathway.size = 100,min.pathway.size = 1,recalculate.pca=T)
pathway_info=r_neutro$misc$pathwayODInfo
pathway_info=pathway_info[order(pathway_info$cz,decreasing=T),]
View(pathway_info)
dispersed_genes=c()

for (i in (substr(rownames(pathway_info),7,10))) {
  u=r_neutro$misc$pwpca[[i]]$xp$rotation
  u=u[unique(rownames(u)),]
  u=(u[order(u,decreasing = F)])
  dispersed_genes[[i]]=u
  print(i)
}

par(las=1,mar = c(4,6,4,2))
IFN_score_neutro = as.numeric(r_neutro$misc$pwpca$`9`$xp$scores) #9th score : corresponds to IFNg score 
IFN_score_neutro = IFN_score_neutro-min(IFN_score_neutro) #Having only positive IFNG score
plot(UMAP_neutro,pch=21,bg=color_convertion(IFN_score_neutro))

IFN_score_neutro = data.frame(Score = IFN_score_neutro , Patient = Patient_count[r$clusters$PCA$community==2])

Mean_IFN_score_patient_neutro = aggregate(IFN_score_neutro$Score,by=list(IFN_score_neutro$Patient),mean)


#Which genes are contributing to IFNg signature in Neutrophils

pdf("IFNG_project/Melanoma_MARSseq/Figures/Contribution_gene_neutrophil_score.pdf")
par(las=1,mar=c(6,6,4,10))
barplot(dispersed_genes$`9`,horiz = T,xlab="Contribution to Neutrophils IFNG signature",cex.lab=1.3,col="dodgerblue3")
dev.off()


#Distribution of IFNg score across patients in neutrophils
n_patient = table(IFN_score_neutro$Patient)
n_patient = n_patient[order(n_patient,decreasing = T)]
top_patients = names(n_patient[1:5]) #Taking only the 5 patients with the highest number of Macrophages


pdf("IFNG_project/Melanoma_MARSseq/Figures/Distribution_IFNG_score_neutrophils.pdf",useDingbats = F,width = 7,height = 5)
p= ggplot(IFN_score_neutro[IFN_score_neutro$Patient%in%top_patients,], aes(x = Score, fill = Patient)) + geom_density(alpha = 0.4) + labs(x = "New x label")
p = p + labs(x = "IFNG score", y = "Density") +theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot(p)
dev.off()



#VI)Linking mean IFNg production and mean IFNg score across patients 

#Taking only patients with at least 15 target cells 


##For Macrophages

n_patient = table(IFN_score_macro$Patient)
n_patient = n_patient[Mean_IFN_score_patient_macro$Group.1]

pdf("IFNG_project/Melanoma_MARSseq/Figures/Correlation_IFN_response_macrophages.pdf",useDingbats = F)
par(las=1)
plot(mean_IFNg_CD3_per_patient[n_patient>15],Mean_IFN_score_patient_macro$x[n_patient>15],pch=21,bg="orange",cex=2,
     xlab = "Mean IFNG expression",ylab = "Mean IFN score in Macrophages",cex.axis=1.2,cex.lab=1.4,xaxs="i",yaxs="i",xlim=c(0,3),ylim=c(0,5))
 cor(mean_IFNg_CD3_per_patient[n_patient>15],Mean_IFN_score_patient_macro$x[n_patient>15])
m = lm(Mean_IFN_score_patient_macro$x[n_patient>15]~mean_IFNg_CD3_per_patient[n_patient>15])
abline(coef(m),lwd=2,lty=2)
summary(m)


##For Neutrophils

n_patient = table(IFN_score_neutro$Patient)
n_patient = n_patient[Mean_IFN_score_patient_neutro$Group.1]

par(las=1)
pdf("IFNG_project/Melanoma_MARSseq/Figures/Correlation_IFN_response_neutrophils.pdf",useDingbats = F)
plot(mean_IFNg_CD3_per_patient[n_patient>15],Mean_IFN_score_patient_neutro$x[n_patient>15],pch=21,bg="orange",cex=2,
     xlab = "Mean IFNG expression",ylab = "Mean IFN score in Neutrophils",cex.axis=1.2,cex.lab=1.4,xaxs="i",yaxs="i",xlim=c(0,3),ylim=c(0,4))
cor(mean_IFNg_CD3_per_patient[n_patient>15],Mean_IFN_score_patient_neutro$x[n_patient>15])
m = lm(Mean_IFN_score_patient_neutro$x[n_patient>15]~mean_IFNg_CD3_per_patient[n_patient>15])
abline(coef(m),lwd=2,lty=2)
anova(m)
dev.off()
