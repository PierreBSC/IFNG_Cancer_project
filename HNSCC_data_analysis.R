library(fifer)
library(pagoda2)
library(Matrix)
library(uwot)
library(igraph)
library(RColorBrewer)
library(ggplot2)

#Analysis of IFNg signature in Cancer cells of Head and Neck Squamous Cell Carcinoma 
#Original paper : "Single-cell transcriptomic analysis of primary and metastatic tumor ecosystems in head and neck cancer" 
#Single-cell technology used : Smart-Seq2 
#Please read first the Melanoma_data_analysis.R script for a good understanding of the approach used to identify the correct IFNg signature

color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","yellow","orange","red"))
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


##I)Pre-processing and loading of the data
##Data file downloaded from the GEO website (GSE123139)
#Link : https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103322&format=file&file=GSE103322%5FHNSCC%5Fall%5Fdata%2Etxt%2Egz

data_raw = read.delim("IFNG_project/HNSCC/GSE103322_HNSCC_all_data.txt",row.names = 1)
Condition = data_raw[1:5,]
Condition = as.data.frame(t(Condition))
data_raw = as(as.matrix(data_raw[-c(1:5),]),"dgCMatrix")

Donor = strsplit(colnames(data_raw),split = "_",fixed = T)
Donor = unlist(lapply(Donor,FUN = function(x) {x[1]}))

##II)Analysis of the whole dataset

#A)Pagoda2 analysis

r<-Pagoda2$new((data_raw),log.scale=F,modelType="raw")
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=50,n.odgenes=1000)
r$makeKnnGraph(k=15,type='PCA',center=T,distance='cosine',n.cores=10)
r$getKnnClusters(method=multilevel.community,type='PCA')
r$getKnnClusters(method=infomap.community,type='PCA',name = "infomap")

r$getDifferentialGenes(type="PCA",clusterType="community",upregulated.only=T,verbose=T,z.threshold=3)
r$getDifferentialGenes(type="PCA",clusterType="infomap",upregulated.only=T,verbose=T,z.threshold=3)

UMAP_plot=umap(r$reductions$PCA,verbose=T,metric="cosine",n_neighbors=15,spread = 2.5)
plot(UMAP_plot,pch=21,bg=string.to.colors(r$clusters$PCA$community)
     ,xaxt="n",yaxt="n",xlab="UMAP1",ylab="UMAP2",bty='n')

#What are the cancer cells ?
plot(UMAP_plot,pch=21,bg=string.to.colors(Condition$`classified  as cancer cell`)
     ,xaxt="n",yaxt="n",xlab="UMAP1",ylab="UMAP2",bty='n')

#Effects of the donor on RNA profile
plot(UMAP_plot,pch=21,bg=string.to.colors(Donor)
     ,xaxt="n",yaxt="n",xlab="UMAP1",ylab="UMAP2",bty='n')

#B)Study of IFNG distribution 

IFNG_positive_cells = r$counts[,"'IFNG'"]>0
table(IFNG_positive_cells,r$clusters$PCA$community)
T_cells = r$clusters$PCA$community==10

IFNG_level = aggregate(data_raw["'IFNG'",T_cells],FUN = mean,by=list(Donor[T_cells]))

Total_IFNG_prod = aggregate(as.matrix(data_raw["'IFNG'",]),by=list(r$clusters$PCA$community),FUN=sum)
Total_IFNG_prod= Total_IFNG_prod$V1
names(Total_IFNG_prod) = 1:length(Total_IFNG_prod)
Total_IFNG_prod = Total_IFNG_prod/sum(Total_IFNG_prod)
Total_IFNG_prod = Total_IFNG_prod[order(Total_IFNG_prod,decreasing = F)]

col_list =brewer.pal(7,"Paired")

pdf("IFNG_project/HNSCC/Figures/IFNG_distribution_CD3.pdf",useDingbats = F)
par(las=1,mar=c(5,7,4,22))
barplot(matrix(Total_IFNG_prod,ncol=1)*100,horiz = F,col=col_list,
        ylab="Total IFNG UMI production (%)",cex.lab=1.3,cex.axis = 1.3)
dev.off()


##III)Analysis of the cancer cells

#In the original publication the authors infered which cells are cancer cells...
#What are the effects of IFNg on those cells ?


#A)Pagoda2 analysis for those  cells

r_cancer_cells<-Pagoda2$new(round(data_raw[,Condition$`classified  as cancer cell`==1]),log.scale=F)
r_cancer_cells$adjustVariance(plot=T,gam.k=2)
r_cancer_cells$calculatePcaReduction(nPcs=100,n.odgenes=1e3)
r_cancer_cells$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine',n.cores=10)
r_cancer_cells$getKnnClusters(method=multilevel.community,type='PCA')
r_cancer_cells$getKnnClusters(method=infomap.community,type='PCA',name = "infomap")

r_cancer_cells$getDifferentialGenes(type="PCA",clusterType="community",upregulated.only=T,verbose=T,z.threshold=3)
r_cancer_cells$getDifferentialGenes(type="PCA",clusterType="infomap",upregulated.only=T,verbose=T,z.threshold=3)

UMAP_cancer_plot=umap(r_cancer_cells$reductions$PCA,verbose=T,metric="cosine",n_neighbors=15,spread = 3.5)
plot(UMAP_cancer_plot,pch=21,bg=string.to.colors(r_cancer_cells$clusters$PCA$community)
     ,xaxt="n",yaxt="n",xlab="UMAP1",ylab="UMAP2",bty='n')

#In the cancer cells : strong donor effect...
plot(UMAP_cancer_plot,pch=21,bg=string.to.colors(Donor[Condition$`classified  as cancer cell`==1])
     ,xaxt="n",yaxt="n",xlab="UMAP1",ylab="UMAP2",bty='n')


#B)Identifying the IFNg signature using the same strategy as for the MARS-seq melanoma data (see the other script)

cor_cancer = cor(as.matrix(r_cancer_cells$counts[,r_cancer_cells$getOdGenes(1000)]),method = "pearson")
Gene_clustering_neutro=hclust(dist(cor_cancer),method = "ward")
Gene_clustering_neutro=cutree(Gene_clustering_neutro,k = 30)

gene_env=c()
for (k in 1:length(unique(Gene_clustering_neutro))) {
  gene_env[[k]]=names(which(Gene_clustering_neutro==k))
}
names(gene_env)=1:30
gene_env_2 <- list2env(gene_env) # convert to an environment

r_cancer_cells$testPathwayOverdispersion(setenv = gene_env_2,type = "counts",verbose = T,
                               plot = F,max.pathway.size = 150,min.pathway.size = 1,recalculate.pca=T)
pathway_info=r_cancer_cells$misc$pathwayODInfo
pathway_info=pathway_info[order(pathway_info$cz,decreasing=T),]
View(pathway_info)
dispersed_genes=c()

for (i in (substr(rownames(pathway_info),7,10))) {
  u=r_cancer_cells$misc$pwpca[[i]]$xp$rotation
  u=u[unique(rownames(u)),]
  u=(u[order(u,decreasing = F)])ma
  dispersed_genes[[i]]=u
  print(i)
}
##Gene cluster 5 : IFNG response

par(las=1,mar = c(4,6,4,2))
IFN_score_neutro = as.numeric(r_cancer_cells$misc$pwpca$`5`$xp$scores)
IFN_score_neutro = IFN_score_neutro-min(IFN_score_neutro)
plot(UMAP_cancer_plot,pch=21,bg=color_convertion(IFN_score_neutro),
     xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
IFN_score = aggregate(IFN_score_neutro,by =list(Donor[Condition$`classified  as cancer cell`==1]),FUN=mean)
rownames(IFNG_level)=IFNG_level$Group.1
plot(IFNG_level[IFN_score$Group.1,"x"],IFN_score$x)
text(IFNG_level[IFN_score$Group.1,"x"],IFN_score$x,labels =IFN_score$Group.1 )

IFNG_score_table = data.frame(Score = IFN_score_neutro, Patient = Donor[Condition$`classified  as cancer cell`==1] )

#C)Visualisation of the results

#Taking only patients with at least 50 cancer cells

Number_cancer_cells = table(Donor[Condition$`classified  as cancer cell`==1])
Good_quality_patients = names(which(Number_cancer_cells>50))

#Distribution of IFNg score in cancer cells across patients
pdf("IFNG_project/HNSCC/Figures/Distribution_IFNG_score_cancer_cells.pdf",width = 10,height = 7,useDingbats = F)
p= ggplot(IFNG_score_table[IFNG_score_table$Patient%in%Good_quality_patients,], aes(x = Score, fill = Patient)) + geom_density(alpha = 0.4) + labs(x = "New x label")
p = p + labs(x = "IFNG score", y = "Density") +theme_bw()  + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                   panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
plot(p)
dev.off()

#Visualisation of IFNg score across UMAP space
pdf("IFNG_project/HNSCC/Figures/HNSCC_Ronan/Figures/UMAP_IFNG_score.pdf",width = 8,height = 8,useDingbats = F)
  plot(UMAP_cancer_plot,pch=21,bg=color_convertion(IFN_score_neutro),
    xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
dev.off()

#Visualisation of donor effect across UMAP space
pdf("IFNG_project/HNSCC/Figures/Figures/UMAP_Patient.pdf",width = 8,height = 8)
 plot(UMAP_cancer_plot,pch=21,bg=string.to.colors(Donor[Condition$`classified  as cancer cell`==1]),
      xaxt="n",yaxt="n",bty="n",xlab="UMAP 1",ylab="UMAP 2")
dev.off()
