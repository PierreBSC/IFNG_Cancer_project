library(fifer)

#Simple script to visualise the results of the Microscopy image analysis
#Please install the fifer package before running it !
#You need to have run the Matlab script Microscopy_Image_analysis_main.m

pdf("IFNG_project/Image_analysis_main/Figures_Biphoton.pdf",width = 9,height = 9,useDingbats = F)
par(las=1)
In_vivo_coloc_score = read.delim("IFNG_project/Image_analysis_main/Colocalisation_in_vivo.txt")
colnames(In_vivo_coloc_score) = c("Colocalisation","Condition")
prism.plots(Colocalisation~Condition,In_vivo_coloc_score,pch=21
            ,bg=string.to.colors(In_vivo_coloc_score$Condition,colors = c("grey","red")),col="black",
            cex=1.3,cex.lab=1.3,cex.axis=1.3,ylab="Colocalisation score",ylim=c(-1,1),
            xaxs="i",yaxs="i",main="Colocalisation score for in-vivo sample")
plotSigBars(Colocalisation~Condition,In_vivo_coloc_score)


In_vitro_coloc_score = read.delim("IFNG_project/Image_analysis_main/Colocalisation_in_vitro.txt")
colnames(In_vitro_coloc_score) = c("Colocalisation","Condition")
prism.plots(Colocalisation~Condition,In_vitro_coloc_score,pch=21
            ,bg=string.to.colors(In_vitro_coloc_score$Condition,colors = c("grey","red")),col="black",
            cex=1.3,cex.lab=1.3,cex.axis=1.3,ylab="Colocalisation score",ylim=c(-1,1),
            xaxs="i",yaxs="i",main="Colocalisation score for in-vitro sample")
plotSigBars(Colocalisation~Condition,In_vitro_coloc_score)

Association_distance = read.delim("IFNG_project/Image_analysis_main/T_cell_distance_effect.txt")
par(las=1)
plot(Association_distance,xlab="Distance to the closest T-cell",ylab="Colocalisation score",
     cex.lab=1.3,cex.axis=1.3,pch=21,bg="grey",xlim=c(0,60),ylim=c(-1,1),xaxs="i",yaxs="i",cex=1.5)
m = lm(Association_distance[,2]~Association_distance[,1])
abline(coef(m),lty=2,lwd=2,col="red")
legend("top",bty="n",legend = "R = 0.008 \nP = 0.931")


plot(Association_distance[Association_distance$Var1_1>0,],xlab="Distance to the closest T-cell",ylab="Colocalisation score",main="T cells removed",
     cex.lab=1.3,cex.axis=1.3,pch=21,bg="grey",xlim=c(0,60),ylim=c(-1,1),xaxs="i",yaxs="i",cex=1.5)
m = lm(Association_distance[,2]~Association_distance[,1],subset = Association_distance$Var1_1>0,)
abline(coef(m),lty=2,lwd=2,col="red")
legend("top",bty="n",legend = "R = -0.010 \nP = 0.965")
dev.off()

