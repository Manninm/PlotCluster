PlotMyCluster<-function(Exp,Procs=4,Bootstraps=25,RowNames=1,skip=0,log=FALSE,SaveMatrices=FALSE,Features='Features',tissue='Tissue',Corr=FALSE,Boxplot=FALSE) {
	library(snow)
	library(pvclust)
	library(tools)	
	library(corrplot)
	library(ggfortify)
  library(phylogram)
  dist<-"euclidian"
  Exp<-read.table(Exp,row.names=RowNames,header=TRUE)
  Exp<-Exp[,-skip]
  print(names(Exp))
	Exp<-Exp[which(rowSums(Exp)>0),]
	cl<-makeCluster(Procs,type="SOCK")
	if (log) {
		Exp<-log2(Exp+1)
		res<-parPvclust(cl,Exp,method.hclust="ward.D", method.dist=dist,nboot=Bootstraps)
		pdf(paste(Features,tissue,"_ward.D_",Bootstraps,".pdf",sep=""), width=30, height=20)
		plot(res)
		dev.off()
		res2<-parPvclust(cl,Exp,method.hclust="ward.D2",method.dist=dist,nboot=Bootstraps)
		pdf(paste(Features,"_",tissue,"_","ward.D2.",Bootstraps,".pdf",sep=""), width=30, height=20)
		plot(res2)
		dev.off()
		stopCluster(cl)
	}
  else{
		res<-parPvclust(cl,Exp,method.hclust="ward.D", method.dist=dist,nboot=Bootstraps)		
		pdf(paste(Features,tissue,"_ward.D_",Bootstraps,".pdf",sep=""), width=30, height=20)		
		plot(res)		
		dev.off()		
		res2<-parPvclust(cl,Exp,method.hclust="ward.D2",method.dist=dist,nboot=Bootstraps)		
		pdf(paste(Features,"_",tissue,"_","ward.D2.",Bootstraps,".pdf",sep=""), width=30, height=20)		
		plot(res2)		
		dev.off()		
		stopCluster(cl)		
	}
	if(SaveMatrices){
		write.dendrogram(as.dendrogram(res),paste(Features,tissue,Bootstraps,'wardD','.gram',sep=''),quote=FALSE,sep='\t')
		write.dendrogram(as.dendrogram(res2),paste(Features,tissue,Bootstraps,'wardD2','.gram',sep=''),quote=FALSE,sep='\t')
	}
	else{
		warning('Proceeding without saving matrices')
	}
	if(Corr){
		co<-cor(Exp)
		pdf(paste("Corrplot_",,tissue,".TPM.pdf",sep=""), width=40, height=40)
		corrplot.mixed(co,upper="circle", lower="number", order="hclust",tl.pos="lt")
		dev.off()
	}
	else{
		warning('Proceeding without Correlation Plot')
	}
	if(Boxplot){
		boxExp<-Exp[which(rowSums(Exp)>0),]
		png(paste("BoxPlot_",Features,tissue,"_",".png",sep=""),width=2000,height=1000)
		par(cex.axis=0.9,mai=c(2.8,0.82,0.82,0.42))
		boxplot(boxExp,las=2,main="BoxploxNormalizedLog2")
		dev.off()
	}
	else{
		warning('Proceeding without BoxPlot')
	}
}

