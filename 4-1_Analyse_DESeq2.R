# Tourne sous R3
options(stringsAsFactors = FALSE)

notAllZero = (rowSums(counts(deseq)) > 0 )
labels=colnames(countdata)

# countsTableNorm=as.data.frame(counts(deseq,normalized=TRUE))
countsTableNorm = as.data.frame(assay(vst(deseq, blind = F)))

#### Comparaison point par point des différents timing ####
time_points = infodata$Condition
comparisons = list(
  "VEG" = unique(time_points[grep("VEG",time_points)]),
  "EARLY" = unique(time_points[grep("EARLY",time_points)]),
  "INTER" = unique(time_points[grep("INTER",time_points)]),
  "LATE" = setdiff(unique(time_points[grep("LATE",time_points)]),unique(time_points[grep("VERY_LATE",time_points)]))
)

regulation=c("Up-regulated","Down-regulated")

significant_up=list()
significant_down=list()

print(paste("Comparaison", RNAi, "et son controle" ))

for(i in names(comparisons)) {
  print(paste(RNAi, "au temps ",i))
  c1=comparisons[[i]][grep("ctrl",comparisons[[i]])]
  c2=comparisons[[i]][grep(RNAi,comparisons[[i]])]
  if(length(c2)==0|length(c1)==0){}else{
    resContrast=results(deseq,contrast=c("Condition",c2, c1))
    
    resContrast=resContrast[notAllZero,]
    
    ##### Répartition des p-value ####
    resContrast_sig = resContrast[ !is.na(resContrast$padj) , ]
    png(paste0(base_img_dir,"padj_density_",condition,"_",i,".png"))
    plot(density(resContrast_sig$padj), 
         main = paste(condition,i,"padj density"),
         xlab = "padj")
    dev.off()
    
    png(paste0(base_img_dir,"padj_log_density_",condition,"_",i,".png"))
    plot(density(-log(resContrast_sig$padj)), 
         main = paste(condition,i,"padj density"),
         xlab = "-log(padj)")
    dev.off()
    
    png(paste0(base_img_dir,"p-value_density_",condition,"_",i,".png"))
    plot(density(resContrast_sig$pvalue), 
         main = paste(condition,i,"p-value density"),
         xlab = "padj")
    dev.off()
    
    png(paste0(base_img_dir,"p-value_log_density_",condition,"_",i,".png"))
    plot(density(-log(resContrast_sig$pvalue)), 
         main = paste(condition,i,"p-value density"),
         xlab = "-log(padj)")
    dev.off()
    
    #### Resulta contrast ####
    
    resContrast_sig = resContrast[ !is.na(resContrast$padj) & resContrast$padj < pvalue , ]
    resContrast_sig = resContrast_sig[ resContrast_sig$log2FoldChange >= log2(FC) | resContrast_sig$log2FoldChange <= log2(1/FC), ]
    
    
    resContrast_sig = as.data.frame(resContrast_sig)
    
    # Initialisation du 1er set de data : sans filtre => aucun gènes supprimés
    datasets=list("NoFilter"=resContrast_sig)
    dname = names(datasets)
    
    #### Créaction des dossiers d'enregistrement ####
    img_dir=paste0(base_img_dir,"/")
    dir.create(img_dir,recursive=T,showWarnings=F)
    
    res_dir=paste0(base_res_dir,"/")
    dir.create(res_dir, showWarnings = FALSE,recursive=T)
    
    
    #### Création du tableau contennat les gènes significativement identifier comme déréguler ####
    res=as.data.frame(datasets[[dname]])
    
    
    res=res[,c("baseMean","log2FoldChange","padj")]
    res$REGULATION=ifelse(res$log2FoldChange>0,"Up-regulated","Down-regulated")
    res.annot= merge(res,annotation,by.x="row.names",by.y="ID")
    res.annot=res.annot[order(res.annot$padj), ]
    colnames(res.annot)[1]="ID"
    
    write.table(res.annot,paste0(res_dir,"DEgenes_",condition,"_",i,"_",dname,".tab"),sep="\t",quote=F,row.names=F)
    
    #### Etude des significatif ####
    
    
    significant_up[[dname]]=c(significant_up[[dname]],rownames(res[res$REGULATION=="Up-regulated",]))
    significant_down[[dname]]=c(significant_down[[dname]],rownames(res[res$REGULATION=="Down-regulated",]))    
    
    
    if(dim(res)[1] !=0) {
      res_vp=as.data.frame(resContrast)
      res_vp$SIGNIFICANT=FALSE
      res_vp[which(is.element(rownames(res_vp),rownames(res))),]$SIGNIFICANT=TRUE
      
      #### Création des volcanoplot pour cette comparaison ####
      png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,".png"), width = 6, height = 6, units = 'in', res = 300)
      plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col=ifelse(res_vp$SIGNIFICANT,"indianred","gray"),xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      dev.off()
      
      ## Tableau correspondant aux donné du volcanoplot ##
      all = length(annotation$ID)
      volca_tab = table(res$REGULATION)
      volca_tab = rbind(volca_tab,volca_tab/all*100)
      write.table(volca_tab, paste0(res_dir,"DEgenes_",condition,"_",i,"_",dname,"_volcanoplot.tab"),sep="\t",quote=F,row.names=F)
      
      # Volcanoplot avec les synonyme des gènes considérer comme significativement dérégulé
      png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,"_annot_synonyms.png"), width = 6, height = 6, units = 'in', res = 300)
      plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col="gray",xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      for(s in select_annotation$ID) {
        if(is.element(s, rownames(res_vp)) & res_vp[s,]$SIGNIFICANT) {
          points(res_vp[s,]$log2FoldChange,-log(res_vp[s,]$padj),col="green")
          text(res_vp[s,]$log2FoldChange+1,-log(res_vp[s,]$padj),select_annotation$NAME[grep(s,select_annotation$ID)])
        }
      }
      dev.off()
      
      #Volcanoplot avec nom choisis
      png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,"_annot_genes.png"), width = 6, height = 6, units = 'in', res = 300,family="ArialMT")
      plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col="gray",xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      for(id in select_ID) {
        #if(res_vp[synonyms[s,]$ID,]$SIGNIFICANT) {
        points(res_vp[id,]$log2FoldChange,-log(res_vp[id,]$padj),col="black")
        text(res_vp[id,]$log2FoldChange+1,-log(res_vp[id,]$padj),select_annotation$NAME[grep(id,select_annotation$ID)])
        #}
      }
      dev.off()
      
      #### Création des heatmap pour cette comparaison ####
      for (r in regulation){
        if(nrow(res[res$REGULATION==r,]) >2) {
          data=countsTableNorm[rownames(res[res$REGULATION==r,]),]
          hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
          
          png(paste(img_dir,"heatmap_",condition,"_",i,"_",dname,"_",r,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
          par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
          heatmap.2(as.matrix(log2(data+1)),
                    Rowv=as.dendrogram(hcGenes), Colv=NULL, 
                    scale="row", labRow="",trace="none", col = hmcol,
                    dendrogram = c("none"),main=paste(i," ",r,"N=",dim(data)[1]))
          dev.off()
        }
      }
      
    }
  }  
}


significant_up[[dname]]=c(significant_up[[dname]],rownames(res[res$REGULATION=="Up-regulated",]))
significant_down[[dname]]=c(significant_down[[dname]],rownames(res[res$REGULATION=="Down-regulated",]))    

print(paste(RNAi, "dataset formating" ))

#### Récupère les ID des gènes dérégulés identififé avec le filtre "dname" ####
significant_up_ids=unique(significant_up[[dname]])
significant_down_ids=unique(significant_down[[dname]])
significant_ids=unique(c(significant_up_ids,significant_down_ids))

#### Création du tableau contenant les gènes dérégulés ####
DEgenes=data.frame(ID=significant_up_ids,REGULATION="Up-regulated")
DEgenes=rbind(DEgenes,data.frame(ID=setdiff(significant_down_ids,significant_up_ids),REGULATION="Down-regulated"))
rownames(DEgenes)=DEgenes$ID
if(length(intersect(significant_up_ids,significant_down_ids))!=0) {
  DEgenes[intersect(significant_up_ids,significant_down_ids),]$REGULATION="Both"
}

write.table(merge(DEgenes,annotation,by="ID"),
            paste0(res_dir,"DEgenes_",condition,"_",dname,".tab"),sep="\t",quote=F,row.names=F,col.names=T)


#### Récuprération des ID des gènes UP et DOWN (exclu les "both") ####
significant_up_ids=DEgenes[DEgenes$REGULATION=="Up-regulated",]$ID
significant_down_ids=DEgenes[DEgenes$REGULATION=="Down-regulated",]$ID

for (r in regulation){
  data=countsTableNorm[significant_up_ids,]
  hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
  
  png(paste(img_dir,"heatmap_",condition,"_",dname,"_",r,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
  par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
  heatmap.2(as.matrix(log2(data+1)),
            Rowv=as.dendrogram(hcGenes), Colv=NULL,
            scale="row", labRow="",trace="none", col = hmcol,
            dendrogram = c("none"),main=paste(r,"N=",dim(data)[1]))
  dev.off()
  
  png(paste(img_dir,"boxplot_",condition,"_",dname,"_",r,".png",sep=""))
  par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
  boxplot(log2(data+1),outline=F,las=2,ylab="Expression level (log2)",lwd=2,cex=1.3,cex.lab=1.3,cex.axis=1.3)
  dev.off()
}


#### AUTOGAMY GENES ####
print(paste("Looking for autogamy genes distribution in ",RNAi, "dataset" ))

profiles=c("Early peak","Intermediate peak","Late peak" ,"Late induction", "Early repression","Late repression" ,"none")
profiles=c("Early repression","Late repression" ,"Early peak","Late peak" ,"Late induction","Intermediate peak","none")

autog_enrichment=data.frame()
autog_enrichment=rbind(autog_enrichment,get_autogamy_enrichment(rownames(wt_autogamy)))
autog_enrichment=rbind(autog_enrichment,get_autogamy_enrichment(significant_ids))
autog_enrichment=rbind(autog_enrichment,get_autogamy_enrichment(significant_up_ids))
autog_enrichment=rbind(autog_enrichment,get_autogamy_enrichment(significant_down_ids))
autog_enrichment[is.na(autog_enrichment)]=0

colnames(autog_enrichment)=c("NB","Autogamy","P Autogamy","Induced","P induced","Repressed","P Repressed",profiles,paste("P",profiles))
rownames(autog_enrichment)=c("ALL","SIG","UP","DOWN")

write.table(t(autog_enrichment),paste0(res_dir,"autog_enrichment.tab"),sep="\t",quote=F,row.names=T)


sink(paste0(res_dir,"autog_enrichment_chi2.tab"))

for(p in profiles) {
  print(paste0(p, " : DE -> SIG, UP, DOWN \n","DE /t", "NB_DE /t", "%DE /t", "ALL /t NB_ALL /t pvalue /t CHI2value /t signif"))
  apply(autog_enrichment[2:4,c(p,"NB")],1,my_chi2,ctl=as.vector(autog_enrichment[1,c(p,"NB")]))
}
sink()

par(xpd=FALSE,mfrow=c(1,1))
png(paste0(img_dir,"barplot_autogamy_proportion.png"),family="ArialMT")
barplot(c(autog_enrichment["ALL","P Autogamy"],autog_enrichment["SIG","P Autogamy"]),
        names.arg = c("ALL", "SIG"),
        col=c("gray","indianred"), border="white",ylim=c(0,100),cex=1.3,cex.axis=1.3,cex.lab=1.3,
        ylab="Developmental gene proportion (%)")
dev.off()

significant=merge(countsTableNorm[significant_ids,],wt_autogamy[,c("ID","EXPRESSION_PROFIL")],all.x=T,by="row.names")

prop=data.frame()
prop=rbind(prop,(table(wt_autogamy$EXPRESSION_PROFIL)/dim(wt_autogamy)[1])[profiles])
prop=rbind(prop,(table(significant$EXPRESSION_PROFIL)/dim(significant)[1])[profiles])
prop=rbind(prop,(table(significant[which(is.element(significant_ids,significant_up_ids)),]$EXPRESSION_PROFIL)/dim(significant[which(is.element(significant_ids,significant_up_ids)),])[1])[profiles])
prop=rbind(prop,(table(significant[which(is.element(significant_ids,significant_down_ids)),]$EXPRESSION_PROFIL)/dim(significant[which(is.element(significant_ids,significant_down_ids)),])[1])[profiles])
colnames(prop)=profiles
rownames(prop)=c("ALL","SIG","UP","DOWN")
prop[is.na(prop)]=0

colors=brewer.pal(length(profiles)-1,"Set1")
colors=rev(colors)
colors[1]="pink"
colors[length(profiles)]="white"
par(mfrow=c(1,1))

png(paste0(img_dir,"barplot_autogamy_cluster_proportion.png"))
barplot(t(prop),col=colors,border="white",cex=1.3,cex.axis=1.3,cex.lab=1.3,ylab="Gene proportion")
legend("topleft",legend=rev(profiles[-length(profiles)]),col=rev(colors[-length(profiles)]),bty="n",pch=15,cex=1.3)
dev.off()


#### IES ####
print(paste("Looking for IES distribution in ",RNAi, "dataset" ))

gene_with_IES_ids=annotation[annotation$NB_IES!=0,]$ID
prop=c(
  length(gene_with_IES_ids)/dim(annotation)[1],
  length(intersect(significant_ids,gene_with_IES_ids))/length(significant_ids),
  length(intersect(significant_up_ids,gene_with_IES_ids))/length(significant_up_ids),
  length(intersect(significant_down_ids,gene_with_IES_ids))/length(significant_down_ids)
)
names(prop)=c("ALL","SIG","UP","DOWN")

png(paste0(img_dir,"barplot_proportion_genes_with_IES.png"))
barplot(prop,col=c("gray","indianred","darkgreen","dodgerblue"),border="white",cex=1.3,cex.axis=1.3,cex.lab=1.3,ylab="Proportion of genes with IES")
dev.off()


