
###############################################
# Reprise des variables et analyses d'Olivier #
###############################################
notAllZero = (rowSums(counts(deseq)) > 0 )
labels=colnames(countdata)

countsTableNorm=as.data.frame(counts(deseq,normalized=TRUE))

#### Moyenne des valeurs de comptage normalisées pour chaque point du timining####
#Utiliser pour les heatmap seulement
meanGeneNormCountsTable = mean_data_tab

#### Comparaison point par point des différents timing ####
time_points = infodata$Condition
comparisons = list(
  "VEG" = unique(time_points[grep("VEG",time_points)]),
  "EARLY" = unique(time_points[grep("EARLY",time_points)]),
  "INTER" = unique(time_points[grep("INTER",time_points)]),
  "LATE" = setdiff(unique(time_points[grep("LATE",time_points)]),unique(time_points[grep("VERY_LATE",time_points)]))
  # "V_LATE" = unique(time_points[grep("VERY_LATE",time_points)])
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
    ######
    
    resContrast_sig = resContrast[ !is.na(resContrast$padj) & resContrast$padj < pvalue , ]
    resContrast_sig = resContrast_sig[ resContrast_sig$log2FoldChange >= log2(FC) | resContrast_sig$log2FoldChange <= log2(1/FC), ]
    
    
    resContrast_sig = as.data.frame(resContrast_sig)
    
    ############ Zone de test de FC adapté aux données CTIP #################
    # manuel_non_DEG=c("PTET.51.1.G0660118", "PTET.51.1.G1790042", "PTET.51.1.G0030302","PTET.51.1.G1740049","PTET.51.1.G0350134","PTET.51.1.G0230191","PTET.51.1.G0210241","PTET.51.1.G0770102","PTET.51.1.G0220178", "PTET.51.1.G0990073", "PTET.51.1.G0480035", "PTET.51.1.G0370168", "PTET.51.1.G0360062", "PTET.51.1.G1460025", "PTET.51.1.G0170355", "PTET.51.1.G0250220", "PTET.51.1.G0950175", "PTET.51.1.G0170354", "PTET.51.1.G0460033", "PTET.51.1.G1110086", "PTET.51.1.G1280115", "PTET.51.1.G1510135", "PTET.51.1.G1630015", "PTET.51.1.G0640197")
    # manuel_DEG_down = c("PTET.51.1.G0540024","PTET.51.1.G0900102","PTET.51.1.G0020380","PTET.51.1.G0120328","PTET.51.1.G0870035","PTET.51.1.G0360066","PTET.51.1.G1010039","PTET.51.1.G0150242","PTET.51.1.G0980137","PTET.51.1.G0170233","PTET.51.1.G0120245","PTET.51.1.G0380022","PTET.51.1.G0230222","PTET.51.1.G0920155","PTET.51.1.G0490162","PTET.51.1.G0070121","PTET.51.1.G0610198","PTET.51.1.G1200062","PTET.51.1.G0680113","PTET.51.1.G1300067","PTET.51.1.G1400105","PTET.51.1.G0030168","PTET.51.1.G0380073","PTET.51.1.G0450225","PTET.51.1.G0620215","PTET.51.1.G0050231","PTET.51.1.G0020217","PTET.51.1.G0350154","PTET.51.1.G0350166","PTET.51.1.G0240239","PTET.51.1.G0010374","PTET.51.1.G0480099","PTET.51.1.G0080368","PTET.51.1.G0210235","PTET.51.1.G0110289","PTET.51.1.G1530110","PTET.51.1.G0260051","PTET.51.1.G1150114","PTET.51.1.G0590028","PTET.51.1.G0340197","PTET.51.1.G0110267","PTET.51.1.G1140146","PTET.51.1.G0350167","PTET.51.1.G0370136","PTET.51.1.G0210213","PTET.51.1.G0220140","PTET.51.1.G0360089","PTET.51.1.G0010451","PTET.51.1.G0020335","PTET.51.1.G0250013","PTET.51.1.G1330044","PTET.51.1.G0060034","PTET.51.1.G0530071")
    # manuel_non_DEG = as.data.frame(resContrast)[is.element(row.names(as.data.frame(resContrast)),manuel_non_DEG),]
    # summary(manuel_non_DEG)
    # manuel_DEG_down = as.data.frame(resContrast)[is.element(row.names(as.data.frame(resContrast)), manuel_DEG_down),]
    # summary(manuel_DEG_down)
    # 
    # log2FC = list(manuel_non_DEG$log2FoldChange, manuel_DEG_down$log2FoldChange)
    # names(log2FC) =c("Non DEG", "DEG")
    # boxplot(log2FC, main = paste("log2FC","\n",i), horizontal=F)
    ##############
    
    # Initialisation du 1er set de data : sans filtre => aucun gènes supprimés
    datasets=list("NoFilter"=resContrast_sig)
    
    # Créaction d'une liste contenant les donné de compatge pour tous les filtres
    
    for(fname in names(Filtering)) {
      datasets[[fname]]=resContrast_sig[setdiff(rownames(resContrast_sig),Filtering[[fname]]),]
      #print(paste(i,fname,nrow(resContrast_sig),length(Filtering[[fname]]),dim(datasets[[fname]])[1]))
    }
    
    # Pour chacun des jeux de données filtrés ou non
    for(dname in names(datasets)) {
      
      #### Créaction des dossiers d'enregistrement ####
      img_dir=paste0(base_img_dir,"/",dname,"/")
      dir.create(img_dir,recursive=T,showWarnings=F)
      
      res_dir=paste0(base_res_dir,"/",dname,"/")
      dir.create(res_dir, showWarnings = FALSE,recursive=T)
      #####
      
      #### Création du tableau contennat les gènes significativement identifier comme déréguler ####
      res=as.data.frame(datasets[[dname]])
      #print(paste(i,dname,dim(res)[1]))
      
      res=res[,c("baseMean","log2FoldChange","padj")]
      res$REGULATION=ifelse(res$log2FoldChange>0,"Up-regulated","Down-regulated")
      res.annot= merge(res,annotation,by.x="row.names",by.y="ID")
      res.annot=res.annot[order(res.annot$padj), ]
      colnames(res.annot)[1]="ID"
      
      write.table(res.annot,paste0(res_dir,"DEgenes_",condition,"_",i,"_",dname,".tab"),sep="\t",quote=F,row.names=F)
      

      
      #####
      
      
      significant_up[[dname]]=c(significant_up[[dname]],rownames(res[res$REGULATION=="Up-regulated",]))
      significant_down[[dname]]=c(significant_down[[dname]],rownames(res[res$REGULATION=="Down-regulated",]))    
      
      
      if(dim(res)[1] !=0) {
        res_vp=as.data.frame(resContrast)
        res_vp$SIGNIFICANT=FALSE
        res_vp[which(is.element(rownames(res_vp),rownames(res))),]$SIGNIFICANT=TRUE
        
        ##### Création des volcanoplot pour cette comparaison ####
        png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,".png"), width = 6, height = 6, units = 'in', res = 300)
        plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col=ifelse(res_vp$SIGNIFICANT,"indianred","gray"),xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
        dev.off()
        
        # Volcanoplot avec les synonyme des gènes considérer comme significativement dérégulé
        png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,"_annot_synonyms.png"), width = 6, height = 6, units = 'in', res = 300)
        plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col="gray",xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
        for(s in annotation_synonyms$ID) {
          if(is.element(s, rownames(res_vp)) & res_vp[s,]$SIGNIFICANT) {
            points(res_vp[s,]$log2FoldChange,-log(res_vp[s,]$padj),col="green")
            text(res_vp[s,]$log2FoldChange+1,-log(res_vp[s,]$padj),annotation_synonyms$NAME[grep(s,annotation_synonyms$ID)])
          }
        }
        dev.off()
        
        ##### Volcanoplot avec nom choisis #####
        png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,"_annot_genes.png"), width = 6, height = 6, units = 'in', res = 300,family="ArialMT")
          plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col="gray",xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
          for(id in select_ID) {
            #if(res_vp[synonyms[s,]$ID,]$SIGNIFICANT) {
            points(res_vp[id,]$log2FoldChange,-log(res_vp[id,]$padj),col="black")
            text(res_vp[id,]$log2FoldChange+1,-log(res_vp[id,]$padj),annotation_synonyms$NAME[grep(id,annotation_synonyms$ID)])
            #}
          }
        dev.off()
        #####
        
        
        ###### Création des heatmap pour cette comparaison ####
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
}


print(paste("Comparaison des dereguler en ", RNAi, "a d'autre donnees" ))

#### Boucle sur tous les filtres utilisés ####
for(dname in unique(c(names(significant_up),names(significant_down)))) {
  
  ##### Créaction des dossiers d'enregistrement ####
  # Crée les dossier s'ils n'ont pas été créer par la boucle précédente
  img_dir=paste0(base_img_dir,"/",dname,"/")
  dir.create(img_dir,recursive=T,showWarnings=F)
  
  res_dir=paste0(base_res_dir,"/",dname,"/")
  dir.create(res_dir,recursive=T,showWarnings=F)
  #####
  
  # Récupère les ID des gènes dérégulés identififé avec le filtre "dname"
  significant_up_ids=unique(significant_up[[dname]])
  significant_down_ids=unique(significant_down[[dname]])
  
  # png(paste(img_dir,"venn_",condition,"_",dname,"_Up_Down_both-regulated.png",sep=""))
  # ggvenn(list("Up-regulated"=significant_up_ids,"Down-regulated"=significant_down_ids), stroke_size = 0.5, set_name_size = 4)
  # dev.off()
  
  ##### Création du tableau contenant les gènes dérégulés ####
  DEgenes=data.frame(ID=significant_up_ids,REGULATION="Up-regulated")
  DEgenes=rbind(DEgenes,data.frame(ID=setdiff(significant_down_ids,significant_up_ids),REGULATION="Down-regulated"))
  rownames(DEgenes)=DEgenes$ID
  if(length(intersect(significant_up_ids,significant_down_ids))!=0) {
    DEgenes[intersect(significant_up_ids,significant_down_ids),]$REGULATION="Both"
  }
  #table(DEgenes$REGULATION)
  write.table(merge(DEgenes,annotation,by="ID"),
              paste0(res_dir,"DEgenes_",condition,"_",dname,".tab"),sep="\t",quote=F,row.names=F,col.names=T)
  #####
  
  # Récuprération des ID des gènes UP et DOWN (exclu les "both") 
  significant_up_ids=DEgenes[DEgenes$REGULATION=="Up-regulated",]$ID
  significant_down_ids=DEgenes[DEgenes$REGULATION=="Down-regulated",]$ID
  
  #### Représentation des gènes up et down reguler par digramme de Venn, boxplot et heat map ####
  # png(paste(img_dir,"venn_",condition,"_",dname,"_Up_and_Down-regulated.png",sep=""))
  # ggvenn(list("Up-regulated"=significant_up_ids,"Down-regulated"=significant_down_ids), stroke_size = 0.5, set_name_size = 4)
  # dev.off()
  
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
  
  #####
  
  significant_ids=unique(c(significant_up_ids,significant_down_ids))
  
  if(dname == "NoFilter" & !is.null(names(Filtering))) {
    
    # for(fname in names(Filtering)) {
    #   png(paste(img_dir,"venn_",condition,"_",dname,"_And_",fname,".png",sep=""))
    #   ggvenn(list("Significant"=significant_ids,"Filtered"=Filtering[[fname]]), stroke_size = 0.5, set_name_size = 4)
    #   dev.off()
    # }
    
    dv=list("Significant"=significant_ids)
    
    for(fname in names(Filtering)) {
      dv[[fname]]=Filtering[[fname]]
    }
    dv[["ExcisionComplexFiltering"]]=NULL
    # png(paste(img_dir,"venn_",condition,"_",dname,"_And_",paste(names(Filtering),collapse="_"),".png",sep=""))
    # ggvenn(dv,simplify=T,  stroke_size = 0.5, set_name_size = 4)
    # dev.off()
    
    dv=list("Significant"=significant_ids,"PGM_Filtering"=pgm_degenes$ID,"Controls"=ctl_rnai_degenes$ID)
    # png(paste(img_dir,"venn_",condition,"_",dname,"_And_PGM_Filtering_Controls.png",sep=""))
    # ggvenn(dv,simplify=T,  stroke_size = 0.5, set_name_size = 4)
    # dev.off()
    
  }
  
  ##################
  # AUTOGAMY GENES #
  ##################
  
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
  #t(autog_enrichment)
  
  write.table(t(autog_enrichment),paste0(res_dir,"autog_enrichment.tab"),sep="\t",quote=F,row.names=T)
  
  # apply(autog_enrichment[2:4,c("Autogamy","NB")],1,my_chi2,ctl=as.vector(autog_enrichment[1,c("Autogamy","NB")]))
  # apply(autog_enrichment[2:4,c("Induced","Autogamy")],1,my_chi2,ctl=as.vector(autog_enrichment[1,c("Induced","Autogamy")]))
  # apply(autog_enrichment[2:4,c("Repressed","Autogamy")],1,my_chi2,ctl=as.vector(autog_enrichment[1,c("Repressed","Autogamy")]))
  

  for(p in profiles) {
    #print(p)
    apply(autog_enrichment[2:4,c(p,"NB")],1,my_chi2,ctl=as.vector(autog_enrichment[1,c(p,"NB")]))
  }


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
  
  #############
  #    IES    #
  #############
  
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
}

