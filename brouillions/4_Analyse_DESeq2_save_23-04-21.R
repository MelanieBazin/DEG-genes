source("2_Mise_en_forme_des_donnees.R")
source("3_Functions.R")

##### Chargement des librairies pour l'analyse #############
library(FactoMineR)
library(factoextra)
library(DESeq2)
library(gplots)
library(seqinr)
library("ggvenn")
library("RColorBrewer")
library("pheatmap")
#####


#### Initatiosation de variable ###
annotation_synonyms = annotation[annotation$SYNONYMS != "",]

Filtering= list()
Filtering= NULL
# ==> "Filtering" n'est pas définie,
# Réponse d'Olivier Pour certaine conditions 
# j'ai par exemple envie de retirer les gènes qui sont DE pendant une manip de silencing,
# ou entre plusieurs manip control... en gros des faux positifs


# Vecteur de couleur pour les heatmap
hmcol = colorRampPalette(brewer.pal(10,"RdBu"))(255)
#hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
hmcol = rev(hmcol)


##### Création des dossier pour ranger les données #############

base_img_dir=paste0("Analyse_DESeq2/",condition,"_", type,"/Images/")
dir.create(base_img_dir,recursive=T,showWarnings=F)

base_res_dir=paste0("Analyse_DESeq2/",condition,"_", type,"/")
dir.create(base_res_dir, recursive=T,showWarnings=F)
#####


##### Analyse DESeq2 ###########

# Filtrage des gènes avec trop peu de compatage (seuil arbitraire)
# countdata = countdata[rowSums(countdata) > 50,]

# Mise en forme des données
deseq = DESeqDataSetFromMatrix(countData = countdata,
                                colData  = infodata,
                                design   = ~ Conditions)


# Analyse DESeq2
deseq = DESeq(deseq)
#####

###############################################
# Reprise des variables et analyses d'Olivier #
###############################################
notAllZero = (rowSums(counts(deseq)) > 0 )
labels=colnames(countdata)


#### Moyenne des valeurs de comptage normalisées pour chaque point du timining####
##EARLY, INTERMEDIATE, LATE ##
time_points = c(paste("CTRL", timing_ctrl, sep = "_" ),paste(condition, timing_rnai, sep = "_" ))

# Extraction des données de comptage de DESeq2
geneNormCountsTable=counts(deseq,normalized=T)
meanGeneNormCountsTable =data.frame(ID=rownames(geneNormCountsTable))

countsTableNorm=as.data.frame(counts(deseq,normalized=TRUE)) #Mise en forme des donnée pour les heatmap

# Calcule des moyennes ==> Pourquoi utiliser cette partie du code au lieu de la fonction collapseReplicates de DESeq2 ?
for(p in unique(time_points)) {
  
  if(length(labels[time_points==p])==1)  { 
    meanGeneNormCountsTable[,p] =geneNormCountsTable[, labels[time_points==p]]
  } else { 
    meanGeneNormCountsTable[,p] =apply(geneNormCountsTable[, labels[time_points==p]],1,mean) 
  }
  
}
rownames(meanGeneNormCountsTable)=meanGeneNormCountsTable$ID
meanGeneNormCountsTable=meanGeneNormCountsTable[,-1]

#### Comparaison point par point des différents timing ####
comparisons = list(
  "VEG" = unique(time_points[grep("VEG",time_points)]),
  "EARLY" = unique(time_points[grep("EARLY",time_points)]),
  "INTER" = unique(time_points[grep("INTER",time_points)]),
  "LATE" = unique(time_points[grep("LATE",time_points)])
)


significant_up=list()
significant_down=list()

for(i in names(comparisons)) {
  c1=comparisons[[i]][1]
  c2=comparisons[[i]][2]
  resContrast=results(deseq,contrast=c("Conditions",c1, c2))
  
  resContrast=resContrast[notAllZero,]
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

  # Initialisation du 1er set de data : sn  as filtre => aucun gènes supprimés
  datasets=list("NoFilter"=resContrast_sig)
  
  # Créaction d'une liste contenant les donné de copatge pour tous les filtres

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
            text(res_vp[s,]$log2FoldChange+1,-log(res_vp[s,]$padj),annotation_synonyms$SYNONYMS[grep(s,annotation_synonyms$ID)])
          }
        }
      dev.off()
      
      ##### Variables "genes" non définie pour le volcanoplot #####
      # png(paste0(img_dir,"volcano_plot_",analysis_name,"_",i,"_",dname,"_annot_genes.png"), width = 6, height = 6, units = 'in', res = 300,family="ArialMT")
      #   plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col="gray",xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      #   for(id in names(genes)) {
      #     #if(res_vp[synonyms[s,]$ID,]$SIGNIFICANT) {
      #     points(res_vp[id,]$log2FoldChange,-log(res_vp[id,]$padj),col="black")
      #     text(res_vp[id,]$log2FoldChange+1,-log(res_vp[id,]$padj),genes[[id]])
      #     #}
      #   }
      # dev.off()
      #####
      
      
      ###### Création des heatmap pour cette comparaison ####
      if(nrow(res[res$REGULATION=="Up-regulated",]) >2) {
        regulation="Up-regulated"
        data=countsTableNorm[rownames(res[res$REGULATION==regulation,]),]
        hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
        
        png(paste(img_dir,"heatmap_",condition,"_",i,"_",dname,"_",regulation,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
          par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
          heatmap.2(as.matrix(log2(data+1)),
                    Rowv=as.dendrogram(hcGenes), Colv=NULL, 
                    scale="row", labRow="",trace="none", col = hmcol,
                    dendrogram = c("none"),main=paste(i," ",regulation,"N=",dim(data)[1]))
        dev.off()
      }
      
      if(nrow(res[res$REGULATION=="Down-regulated",]) >2) {
        regulation="Down-regulated"
        data=countsTableNorm[rownames(res[res$REGULATION==regulation,]),]
        hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
        
        png(paste(img_dir,"heatmap_",condition,"_",i,"_",dname,"_",regulation,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
          par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
          heatmap.2(as.matrix(log2(data+1)),
                    Rowv=as.dendrogram(hcGenes), Colv=NULL, 
                    scale="row", labRow="", trace="none", col = hmcol,
                    dendrogram = c("none"),main=paste(i," ",regulation,"N=",dim(data)[1]))
        dev.off()  
      }
      #####
    }
  }# Fin de la boucle pour chaque filtre  
}# Fin de la boucle pour chaque conditions ###




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
  
  pdf(paste(img_dir,"venn_",condition,"_",dname,"_Up_Down_both-regulated.pdf",sep=""))
    ggvenn(list("Up-regulated"=significant_up_ids,"Down-regulated"=significant_down_ids), stroke_size = 0.5, set_name_size = 4)
  dev.off()
  
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
  
  #### Représentation des gènes up et down reguler par digramme d eVenn, boxplot et heat map ####
  pdf(paste(img_dir,"venn_",condition,"_",dname,"_Up_and_Down-regulated.pdf",sep=""))
    ggvenn(list("Up-regulated"=significant_up_ids,"Down-regulated"=significant_down_ids), stroke_size = 0.5, set_name_size = 4)
  dev.off()
  
  
  regulation="Up-regulated"
  data=countsTableNorm[significant_up_ids,]
  hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
  
  png(paste(img_dir,"heatmap_",condition,"_",dname,"_",regulation,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
    par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
    heatmap.2(as.matrix(log2(data+1)),
              Rowv=as.dendrogram(hcGenes), Colv=NULL,
              scale="row", labRow="",trace="none", col = hmcol,
              dendrogram = c("none"),main=paste(regulation,"N=",dim(data)[1]))
  dev.off()
  
  pdf(paste(img_dir,"boxplot_",condition,"_",dname,"_",regulation,".pdf",sep=""))
    par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
    boxplot(log2(data+1),outline=F,las=2,ylab="Expression level (log2)",lwd=2,cex=1.3,cex.lab=1.3,cex.axis=1.3)
  dev.off()
  
  
  regulation="Down-regulated"
  data=countsTableNorm[significant_down_ids,]
  hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
  
  png(paste(img_dir,"heatmap_",condition,"_",dname,"_",regulation,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
    par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
    heatmap.2(as.matrix(log2(data+1)),
              Rowv=as.dendrogram(hcGenes), Colv=NULL, 
              scale="row", labRow="",trace="none", col = hmcol,
              dendrogram = c("none"),main=paste(regulation,"N=",dim(data)[1]))
  dev.off()
  
  pdf(paste(img_dir,"boxplot_",condition,"_",dname,"_",regulation,".pdf",sep=""))
    par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
    boxplot(log2(data+1),outline=F,las=2,ylab="Expression level (log2)",lwd=2,cex=1.3,cex.lab=1.3,cex.axis=1.3)
  dev.off()
  #####
  
  
  significant_ids=unique(c(significant_up_ids,significant_down_ids))
  
  if(dname == "NoFilter" & !is.null(names(Filtering))) {

    for(fname in names(Filtering)) {
      pdf(paste(img_dir,"venn_",analysis_name,"_",dname,"_And_",fname,".pdf",sep=""))
        ggvenn(list("Significant"=significant_ids,"Filtered"=Filtering[[fname]]), stroke_size = 0.5, set_name_size = 4)
      dev.off()
    }
    
    dv=list("Significant"=significant_ids)
    
    for(fname in names(Filtering)) {
      dv[[fname]]=Filtering[[fname]]
    }
    dv[["ExcisionComplexFiltering"]]=NULL
    pdf(paste(img_dir,"venn_",analysis_name,"_",dname,"_And_",paste(names(Filtering),collapse="_"),".pdf",sep=""))
      ggvenn(dv,simplify=T,  stroke_size = 0.5, set_name_size = 4)
    dev.off()
    
    dv=list("Significant"=significant_ids,"PGM_Filtering"=pgm_degenes$ID,"Controls"=ctl_rnai_degenes$ID)
    pdf(paste(img_dir,"venn_",analysis_name,"_",dname,"_And_PGM_Filtering_Controls.pdf",sep=""))
      ggvenn(dv,simplify=T,  stroke_size = 0.5, set_name_size = 4)
    dev.off()
    
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
    print(p)
    apply(autog_enrichment[2:4,c(p,"NB")],1,my_chi2,ctl=as.vector(autog_enrichment[1,c(p,"NB")]))
  }
  
  
  par(xpd=FALSE,mfrow=c(1,1))
  pdf(paste0(img_dir,"barplot_autogamy_proportion.pdf"),width=2,family="ArialMT")
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
  
  pdf(paste0(img_dir,"barplot_autogamy_cluster_proportion.pdf"),width=5)
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
  
  pdf(paste0(img_dir,"barplot_proportion_genes_with_IES.pdf"),width=5)
    barplot(prop,col=c("gray","indianred","darkgreen","dodgerblue"),border="white",cex=1.3,cex.axis=1.3,cex.lab=1.3,ylab="Proportion of genes with IES")
  dev.off()
}
