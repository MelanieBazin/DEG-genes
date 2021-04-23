source("2_Mise_en_forme_des_donnees.R")

annotation_synonyms = annotation[annotation$SYNONYMS != "",]

##### Création des dossier pour ranger les données #############

base_img_dir=paste0("Analyse_DESeq2/",condition,"_", type,"/Images/")
dir.create(base_img_dir,recursive=T,showWarnings=F)

base_res_dir=paste0("Analyse_DESeq2/",condition,"_", type,"/")
dir.create(base_res_dir, showWarnings = FALSE,recursive=T)


##### Chargement des librairies pour l'analyse #############
library(FactoMineR)
library(factoextra)
library(DESeq2)

##### Analyse DESeq2 ###########
  
# Filtration des gènes avec trop peu de compatage (seuil arbitraire)
# countdata = countdata[rowSums(countdata) > 50,]

# Mise en forme des données
deseq = DESeqDataSetFromMatrix(countData = countdata,
                                colData  = infodata,
                                design   = ~ Conditions)


# Analyse DESeq2
deseq = DESeq(deseq)

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

## Comparerle controle au RNAi
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
  
  # ==> "Filtering" n'est pas définie,
  # Réponse d'Olivier Pour certaine conditions 
  # j'ai par exemple envie de retirer les gènes qui sont DE pendant une manip de silencing,
  # ou entre plusieurs manip control... en gros des faux positifs
  
  # for(fname in names(Filtering)) {  
  # 
  #   datasets[[fname]]=resContrast_sig[setdiff(rownames(resContrast_sig),Filtering[[fname]]),]
  #   #print(paste(i,fname,nrow(resContrast_sig),length(Filtering[[fname]]),dim(datasets[[fname]])[1]))
  # }
  
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
      png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,".png"), width = 6, height = 6, units = 'in', res = 300,family="ArialMT")
      plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col=ifelse(res_vp$SIGNIFICANT,"indianred","gray"),xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      dev.off()
      
      # Volcanoplot avec les synonyme des gènes considérer comme significativement dérégulé
      png(paste0(img_dir,"volcano_plot_",condition,"_",i,"_",dname,"_annot_synonyms.png"), width = 6, height = 6, units = 'in', res = 300,family="ArialMT")
      plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col="gray",xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      for(s in annotation_synonyms$ID) {
        if(is.element(s, rownames(res_vp)) & res_vp[s,]$SIGNIFICANT) {
          points(res_vp[s,]$log2FoldChange,-log(res_vp[s,]$padj),col="green")
          text(res_vp[s,]$log2FoldChange+1,-log(res_vp[s,]$padj),annotation_synonyms$SYNONYMS[grep(s,annotation_synonyms$ID)])
        }
      }
      dev.off()
      
      # Variables "genes" non définie pour le volcanoplot
      # png(paste0(img_dir,"volcano_plot_",analysis_name,"_",i,"_",dname,"_annot_genes.png"), width = 6, height = 6, units = 'in', res = 300,family="ArialMT")
      # plot(res_vp$log2FoldChange,-log(res_vp$padj),log="y",col="gray",xlab=paste0("log2(",c1,"/",c2,")"),ylab="-log(p-value)",pch=20,main=i,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      # for(id in names(genes)) {
      #   #if(res_vp[synonyms[s,]$ID,]$SIGNIFICANT) {
      #   points(res_vp[id,]$log2FoldChange,-log(res_vp[id,]$padj),col="black")
      #   text(res_vp[id,]$log2FoldChange+1,-log(res_vp[id,]$padj),genes[[id]])
      #   #}
      # }
      # dev.off()
      #####
      
      
      ##### Création des heat map pour cette comparaison ####
      if(length(rownames(res[res$REGULATION=="Up-regulated",])) >2) {
        regulation="Up-regulated"
        data=countsTableNorm[rownames(res[res$REGULATION==regulation,]),]
        hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
        
        png(paste(img_dir,"heatmap_",analysis_name,"_",i,"_",dname,"_",regulation,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
        par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
        heatmap.2(as.matrix(log2(data+1)),Rowv=as.dendrogram(hcGenes), Colv=NULL, scale="row", labRow="", col=hmcol,trace="none",dendrogram = c("none"),main=paste(i," ",regulation,"N=",dim(data)[1]))
        dev.off()
      }
      
      if(length(rownames(res[res$REGULATION=="Down-regulated",])) >2) {
        regulation="Down-regulated"
        data=countsTableNorm[rownames(res[res$REGULATION==regulation,]),]
        hcGenes=hclust(as.dist(1-cor(t(log2(data+1)), method="pearson")), method="complete")
        
        png(paste(img_dir,"heatmap_",analysis_name,"_",i,"_",dname,"_",regulation,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
        par(mar=c(8.1, 4.1, 4.1, 4.1), xpd=TRUE)
        heatmap.2(as.matrix(log2(data+1)),Rowv=as.dendrogram(hcGenes), Colv=NULL, scale="row", labRow="", col=hmcol,trace="none",dendrogram = c("none"),main=paste(i," ",regulation,"N=",dim(data)[1]))
        dev.off()  
      }
    }
  }  
  
}


