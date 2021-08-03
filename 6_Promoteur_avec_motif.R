options(stringsAsFactors = FALSE)
library(seqinr)
library(stringr) 

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")

colors=c("purple", "orange", "red3","forestgreen", "royalblue","deeppink", "grey")
colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")

# Création des filtres
all_deg = read.table(paste0(path,"Resumer_DEgenes.tab"), header = T)[,c(1:9,14,17,20,23,26,29)]
up = all_deg$ID[which(all_deg$PGM_LATE_REGULATION == "Up-regulated" 
                      & all_deg$KU80c_LATE_REGULATION== "Up-regulated" 
                      & all_deg$XRCC4_LATE_REGULATION == "Up-regulated")]
ctip_down = all_deg$ID[unique(grep("Down-regulated",all_deg$CTIP_EARLY_REGULATION),
                              grep("Down-regulated",all_deg$CTIP_INTER_REGULATION))]
up_ezl1 = all_deg$ID[which(all_deg$EZL1_LATE_REGULATION == "Up-regulated" 
                           & all_deg$PGM_LATE_REGULATION == "Up-regulated" 
                           & all_deg$KU80c_LATE_REGULATION== "Up-regulated" 
                           & all_deg$XRCC4_LATE_REGULATION == "Up-regulated")]
selection = intersect(up, ctip_down)
selection_ezl1 = intersect(up_ezl1, ctip_down)

filtres_1 = list(
  TOUS = annotation,
  CTIP = annotation[is.element(annotation$ID, ctip_down),],
  PGM = annotation[is.element(annotation$ID, all_deg$ID[which(all_deg$PGM_LATE_REGULATION == "Up-regulated")]),],
  KU80c = annotation[is.element(annotation$ID, all_deg$ID[which(all_deg$KU80c_LATE_REGULATION == "Up-regulated")]),],
  XRCC4 = annotation[is.element(annotation$ID, all_deg$ID[which(all_deg$XRCC4_LATE_REGULATION == "Up-regulated")]),],
  EZL1 = annotation[is.element(annotation$ID, all_deg$ID[which(all_deg$EZL1_LATE_REGULATION == "Up-regulated")]),],
  UP = annotation[is.element(annotation$ID, up),],
  SELECT = annotation[is.element(annotation$ID, selection),],
  UP_EZL1 = annotation[is.element(annotation$ID, up_ezl1),],
  SELECT_EZL1 = annotation[is.element(annotation$ID, selection_ezl1),]
)

path_prom = "./DATA/Promoteur/"
path_save = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/"

files_prom = list.files(path_prom, pattern = ".fa")
files_motif = list.files(path_save, pattern = ".gff")
# 
# f = files_motif[1]
# sens = "+"

for (f in files_motif){
  print(f)
  p_name = sub(".gff","",sub(".gff3","",f))
  
  motif = read.table(paste0(path_save,f), header=F, sep="\t")
  colnames(motif)=colnamesgff3
  
  pdf(paste0(path_save,p_name,"_Repartion_motif_sans_filtre.pdf"))
  hist(motif$START-150,
       breaks = nrow(motif),
       main = paste(p_name, "\n","Position du motif sans filtres"))
  boxplot(motif$START-150, horizontal = T)
  dev.off()
  
  
  # Enlever les gènes qui ont plusieur motif et les contabilisés comme ayant un seul motif dans un sens
  motif_unique = motif[!duplicated(motif$ID),]
  duplicated = motif[is.element(motif$ID, unique(motif$ID[duplicated(motif$ID)])),]
  for(n in unique(duplicated$ID)){
    dup = duplicated[grep(n, duplicated$ID),]
    tab = table(dup$STRAND)
    if(is.null(tab[1])){
      sens = "+"
    }else if(is.null(tab[2])){
      sens = "-"
    }else{
      sens = "both"
    }
    motif_unique$STRAND[grep(n, motif_unique$ID)] = sens
  }
  
  
  
#### Boucle sur les filtres 1 : parmis les gènes ... combien ont un motif ####
  tab = data.frame(row.names = c("+","-","both", "none"))
  for(f1 in names(filtres_1)){
    print(f1)
    motif_filtre = motif_unique[is.element(motif_unique$ID, filtres_1[[f1]]$ID),]
    
    tmp = c()
    for(profils in c("+","-","both")){
      tmp = c(tmp,sum(grepl(profils, motif_unique$STRAND)))
    }
    tmp = c(tmp, length(setdiff(filtres_1[[f1]]$ID,motif_unique$ID)))
    tab = cbind(tab, tmp)
    colnames(tab)[ncol(tab)]=f1
  }

  write.table(tab, paste0(path_save,p_name,"_Filtres_dans_les_prom.tab"), row.names = T, sep = "\t")
  
  pdf(paste0(path_save, p_name, "_Filtres_dans_les_prom.pdf"))
  barplot(as.matrix(tab),
          col=c("blue", "red", "green", "grey"),border="white",
          cex=0.5,cex.axis=0.5,cex.lab=0.5,
          ylab="Gene number",
          main = paste("Promoteurs contenant le motif"))
  legend("topright",legend=rev(rownames(tab)),col=rev(c("blue", "red","green", "grey")),bty="n",pch=15,cex=0.5)
  
  for (c in 1:ncol(tab)){
    tab[,c] = tab[,c]/sum(tab[,c])*100
  }
  
  barplot(as.matrix(tab),
          col=c("blue", "red", "green","grey"),border="white",
          cex=0.5,cex.axis=0.5,cex.lab=0.5,
          ylab="Gene proportion (%)",
          main = paste("Promoteurs contenant le motif"))
  legend("topright",legend=rev(rownames(tab)),col=rev(c("blue", "red", "green","grey")),bty="n",pch=15,cex=0.5)
  dev.off()
  
  
  #### Boucle sur les filtres2 : parmis les gènes avec motif combien sont... ####  
  pdf(paste0(path_save, p_name, "_Filtres_parmis_les_motifs.pdf"))
  for(sens in c(unique(motif_unique$STRAND), "all")){
    print(sens)
    if (sens == "all"){}else{
      motif_sens = motif_unique[which(motif_unique$STRAND == sens),]
    }  
    
    filtres_2 = list(
      TOUS = motif_sens,
      CTIP = motif_sens[is.element(motif_sens$ID, ctip_down),],
      PGM = motif_sens[is.element(motif_sens$ID, all_deg$ID[which(all_deg$PGM_LATE_REGULATION == "Up-regulated")]),],
      KU80c = motif_sens[is.element(motif_sens$ID, all_deg$ID[which(all_deg$KU80c_LATE_REGULATION == "Up-regulated")]),],
      XRCC4 = motif_sens[is.element(motif_sens$ID, all_deg$ID[which(all_deg$XRCC4_LATE_REGULATION == "Up-regulated")]),],
      EZL1 = motif_sens[is.element(motif_sens$ID, all_deg$ID[which(all_deg$EZL1_LATE_REGULATION == "Up-regulated")]),],
      UP = motif_sens[is.element(motif_sens$ID, up),],
      SELECT = motif_sens[is.element(motif_sens$ID, selection),],
      UP_EZL1 = motif_sens[is.element(motif_sens$ID, up_ezl1),],
      SELECT_EZL1 = motif_sens[is.element(motif_sens$ID, selection_ezl1),]
    )
  
    
    tab = data.frame(row.names = unique(annotation$EXPRESSION_PROFIL))
    for(f2 in names(filtres_2)){
      print(f2)
      motif_filtre2 = filtres_2[[f2]] 
      if(nrow(motif_filtre2)>0){
        # Position du motif sur les promoteurs
        hist(motif_filtre2$START-150,
             breaks = nrow(motif_filtre2),
             main = paste(p_name, "\n","Position du motif", sens, f2))
        boxplot(motif_filtre2$START-150, horizontal = T)
        
        # Proportion de motifs en fonction des profils d'expression
        motif_annot = annotation[is.element(annotation$ID, unique(motif_filtre2$ID)),]
        tmp = c()
        for(profils in unique(annotation$EXPRESSION_PROFIL)){
          tmp = c(tmp,sum(grepl(profils, motif_annot$EXPRESSION_PROFIL)))
        }
        tab = cbind(tab, tmp)
        colnames(tab)[ncol(tab)]=f2
      }
    }
    tab = tab[order(rownames(tab)),]
    write.table(tab, paste0(path_save,p_name,"_Filtres_parmis_les_motifs",sens,".tab"), row.names = T, sep = "\t")
    
    barplot(as.matrix(tab),
            col=colors,border="white",
            cex=0.5,cex.axis=0.5,cex.lab=0.5,
            ylab="Gene number",
            main = paste("Promoteurs contenant le motif", sens , f2))
    legend("topright",legend=rev(rownames(tab)),col=rev(colors),bty="n",pch=15,cex=0.5)
    
    for (c in 1:ncol(tab)){
      tab[,c] = tab[,c]/sum(tab[,c])*100
    }
    
    barplot(as.matrix(tab),
            col=colors,border="white",
            cex=0.5,cex.axis=0.5,cex.lab=0.5,
            ylab="Gene proportion (%)",
            main = paste("Promoteurs contenant le motif", sens , f2))
    legend("topright",legend=rev(rownames(tab)),col=rev(colors),bty="n",pch=15,cex=0.5)
    
    
  }
  dev.off()
}  

