options(stringsAsFactors = FALSE)
library(stringr) 

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
up = read.table(paste0(path,"Resumer_DEgenes_selection_UP.tab"), header = T)
down = read.table(paste0(path,"Resumer_DEgenes_DOWN.tab"), header = T)
all_deg = read.table(paste0(path,"Resumer_DEgenes.tab"), header = T)

colors=c("purple", "orange", "red3","forestgreen", "royalblue","deeppink", "grey")

# path = paste0(path,"Motif_20ntmax/parmis tout les genes/")
path = paste0(path,"Motif_intermed/Parmis tous les gènes/FIMO/")
dir.create(path, recursive=T,showWarnings=F)

# depuis les listes sortie par FIMO des différents motif selon le jeu de gene controle (1-5 liste de genes randomiser des 816 genes et la liste de touts les genes)
files = list.files(path, pattern = ".tsv")
for (f in files){
  # Filtre sur le nombre de motif présent dans les promoteurs
  prom_with_motif = read.table(paste0(path,f), header=T, sep="\t")
  
  tab_annot = merge(annotation, prom_with_motif, by.x = "ID", by.y = "sequence_name")
  write.table(tab_annot, paste0(path, "Motif_tab_",sub(".tsv","",f),".tab"), sep = "\t")  
  
  multi = prom_with_motif[is.element(prom_with_motif$sequence_name, unique(prom_with_motif$sequence_name[duplicated(prom_with_motif$sequence_name)])),]
  unique = setdiff(prom_with_motif$sequence_name,unique(multi$sequence_name))
  unique = prom_with_motif[is.element(prom_with_motif$sequence_name, unique),]
  prom_with_motif_nonredondant = prom_with_motif[!duplicated(prom_with_motif$sequence_name),]
  
  nb_motif = list(prom_with_motif, multi, unique, prom_with_motif_nonredondant)
  names(nb_motif) = c("tous","multi", "unique", "nonredondant")
  
  pdf(paste0(path, "Motif_",sub(".tsv","",f),".pdf"))
  for (n in names(nb_motif)){
    tab1 = nb_motif[[n]]
    
    # Filter sur les catégores de gènes ayant le motif
    ALL = annotation$ID
    MOTIF = ALL[is.element(ALL, tab1$sequence_name)]
    UP = MOTIF[is.element(MOTIF,up$ID)]
    UP_plus = UP[is.element(UP,unique(tab1$sequence_name[tab1$strand == "+"]))]
    UP_moins = UP[is.element(UP,unique(tab1$sequence_name[tab1$strand == "-"]))]
    
    filtre = list(ALL, MOTIF, UP, UP_plus, UP_moins)
    names(filtre) = c("ALL", "MOTIF", "UP", "UP_plus", "UP_moins")
    

    tab = data.frame(row.names = unique(annotation$EXPRESSION_PROFIL))
    for (i in names(filtre)){
      #### Definir le positionnement des motifs ####
      prom_with_motif_temp = prom_with_motif[which(is.element(prom_with_motif$sequence_name, filtre[[i]])),]
      hist(prom_with_motif_temp$start-150,
           breaks = nrow(tab),
           main = paste("Position du motif depuis TSS parmis", i, "de", n))
      
      #### Definir la proportion de gènes portant le motif ####
      annotation_temp = annotation[which(is.element(annotation$ID, filtre[[i]])),]
      tmp = c()
      for (j in unique(annotation$EXPRESSION_PROFIL)){
        tmp = c(tmp, length(grep(j, annotation_temp$EXPRESSION_PROFIL)))
      }
      tab = cbind(tab, tmp)
      colnames(tab)[ncol(tab)] = i
      
    }
    tab = tab[c(6,4,2,5,3,7,1),]

    barplot(as.matrix(tab),
            col=colors,border="white",
            cex=1,cex.axis=1,cex.lab=1,
            ylab="Gene number",
            main = paste("Promoteurs contenant le motif dans", n))
    legend("topright",legend=rev(rownames(tab)),col=rev(colors),bty="n",pch=15,cex=1.3)
    
    for (c in 1:ncol(tab)){
      tab[,c] = tab[,c]/sum(tab[,c])*100
    }
    
    barplot(as.matrix(tab),
            col=colors,border="white",
            cex=1,cex.axis=1,cex.lab=1,
            ylab="Gene proportion (%)",
            main = paste("Promoteurs contenant le motif dans", n))
    legend("topright",legend=rev(rownames(tab)),col=rev(colors),bty="n",pch=15,cex=1.3)
    
    # Parmis les gènes qui ont le motif, regarder les différentes cathégories ####
    filter2 = list(
      ALL = annotation$ID,
      MOTIF = tab1$sequence_name,
      PLUS = tab1$sequence_name[is.element(UP,unique(tab1$sequence_name[tab1$strand == "+"]))],
      MOINS = tab1$sequence_name[is.element(UP,unique(tab1$sequence_name[tab1$strand == "-"]))]
    )
    
    tab2 = data.frame(row.names = c("UP","DOWN","NA")) # Regarder la proportion de gènes UP et DOWN parmis les gènes porteur du motif
    tab3 = data.frame(row.names = unique(annotation$EXPRESSION_PROFIL)) # Regarder la proportion des profils d'expression des gènes porteur du motif
    for(f2 in names(filter2)){
      
      all_deg_temp = all_deg[which(is.element(all_deg$ID, filter2[[f2]])),]
      up_temp = up[which(is.element(up$ID, filter2[[f2]])),]
      down_temp = down[which(is.element(down$ID, filter2[[f2]])),]
      
      tmp2 = c(nrow(up_temp),nrow(down_temp),nrow(all_deg_temp)-nrow(up_temp)-nrow(down_temp))
      tab2 = cbind(tab2, tmp2)
      colnames(tab2)[ncol(tab2)] = f2
      
      annotation_temp2 = annotation[which(is.element(annotation$ID, filter2[[f2]])),]
      tmp3 = c()
      for (j in unique(annotation$EXPRESSION_PROFIL)){
        tmp3 = c(tmp3, length(grep(j, annotation_temp2$EXPRESSION_PROFIL)))
      }
      tab3 = cbind(tab3, tmp3)
      colnames(tab3)[ncol(tab3)] = f2

    }
    tab3 = tab3[c(6,4,2,5,3,7,1),]
    
    barplot(as.matrix(tab2),
            col=colors,border="white",
            cex=1,cex.axis=1,cex.lab=1,
            ylab="Motif number",
            main = paste("Motifs dans", n))
    legend("topright",legend=rownames(tab2), col=colors, bty="n",pch=15,cex=1.3)
    
    for (c in 1:ncol(tab2)){
      tab2[,c] = tab2[,c]/sum(tab2[,c])*100
    }
    
    barplot(as.matrix(tab2),
            col=colors,border="white",
            cex=1,cex.axis=1,cex.lab=1,
            ylab="Motif proportion (%)",
            main = paste("Motifs dans", n))
    legend("topright",legend=rownames(tab2), col=colors, bty="n",pch=15,cex=1.3)
    
    barplot(as.matrix(tab3),
            col=colors,border="white",
            cex=1,cex.axis=1,cex.lab=1,
            ylab="Motif number",
            main = paste("Motifs dans", n))
    legend("topright",legend=rev(rownames(tab3)),col=rev(colors),bty="n",pch=15,cex=1.3)
    
    for (c in 1:ncol(tab3)){
      tab3[,c] = tab3[,c]/sum(tab3[,c])*100
    }
    
    barplot(as.matrix(tab3),
            col=colors,border="white",
            cex=1,cex.axis=1,cex.lab=1,
            ylab="Motif number",
            main = paste("Motifs dans", n))
    legend("topright",legend=rev(rownames(tab3)),col=rev(colors),bty="n",pch=15,cex=1.3)
    
    
  }
  dev.off()
}
