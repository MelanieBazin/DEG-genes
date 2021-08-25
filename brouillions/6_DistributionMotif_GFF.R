options(stringsAsFactors = FALSE)
library(seqinr)
library(stringr) 
# 
# f = files[2]
# n = names(nb_motif)[1]
# i = names(filtre)[1]
# f2 = names(filter2)[2]

type = "FIMO" #fuzznuc ou FIMO

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
path_fimo_files = paste0(path, "MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/",type,"/")

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
up = read.table(paste0(path,"Resumer_DEgenes_selection_UP.tab"), header = T)
down = read.table(paste0(path,"Resumer_DEgenes_DOWN.tab"), header = T)
all_deg = read.table(paste0(path,"Resumer_DEgenes.tab"), header = T)

colors=c("purple", "orange", "red3","forestgreen", "royalblue","deeppink", "grey")
colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")

path = paste0(path,"MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/",type,"/GFF/")
dir.create(path, recursive=T,showWarnings=F)

# depuis les listes sortie des différents motif selon le jeu de gene controle (1-5 liste de genes randomiser des 816 genes et la liste de touts les genes)
if (type == "fuzznuc"){
  files = list.files("./DATA/Promoteur/", pattern = ".gff")
}else if (type == "FIMO"){
  files = list.files(path_fimo_files, pattern = ".gff")
}

for (f in files){
  print(f)
  # Filtre sur le nombre de motif présent dans les promoteurs
  if (type == "fuzznuc"){
    prom_with_motif = read.table(paste0("./DATA/Promoteur/",f), header=F, sep="\t")
  }else if (type == "FIMO"){
    prom_with_motif = read.table(paste0(path_fimo_files,f), header=F, sep="\t")
  }
  colnames(prom_with_motif)=colnamesgff3
  write.table(prom_with_motif, paste0(path,"Motifs_TSS"),sep="\t",row.names=F,quote=F)
  
  tab_annot = merge(annotation, prom_with_motif, by = "ID")
  write.table(tab_annot, paste0(path, "Motif_tab_",sub(".gff3","",f),".tab"), sep = "\t")  
  
  multi = prom_with_motif[is.element(prom_with_motif$ID, unique(prom_with_motif$ID[duplicated(prom_with_motif$ID)])),]
  unique = setdiff(prom_with_motif$ID,unique(multi$ID))
  unique = prom_with_motif[is.element(prom_with_motif$ID, unique),]
  prom_with_motif_nonredondant = prom_with_motif[!duplicated(prom_with_motif$ID),]
  
  nb_motif = list(prom_with_motif, multi, unique, prom_with_motif_nonredondant)
  names(nb_motif) = c("tous","multi", "unique", "nonredondant")
  
  pdf(paste0(path, "Motif_",sub(".gff","",f),".pdf"))
  for (n in names(nb_motif)){
    print(n)
    if (!is.null(nrow(nb_motif[[n]]))){
    
    tab1 = nb_motif[[n]]
    
    # Filter sur les catégores de gènes ayant le motif
    ALL = annotation$ID
    MOTIF = ALL[is.element(ALL, tab1$ID)]
    UP = MOTIF[is.element(MOTIF,up$ID)]
    UP_plus = UP[is.element(UP,unique(tab1$ID[tab1$STRAND == "+"]))]
    UP_moins = UP[is.element(UP,unique(tab1$ID[tab1$STRAND == "-"]))]
    
    filtre = list(ALL, MOTIF, UP, UP_plus, UP_moins)
    names(filtre) = c("ALL", "MOTIF", "UP", "UP_plus", "UP_moins")
    
    
    tab = data.frame(row.names = unique(annotation$EXPRESSION_PROFIL))
    for (f1 in names(filtre)){
      print(f1)
      if (!is.null(nrow(filtre[[f1]]))){
      #### Definir le positionnement des motifs ####
      prom_with_motif_temp = prom_with_motif[which(is.element(prom_with_motif$ID, filtre[[f1]])),]
      hist(prom_with_motif_temp$START-150,
           breaks = nrow(tab),
           main = paste("Position du motif depuis TSS parmis", f1, "de", n))
      }
      #### Definir la proportion de gènes portant le motif ####
      annotation_temp = annotation[which(is.element(annotation$ID, filtre[[f1]])),]
      tmp = c()
      for (j in unique(annotation$EXPRESSION_PROFIL)){
        tmp = c(tmp, length(grep(j, annotation_temp$EXPRESSION_PROFIL)))
      }
      tab = cbind(tab, tmp)
      colnames(tab)[ncol(tab)] = f1
      
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
      MOTIF = tab1$ID,
      PLUS = tab1$ID[is.element(UP,unique(tab1$ID[tab1$STRAND == "+"]))],
      MOINS = tab1$ID[is.element(UP,unique(tab1$ID[tab1$STRAND == "-"]))]
    )
    
    tab2 = data.frame(row.names = c("UP","DOWN","NA")) # Regarder la proportion de gènes UP et DOWN parmis les gènes porteur du motif
    tab3 = data.frame(row.names = unique(annotation$EXPRESSION_PROFIL)) # Regarder la proportion des profils d'expression des gènes porteur du motif
    for(f2 in names(filter2)){
      print(f2)
      if (!is.null(nrow(filter2[[f2]]))){
      
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
    
    
  }}
  dev.off()
  }
  }

###### Notes #####
#get_sequence_from_gff.pl -gff ~/DATA/ptetraurelia_mac_51_annotation_v2.0.gff3 -genome /data/PARAMECIUM/REFERENCES/tetraurelia/ptetraurelia_mac_51.fa -upstream 150 -type CDS -out IN_MAC_upstream_150nt_CDS_start.tab > IN_MAC_upstream_150nt_CDS_start.fa
# promotors=read.table("SEARCH MOTIF/IN_MAC_upstream_150nt_CDS_start.tab",header=F, sep="\t")
# colnames(promotors)=colnamesgff3[1:4]

# Motif = AAAATC[AC]TT[AGT][AT]A[AT]TA[AT]TT

# fuzznuc -sequence IN_MAC_upstream_150nt_CDS_start.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]TA[AT]TT -outfile IN_MAC_Motif_upstream_150nt_CDS_start.fuzznuc.gff3 -complement -rformat2 gff
# fuzznuc -sequence IN_MAC_upstream_150nt_TSS.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]TA[AT]TT -outfile IN_MAC_Motif_upstream_150nt_TSS.fuzznuc.gff3 -complement -rformat2 gff
# fuzznuc -sequence IN_MAC_IES_upstream_150nt_CDS_start.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]TA[AT]TT -outfile IN_MAC_IES_Motif_upstream_150nt_CDS_start.fuzznuc.gff3 -complement -rformat2 gff
# fuzznuc -sequence IN_MAC_IES_upstream_150nt_TSS.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]TA[AT]TT -outfile IN_MAC_Motif_IES_upstream_150nt_TSS.fuzznuc.gff3 -complement -rformat2 gff
# fuzznuc -sequence PromUP_inter.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]TA[AT]TT -outfile IN_PromUP_TSS_Motif_upstream_150nt_TSS.fuzznuc.gff3 -complement -rformat2 gff

