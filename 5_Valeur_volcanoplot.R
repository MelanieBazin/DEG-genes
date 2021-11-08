

ALL = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/Resumer_DEgenes.tab"), sep = "\t", h=T )
ALL = ALL[-41534,]
all = nrow(ALL)

PGM = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","PGM","/NoFilter/DEgenes_tout_LATE_NoFilter.tab"), sep = "\t", h=T )
KU80c = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","KU80c","/NoFilter/DEgenes_tout_LATE_NoFilter.tab"), sep = "\t", h=T)
XRCC4 = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","XRCC4","/NoFilter/DEgenes_tout_LATE_NoFilter.tab"), sep = "\t", h=T )
EZL1 = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","EZL1","/NoFilter/DEgenes_tout_LATE_NoFilter.tab"), sep = "\t", h=T )

UP = ALL[ALL$PGM_LATE_REGULATION == "Up-regulated" &
           ALL$KU80c_LATE_REGULATION == "Up-regulated" &
           ALL$XRCC4_LATE_REGULATION == "Up-regulated" &
           ALL$EZL1_LATE_REGULATION == "Up-regulated",]

DOWN = ALL[ALL$PGM_LATE_REGULATION == "Down-regulated" &
             ALL$KU80c_LATE_REGULATION == "Down-regulated" &
             ALL$XRCC4_LATE_REGULATION == "Down-regulated" &
             ALL$EZL1_LATE_REGULATION == "Down-regulated",]


CTIP_e = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","CTIP","/NoFilter/DEgenes_tout_EARLY_NoFilter.tab"), sep = "\t", h=T )
CTIP_i = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","CTIP","/NoFilter/DEgenes_tout_INTER_NoFilter.tab"), sep = "\t", h=T )

CTIP_DOWN = ALL[ALL$CTIP_EARLY_REGULATION == "Down-regulated" |
                  ALL$CTIP_INTER_REGULATION == "Down-regulated",]

CTIP_UP = ALL[ALL$CTIP_EARLY_REGULATION == "Up-regulated" |
                ALL$CTIP_INTER_REGULATION == "Up-regulated",]

MOTIF = read.table("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/FIMO2_IN_MAC_CDS.tsv",sep = "\t", h=T)
MOTIF$start =  MOTIF$start-150
MOTIF$stop = MOTIF$stop-150

FILTRE_ARNi = list(
  ALL = ALL,
  PGM_UP = PGM[PGM$REGULATION == "Up-regulated",],
  KU80c_UP = KU80c[KU80c$REGULATION == "Up-regulated",],
  XRCC4_UP = XRCC4[XRCC4$REGULATION == "Up-regulated",],
  EZL1_UP = EZL1[EZL1$REGULATION == "Up-regulated",],
  UP = UP,
  CTIP_DOWN = CTIP_DOWN,
  UP_DOWN = ALL[which(is.element(ALL$ID,intersect(UP$ID, CTIP_DOWN$ID))),],
  PGM_DOWN = PGM[PGM$REGULATION == "Down-regulated",],
  KU80c_DOWN = KU80c[KU80c$REGULATION == "Down-regulated",],
  XRCC4_DOWN = XRCC4[XRCC4$REGULATION == "Down-regulated",],
  EZL1_DOWN = EZL1[EZL1$REGULATION == "Down-regulated",],
  DOWN = DOWN,
  CTIP_UP = CTIP_UP,
  DOWN_UP = ALL[which(is.element(ALL$ID,intersect(DOWN$ID, CTIP_UP$ID))),]
)


ALL2 = ALL
colnames(ALL2)[1] = "sequence_name"
FILTRE_MOTIF = list(
  ALL = ALL2,
  MOTIF = MOTIF,
  MOTIF_pos= MOTIF[MOTIF$start < -50 & MOTIF$start > -80,],
  MOTIF_plus = MOTIF[MOTIF$strand == "+",],
  MOTIF_moins = MOTIF[MOTIF$strand == "-",],
  MOTIF_pos_plus = MOTIF[MOTIF$strand == "+" & MOTIF$start < -50 & MOTIF$start > -80,],
  MOTIF_pos_moins = MOTIF[MOTIF$strand == "-" & MOTIF$start < -50 & MOTIF$start > -80,]
  
)
#### FONCTIONS ####
row_order = c("Early peak", "Intermediate peak", "Late peak", "Early repression" ,"Late induction", "Late repression", "none" )
colors = c("purple3","red2","chartreuse4","dodgerblue3","deeppink","darkorange","snow3")
save_path = "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/Plots/"

  
### Calcule des nombre de gène up et down en rouge sur les volcanoplot ####
# tab = table(PGM$REGULATION)
# print("PGM")
# print(rbind(tab,tab/all*100))
# 
# tab = table(KU80c$REGULATION)
# print("KU80c")
# print(rbind(tab,tab/all*100))
# 
# tab = table(XRCC4$REGULATION)
# print("XRCC4")
# print(rbind(tab,tab/all*100))
# 
# tab = table(EZL1$REGULATION)
# print("EZL1")
# print(rbind(tab,tab/all*100))
# 
# tab = table(CTIP_e$REGULATION)
# print("CTIP_e")
# print(rbind(tab,tab/all*100))
# 
# tab = table(CTIP_i$REGULATION)
# print("CTIP_i")
# print(rbind(tab,tab/all*100))


### Calcule des nombres de gènes dans chaque RNAi par profils ###
# f_motif = names(FILTRE_MOTIF)[1]
# f_arni = names(FILTRE_ARNi)[2]

sink(file = NULL)
t = table(ALL$EXPRESSION_PROFIL)
TAB= array(dim = c(length(t),length(names(FILTRE_MOTIF))*length(names(FILTRE_ARNi))))
rownames(TAB) = names(t)


col = 0
name = c()
for(f_motif in names(FILTRE_MOTIF)){
  dir.create(paste0(save_path,f_motif,"/"),recursive=T,showWarnings=F)
  tab1 = FILTRE_MOTIF[[f_motif]]

  for (f_arni in names(FILTRE_ARNi)){
    tab2 = FILTRE_ARNi[[f_arni]]
    tab = tab2[which(is.element(tab2$ID, tab1$sequence_name)),]
    t= table(tab$EXPRESSION_PROFIL)
    col = col+1
    
    for(r in rownames(TAB)){
      l = grep(r, names(t))
      if (length(l)>0){
        TAB[r,col] = t[l]
      }else{
        TAB[r,col] = 0
      }
    }
    
    name = c(name, paste(f_arni,f_motif, sep = "_"))

    if (f_arni != names(FILTRE_ARNi)[1]){
      sink(file = paste0(save_path,f_motif,"/",f_arni,"_Profil_Chi2.txt"))
      print(paste(f_motif,f_arni))
      sTAB = cbind(TAB[,1], TAB[,col])
      colnames(sTAB) = c(names(FILTRE_ARNi)[1],f_arni)
      sTAB = sTAB[row_order,]
      pvalue = c()
      star = c()
      for (p in rownames(sTAB)){
        stab = matrix(data = NA, nrow = 2, ncol = 2)
        stab[,1] = c(sTAB[p,2],sTAB[p,"ALL"]) 
        stab[,2] = c(sum(sTAB[,2])-sTAB[p,2], sum(sTAB[,"ALL"])-sTAB[p,"ALL"])
        
        chi2 = chisq.test(stab)
        chi2$pvalue = c(pvalue, chi2$p.value)
        print(p)
        print(chi2)
        if (is.na(chi2$p.value)){
          star = c(star," ")
        }else if (chi2$p.value < 1e-200){
          star = c(star,"****")
        }else if (chi2$p.value < 1e-100){
          star = c(star,"***")
        }else if (chi2$p.value < 1e-20){
          star = c(star,"**")
        }else if (chi2$p.value < 1e-10){
          star = c(star,"*")
        }else {
          star = c(star,"ns")
        }
      }  
      sTAB = cbind(sTAB, pvalue, star)
      write.table(sTAB, paste0(save_path,f_motif,"/",f_arni,"_Profils_stat.tab"), sep = "\t")
      
      # Faire les représentataions graphique
      data_plot = cbind(TAB[,1]/sum(TAB[,1])*100, TAB[,col]/TAB[,1]*100)
      colnames(data_plot) = c(names(FILTRE_ARNi)[1],f_arni)
      data_plot = data_plot[row_order,]
      
      png(paste0(save_path,f_motif,"/",f_arni,"_Profils_barplot.png"),width = 800, height = 500)
      barplot(t(data_plot),
              col = c("#353436","#1b98e0"),
              beside = TRUE,
              main = f_arni,
              ylab = "% de gènes")
      
      legend("topleft",
             legend = c("Reference","ARNi"),
             fill = c("#353436","#1b98e0"),
             bty = "n")
      
      dev.off()

      if (sum(TAB[,col]) > 0 ){
      ## Diagame emplié ##
      data_plot = cbind(TAB[,1]/sum(TAB[,1])*100, TAB[,col]/sum(TAB[,col])*100)
      colnames(data_plot) = c(names(FILTRE_ARNi)[1],f_arni)
      data_plot = data_plot[row_order,]      
      
      png(paste0(save_path,f_motif,"/",f_arni,"_Profils_barplot_pile.png"),width = 200, height = 500)
       barplot(data_plot,
              col = colors,
              main = f_arni,
              ylab = "% de gènes")

      dev.off()
      }
      sink(file = NULL)
      }

  }
  
  
}



colnames(TAB) = name

write.table(TAB, paste0(save_path,"Expression_profil_FIMO2_IN_MAC_CDS_table.tab"), sep = "\t")

png(paste0(save_path,"/Legend_pile_Profils.png"),width = 500, height = 500)
plot(NULL)
legend("topleft",
       legend = rev(row_order),
       fill = rev(colors),
       bty = "n")
dev.off()


#### Regarder le nombre de gènes IES + parmis les différnets cathégories ###
sink(file = NULL)
TAB= array(dim = c(2,length(names(FILTRE_MOTIF))*length(names(FILTRE_ARNi))))
rownames(TAB) = c("SANS_IES","AVEC_IES")

# f_motif = "ALL"
# f_arni = "ALL"
col = 0
name = c()
for(f_motif in names(FILTRE_MOTIF)){
  tab1 = FILTRE_MOTIF[[f_motif]]
  
  for (f_arni in names(FILTRE_ARNi)){
    tab2 = FILTRE_ARNi[[f_arni]]
    tab = tab2[which(is.element(tab2$ID, tab1$sequence_name)),]
    
    col = col+1
    
    TAB["SANS_IES",col] = sum(tab$NB_IES ==0)
    TAB["AVEC_IES",col] = sum(tab$NB_IES !=0)

    name = c(name, paste(f_arni,f_motif, sep = "_"))
    
    if (f_arni != names(FILTRE_ARNi)[1]){
      sink(file = paste0(save_path,f_motif,"/",f_arni,"_IES_Chi2.txt"))
      print(paste(f_motif,f_arni))
      
      sTAB = cbind(TAB[,1], TAB[,col])
      colnames(sTAB) = c(names(FILTRE_ARNi)[1],f_arni)
      pvalue = c()
      star = c()
      for (p in rownames(sTAB)){
        stab = matrix(data = NA, nrow = 2, ncol = 2)
        stab[,1] = c(sTAB[p,2],sTAB[p,"ALL"]) 
        stab[,2] = c(sum(sTAB[,2])-sTAB[p,2], sum(sTAB[,"ALL"])-sTAB[p,"ALL"])
        
        chi2 = chisq.test(stab)
        chi2$pvalue = c(pvalue, chi2$p.value)
        print(p)
        print(chi2)
        if (chi2$p.value == 0 | is.na(chi2$p.value)){
          star = c(star," ")
        }else if (chi2$p.value < 1e-200){
          star = c(star,"****")
        }else if (chi2$p.value < 1e-100){
          star = c(star,"***")
        }else if (chi2$p.value < 1e-20){
          star = c(star,"**")
        }else if (chi2$p.value < 1e-10){
          star = c(star,"*")
        }else {
          star = c(star,"ns")
        }
      }  
      sTAB = cbind(sTAB, pvalue, star)
      write.table(sTAB, paste0(save_path,f_motif,"/",f_arni,"_IES_stat.tab"), sep = "\t")
      
      # Faire les représentataions graphique
      data_plot = cbind(TAB[,1]/sum(TAB[,1])*100, TAB[,col]/TAB[,1]*100)
      colnames(data_plot) = c(names(FILTRE_ARNi)[1],f_arni)
      
      png(paste0(save_path,f_motif,"/",f_arni,"_IES_barplot.png"),width = 400, height = 500)
      barplot(t(data_plot),
              col = c("#353436","#1b98e0"),
              beside = TRUE,
              main = f_arni,
              ylab = "% de gènes")
      
      legend("topright",
             legend = c("Reference","ARNi"),
             fill = c("#353436","#1b98e0"),
             bty = "n")
      
      dev.off()
      
      if (sum(TAB[,col]) > 0 ){
        ## Diagame emplié ##
        data_plot = cbind(TAB[,1]/sum(TAB[,1])*100, TAB[,col]/sum(TAB[,col])*100)
        colnames(data_plot) = c(names(FILTRE_ARNi)[1],f_arni)
     
        
        png(paste0(save_path,f_motif,"/",f_arni,"_IES_barplot_pile.png"),width = 200, height = 500)
        barplot(data_plot,
                col = c("#353436","#1b98e0"),
                main = f_arni,
                ylab = "% de gènes")
        
        dev.off()
      }
      sink(file = NULL)
    }
    
    
    
  }
}

colnames(TAB) = name
write.table(TAB, paste0(save_path,"IES_FIMO2_IN_MAC_CDS_table.tab"), sep = "\t")

png(paste0(save_path,"/Legend_pile_IES.png"),width = 500, height = 500)
plot(NULL)
legend("topleft",
       legend = rev(rownames(data_plot)),
       fill = rev(c("#353436","#1b98e0")),
       bty = "n")
dev.off()



#### Récupérer les motifs identifiés ####
# library(seqinr)
# motifs = as.list(MOTIF$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_identifie_FIMO2_IM_MAC_CDS.fa" )
# 
# 
# 
# MOTIF_position = MOTIF[MOTIF$start < -50 & MOTIF$start > -80,]
# motifs = as.list(MOTIF_position$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF_position$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_position_FIMO2_IM_MAC_CDS.fa" )
# 
# 
# MOTIF_pos_plus = MOTIF_position[MOTIF_position$strand == "+",]
# motifs = as.list(MOTIF_pos_plus$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF_pos_plus$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_plus_FIMO2_IM_MAC_CDS.fa" )
# 
# MOTIF_pos_moins = MOTIF_position[MOTIF_position$strand == "-",]
# motifs = as.list(MOTIF_pos_moins$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF_pos_moins$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_moins_FIMO2_IM_MAC_CDS.fa" )
# 
# MOTIF_pos_moins_UP = MOTIF_pos_moins[which(is.element(MOTIF_pos_moins$sequence_name, UP$ID)),]
# motifs = as.list(MOTIF_pos_moins_UP$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF_pos_moins_UP$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_moins_UP_FIMO2_IM_MAC_CDS.fa" )
# 
# MOTIF_pos_plus_UP = MOTIF_pos_plus[which(is.element(MOTIF_pos_plus$sequence_name, UP$ID)),]
# motifs = as.list(MOTIF_pos_plus_UP$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF_pos_plus_UP$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_plus_UP_FIMO2_IM_MAC_CDS.fa" )
# 
# 
# MOTIF_pos_moins_DOWN = MOTIF_pos_moins[which(is.element(MOTIF_pos_moins$sequence_name, DOWN$ID)),]
# motifs = as.list(MOTIF_pos_moins_DOWN$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF_pos_moins_DOWN$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_moins_DOWN_FIMO2_IM_MAC_CDS.fa" )
# 
# MOTIF_pos_plus_DOWN = MOTIF_pos_plus[which(is.element(MOTIF_pos_plus$sequence_name, DOWN$ID)),]
# motifs = as.list(MOTIF_pos_plus_DOWN$matched_sequence)
# write.fasta(sequences = motifs, names = MOTIF_pos_plus_DOWN$sequence_name, 
#             file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_plus_DOWN_FIMO2_IM_MAC_CDS.fa" )

