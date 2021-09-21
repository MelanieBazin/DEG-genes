

ALL = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/Resumer_DEgenes.tab"), sep = "\t", h=T )
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
  UP_DOWN = ALL[which(is.element(ALL$ID,intersect(DOWN$ID, CTIP_UP$ID))),]
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
# f_arni = names(FILTRE_ARNi)[1]


t = table(ALL$EXPRESSION_PROFIL)
TAB= array(dim = c(length(t),length(names(FILTRE_MOTIF))*length(names(FILTRE_ARNi))))
rownames(TAB) = names(t)


col = 0
name = c()
for(f_motif in names(FILTRE_MOTIF)){
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

  }
}

colnames(TAB) = name
write.table(TAB, "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/Expression_profil_FIMO2_IN_MAC_CDS_table.tab", sep = "\t")

#### Regarder le nombre de gènes IES + parmis les différnets cathégories ###

TAB= array(dim = c(2,length(names(FILTRE_MOTIF))*length(names(FILTRE_ARNi))))
rownames(TAB) = c("SANS_IES","AVEC_IES")


col = 0
name = c()
for(f_motif in names(FILTRE_MOTIF)){
  tab1 = FILTRE_MOTIF[[f_motif]]
  
  for (f_arni in names(FILTRE_ARNi)){
    tab2 = FILTRE_ARNi[[f_arni]]
    tab = tab2[which(is.element(tab2$ID, tab1$sequence_name)),]
    
    col = col+1
    
    TAB["SANS_IES",col] = sum(grepl(0, tab$NB_IES))
    TAB["AVEC_IES",col] = sum(!grepl(0, tab$NB_IES))

    name = c(name, paste(f_arni,f_motif, sep = "_"))
    
  }
}

colnames(TAB) = name
write.table(TAB, "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/IES_FIMO2_IN_MAC_CDS_table.tab", sep = "\t")



#### Récupérer les motifs identifiés ####
library(seqinr)
motifs = as.list(MOTIF$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_identifie_FIMO2_IM_MAC_CDS.fa" )



MOTIF_position = MOTIF[MOTIF$start < -50 & MOTIF$start > -80,]
motifs = as.list(MOTIF_position$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF_position$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_position_FIMO2_IM_MAC_CDS.fa" )


MOTIF_pos_plus = MOTIF_position[MOTIF_position$strand == "+",]
motifs = as.list(MOTIF_pos_plus$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF_pos_plus$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_plus_FIMO2_IM_MAC_CDS.fa" )

MOTIF_pos_moins = MOTIF_position[MOTIF_position$strand == "-",]
motifs = as.list(MOTIF_pos_moins$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF_pos_moins$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_moins_FIMO2_IM_MAC_CDS.fa" )

MOTIF_pos_moins_UP = MOTIF_pos_moins[which(is.element(MOTIF_pos_moins$sequence_name, UP$ID)),]
motifs = as.list(MOTIF_pos_moins_UP$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF_pos_moins_UP$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_moins_UP_FIMO2_IM_MAC_CDS.fa" )

MOTIF_pos_plus_UP = MOTIF_pos_plus[which(is.element(MOTIF_pos_plus$sequence_name, UP$ID)),]
motifs = as.list(MOTIF_pos_plus_UP$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF_pos_plus_UP$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_plus_UP_FIMO2_IM_MAC_CDS.fa" )


MOTIF_pos_moins_DOWN = MOTIF_pos_moins[which(is.element(MOTIF_pos_moins$sequence_name, DOWN$ID)),]
motifs = as.list(MOTIF_pos_moins_DOWN$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF_pos_moins_DOWN$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_moins_DOWN_FIMO2_IM_MAC_CDS.fa" )

MOTIF_pos_plus_DOWN = MOTIF_pos_plus[which(is.element(MOTIF_pos_plus$sequence_name, DOWN$ID)),]
motifs = as.list(MOTIF_pos_plus_DOWN$matched_sequence)
write.fasta(sequences = motifs, names = MOTIF_pos_plus_DOWN$sequence_name, 
            file.out =  "Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/Motif_pos_plus_DOWN_FIMO2_IM_MAC_CDS.fa" )

