
all = 41533
PGM = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","PGM","/NoFilter/DEgenes_tout_LATE_NoFilter.tab"), sep = "\t", h=T )
KU80c = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","KU80c","/NoFilter/DEgenes_tout_LATE_NoFilter.tab"), sep = "\t", h=T)
XRCC4 = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","XRCC4","/NoFilter/DEgenes_tout_LATE_NoFilter.tab"), sep = "\t", h=T )

CTIP_e = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","CTIP","/NoFilter/DEgenes_tout_EARLY_NoFilter.tab"), sep = "\t", h=T )
CTIP_i = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/DESeq/","CTIP","/NoFilter/DEgenes_tout_INTER_NoFilter.tab"), sep = "\t", h=T )


tab = table(PGM$REGULATION)
print("PGM")
print(rbind(tab,tab/all*100))

tab = table(KU80c$REGULATION)
print("KU80c")
print(rbind(tab,tab/all*100))

tab = table(XRCC4$REGULATION)
print("XRCC4")
print(rbind(tab,tab/all*100))


tab = table(CTIP_e$REGULATION)
print("CTIP_e")
print(rbind(tab,tab/all*100))

tab = table(CTIP_i$REGULATION)
print("CTIP_i")
print(rbind(tab,tab/all*100))
