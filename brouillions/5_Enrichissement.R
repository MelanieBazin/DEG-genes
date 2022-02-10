source("5_Filtres_liste_genes.R")
path = "./Analyse/2021-07-09_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/test/"

nb_sim_corr = 100

all = nrow(ALL)

# Hestimation des enrichissement par la loi hypèregéométrique
for (f1 in names(FILTRE_ARNi)){
  gene_list = FILTRE_ARNi[[f1]]
  
  pVal_tab = Enrichissement(annotation, gene_list, nb_sim_corr)
  write.table(pVal_tab, paste0(path, f1, "_adj_pVal.tab"), sep="\t",row.names=T,quote=F)
  
  
  
}

for (f2 in names(FILTRE_MOTIF)){
  gene_list = FILTRE_MOTIF[[f2]]
  
  pVal_tab = Enrichissement(annotation, gene_list, nb_sim_corr)
  write.table(pVal_tab, paste0(path, f2, "_adj_pVal.tab"), sep="\t",row.names=T,quote=F)
  
}

# Hestimation des enrichissement par chi2