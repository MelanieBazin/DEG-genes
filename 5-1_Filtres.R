all_genes = annotation$ID
inter_genes = annotation$ID[annotation$EXPRESSION_PROFIL == "Intermediate peak"]
early_genes = annotation$ID[annotation$EXPRESSION_PROFIL == "Early peak"]
auto_genes = annotation$ID[annotation$EXPRESSION_PROFIL != "none"]

#### Ouverture tableau de données des gènes DEG ####
path = paste0("./Analyse/",file_name, "/", condition,"/")
TAB = list()
for (R in RNAi){
  if (R == "CTIP"){
    TAB = c(TAB, CTIP_early = list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_EARLY_NoFilter.tab"), header=T,sep="\t",quote='')))
    TAB = c(TAB, CTIP_inter = list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_INTER_NoFilter.tab"), header=T,sep="\t",quote='')))
  }else {
    TAB = c(TAB, list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_LATE_NoFilter.tab"), header=T,sep="\t",quote='')))
  }
}
names(TAB)[3:length(names(TAB))]= RNAi[-1]

#### génération des filtres ####
# Les gènes UP dérégulés en PGM, KU80c et/OU XRCC4
UP_PKX = list()
for (R in RNAi[-1]){
  tab = TAB[[R]]
  tab = tab[tab$REGULATION == "Up-regulated", "ID"]
  UP_PKX = c(UP_PKX, list(tab))
}
names(UP_PKX) = paste("UP", RNAi[-1], sep = "_")

up_pkx = intersect(UP_PKX[["UP_XRCC4"]], intersect(UP_PKX[["UP_PGM"]],UP_PKX[["UP_KU80c"]]))

# Les gènes NOT UP dérégulés en PGM, KU80c et/OU XRCC4
not_UP_PKX = list()
for (l in names(UP_PKX)){
  up_genes = UP_PKX[[l]]
  not_UP_PKX = c(not_UP_PKX, list(setdiff(all_genes, up_genes)))
}
names(not_UP_PKX) = paste("not",names(UP_PKX), sep = "_")

# Les gènes DOWN dérégulés en PGM, KU80c et/OU XRCC4
DOWN_PKX = list()
for (R in RNAi[-1]){
  tab = TAB[[R]]
  tab = tab[tab$REGULATION == "Down-regulated", "ID"]
  DOWN_PKX = c(DOWN_PKX, list(tab))
}
names(DOWN_PKX) = paste("DOWN", RNAi[-1], sep = "_")

# Les gènes NOT dérégulés en PGM, KU80c et/OU XRCC4
not_DE_PKX = list()
for (l in names(UP_PKX)){
  up_genes = UP_PKX[[l]]
  down_genes = DOWN_PKX[[l]]
  not_DE_PKX = c(not_DE_PKX, list(setdiff(all_genes, unique(down_genes,up_genes))))
}
names(not_DE_PKX) = paste("not_DE",RNAi[-1], sep = "_")

# Les gènes DOWN dérégulés en CTIP early ou inter
CTIP = rbind(TAB[["CTIP_early"]],TAB[["CTIP_inter"]])

DOWN_C = unique(CTIP[CTIP$REGULATION == "Down-regulated", "ID"])
up_pkx_down_c = intersect(up_pkx, DOWN_C)

# Les gènes UP dérégulés en CTIP
UP_C = unique(CTIP[CTIP$REGULATION == "Up-regulated", "ID"])

# Les gènes NOT DE dérégulés en CTIP
not_DE_C = setdiff(all_genes, unique(DOWN_C,UP_C))

