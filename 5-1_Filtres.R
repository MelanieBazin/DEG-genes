#♪ Récupéartion des données annotation
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
annotation = annotation[,c(1,3:5,13,6:11,2)]
rownames(annotation)=annotation$ID

AUTOGAMY = list(all_genes = annotation$ID,
                 autogamy = annotation$ID[annotation$EXPRESSION_PROFIL != "none"],
                 inter_peak = annotation$ID[annotation$EXPRESSION_PROFIL == "Intermediate peak"],
                 early_peak = annotation$ID[annotation$EXPRESSION_PROFIL == "Early peak"])

# Récupération des données tuboID
TurboPGM = read.table("./DATA/TurboID/2114003-Pgm-ProteinMeasurements.txt",header=T,sep="\t")
TurboPGML4 = read.table("./DATA/TurboID/2114003-PgmL4-ProteinMeasurements.txt",header=T,sep="\t",quote='')

TURBO = list(turbo_PGM = annotation$ID[which(is.element(annotation$PROTEIN_NAME, TurboPGM$PROTEIN_NAME))] ,
             turbo_PGML4 = annotation$ID[which(is.element(annotation$PROTEIN_NAME, TurboPGML4$PROTEIN_NAME))],
             turbo_OU = annotation$ID[which(is.element(annotation$PROTEIN_NAME, unique(TurboPGM$PROTEIN_NAME, TurboPGML4$PROTEIN_NAME)))])

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

up_pkx = intersect(UP_PKX$UP_XRCC4, intersect(UP_PKX$UP_PGM,UP_PKX$UP_KU80c))

# Les gènes NOT UP dérégulés en PGM, KU80c et/OU XRCC4
not_UP_PKX = list()
for (l in names(UP_PKX)){
  up_genes = UP_PKX[[l]]
  not_UP_PKX = c(not_UP_PKX, list(AUTOGAMY$all_genes[!is.element(AUTOGAMY$all_genes, up_genes)]))
}
names(not_UP_PKX) = paste("not",names(UP_PKX), sep = "_")



# Les gènes DOWN dérégulés en CTIP early ou inter
CTIP = rbind(TAB$CTIP_early,TAB$CTIP_inter)

stdCTIP = list(DOWN_CTIP = unique(CTIP[CTIP$REGULATION == "Down-regulated", "ID"]),
               UP_ALL = up_pkx,
               DOWN_UP = intersect(unique(CTIP[CTIP$REGULATION == "Down-regulated", "ID"]), up_pkx))


not_stdCTIP = list()
for (l in names(stdCTIP)){
  up_genes = stdCTIP[[l]]
  not_stdCTIP = c(not_stdCTIP, list(AUTOGAMY$all_genes[!is.element(AUTOGAMY$all_genes, up_genes)]))
}
names(not_stdCTIP) = paste("not",names(stdCTIP), sep = "_")

Crossinglist <- function (list1, list2){
  LIST = list()
  list_name = c()
  for (l in names(list1)){
    for (m in names(list2)){
      LIST = c(LIST, list(intersect(list1[[l]],list2[[m]])))
      list_name = c(list_name, paste(names(list1[l]),names(list2[m]), sep = "x"))
    }
  }
  names(LIST)=list_name
  
  return(LIST)
}

