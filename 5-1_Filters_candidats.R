####
# Definition of list of filter
# -> Analyse the deregulated genes
# -> Compare the different KD
####

AUTOGAMY = list()
for (R in unique(annotation$EXPRESSION_PROFIL)){
  tab = annotation[annotation$EXPRESSION_PROFIL == R,1]
  AUTOGAMY = c(AUTOGAMY, list(tab))
}
names(AUTOGAMY) = unique(annotation$EXPRESSION_PROFIL)


#### Filter candidates genes ####
# Loaclised the data
path = paste0("./Analyse/",date,"_DESeq2_analysis/")
path = paste0(path,
              list.files(path)[grep(paste0(condition,"_FC-", FC, "_pval-", pvalue),list.files(path))],"/")

# Definition of the time course that will be analysed
infodata = read.table(paste0(path, condition,"_infodata_collapse.tab"))
RNAi = unique(infodata$KnockDown)
RNAi = RNAi[RNAi != "ctrl"]

# Creation of a list contain the table used for the analysis 
DEtab_list = list()
for (R in RNAi){
  if (R == "CTIP"){
    DEtab_list = c(DEtab_list, CTIP_early = list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_EARLY.tab"), header=T,sep="\t",quote='')))
    DEtab_list = c(DEtab_list, CTIP_inter = list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_INTER.tab"), header=T,sep="\t",quote='')))
  }else {
    DEtab_list = c(DEtab_list, list(read.table(paste0(path,"DESeq/",R,"/DEgenes_",condition,"_LATE.tab"), header=T,sep="\t",quote='')))
  }
}
names(DEtab_list)[3:length(names(DEtab_list))]= RNAi[-1]

# UP deregulated genes  in PGM, KU80c &XRCC4 KD
UP_filter = list()
for (R in RNAi[-1]){
  tab = DEtab_list[[R]]
  tab = tab[tab$REGULATION == "Up-regulated", "ID"]
  UP_filter = c(UP_filter, list(tab))
}
rm(R, tab)
names(UP_filter) = paste("UP", RNAi[-1], sep = "_")
UP_filter = c(UP_filter,
              UP_PKX = list(intersect(UP_filter$UP_XRCC4, 
                                      intersect(UP_filter$UP_PGM,
                                                UP_filter$UP_KU80c))))
# DOWN regulated in PGM, KU80c &XRCC4 KD
DOWN_pkx_filter = list()
for (R in RNAi[-1]){
  tab = DEtab_list[[R]]
  tab = tab[tab$REGULATION == "Down-regulated", "ID"]
  DOWN_pkx_filter  = c(DOWN_pkx_filter , list(tab))
}
rm(R, tab)
names(DOWN_pkx_filter ) = paste("DOWN", RNAi[-1], sep = "_")
DOWN_pkx_filter  = c(DOWN_pkx_filter ,
                     DOWN_pkx_filter  = list(intersect(DOWN_pkx_filter$DOWN_XRCC4, 
                                                       intersect(DOWN_pkx_filter$DOWN_PGM,
                                                                 DOWN_pkx_filter$DOWN_KU80c))))

# DOWN deregulated genes in CTIP KD
DOWN_filter = list(D_CTIP_early = DEtab_list[["CTIP_early"]][DEtab_list[["CTIP_early"]]$REGULATION == "Down-regulated", "ID"],
                   D_CTIP_inter = DEtab_list[["CTIP_inter"]][DEtab_list[["CTIP_inter"]]$REGULATION == "Down-regulated", "ID"])

DOWN_filter = c(DOWN_filter, 
                DOWN_CTIP = list(unique(c(DOWN_filter[["D_CTIP_early"]], 
                                          DOWN_filter[["D_CTIP_inter"]]))))

# UP deregulated genes in CTIP KD
UP_C_filter = list(UP_CTIP_early = DEtab_list[["CTIP_early"]][DEtab_list[["CTIP_early"]]$REGULATION == "Up-regulated", "ID"],
                   UP_CTIP_inter = DEtab_list[["CTIP_inter"]][DEtab_list[["CTIP_inter"]]$REGULATION == "Up-regulated", "ID"])

UP_C_filter = c(UP_C_filter, 
                UP_CTIP = list(unique(c(UP_C_filter[["UP_CTIP_early"]], 
                                        UP_C_filter[["UP_CTIP_inter"]]))))

# Genes of interest
DE_genes = list(UP_PKX = UP_filter$UP_PKX,
                D_CTIP = DOWN_filter$DOWN_CTIP,
                UPpkx_Dc = intersect(UP_filter$UP_PKX, DOWN_filter$DOWN_CTIP)) 

candidats = intersect(DE_genes$UPpkx_Dc, AUTOGAMY$`Intermediate peak`)




