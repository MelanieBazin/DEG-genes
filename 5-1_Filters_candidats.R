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
              list.files(path)[grep(condition,list.files(path))],"/")

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

# UP deregulated genes
UP_filter = list()
for (R in RNAi[-1]){
  tab = DEtab_list[[R]]
  tab = tab[tab$REGULATION == "Up-regulated", "ID"]
  UP_filter = c(UP_filter, list(tab))
}
rm(R, tab)
names(UP_filter) = paste("UP", RNAi[-1], sep = "_")
UP_filter = c(UP_filter, 
              UP_PK =list(intersect(UP_filter$UP_PGM,
                                    UP_filter$UP_KU80c)),
              UP_PX =list(intersect(UP_filter$UP_PGM,
                                    UP_filter$UP_XRCC4)),
              UP_XK =list(intersect(UP_filter$UP_XRCC4,
                                    UP_filter$UP_KU80c)),
              UP_PKX = list(intersect(UP_filter$UP_XRCC4, 
                                      intersect(UP_filter$UP_PGM,
                                                UP_filter$UP_KU80c))))

# DOWN deregulated genes
DOWN_filter = list(CTIP_early = DEtab_list[["CTIP_early"]][DEtab_list[["CTIP_early"]]$REGULATION == "Down-regulated", "ID"],
                   CTIP_inter = DEtab_list[["CTIP_inter"]][DEtab_list[["CTIP_inter"]]$REGULATION == "Down-regulated", "ID"])

DOWN_filter = c(DOWN_filter, 
                DOWN_CTIP = list(unique(c(DOWN_filter[["CTIP_early"]], 
                                          DOWN_filter[["CTIP_inter"]]))))

# Genes of interest
DE_genes = list(UP_PKX = UP_filter$UP_PKX,
                D_CTIP = DOWN_filter$DOWN_CTIP,
                UPpkx_Dc = intersect(UP_filter$UP_PKX, DOWN_filter$DOWN_CTIP)) 

candidats = intersect(DE_genes$UPpkx_Dc, AUTOGAMY$`Intermediate peak`)




