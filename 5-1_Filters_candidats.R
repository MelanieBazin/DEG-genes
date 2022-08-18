####
# Definition of list of filter ####
# -> Analyse the deregulated genes
# -> Compare the different KD
####

GENES_ID = list(all_genes = annotation$ID,
                autogamy = annotation$ID[annotation$EXPRESSION_PROFIL != "none"],
                inter_peak = annotation$ID[annotation$EXPRESSION_PROFIL == "Intermediate peak"])

#### Filter candidates genes ####
# UP deregulated genes
UP_filter = list()
for (R in RNAi[-1]){
  tab = DEtab_list[[R]]
  tab = tab[tab$REGULATION == "Up-regulated", "ID"]
  UP_filter = c(UP_filter, list(tab))
}
rm(R)
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
                DOWN_CTIP = list(unique(DOWN_filter[["CTIP_early"]], 
                                        DOWN_filter[["CTIP_inter"]])))

# Genes of interest
DE_genes = list(UP_PKX = UP_filter$UP_PKX,
                DOWN_CTIP = DOWN_filter$DOWN_CTIP,
                UPpkx_DOWN_c = intersect(upPKX, DE_genes$DOWN_CTIP)) 

candidats = intersect(upPKX_downC, GENES_ID$inter_peak)

#### Cross the UP_filter with the DOWN_filter
UP_PKX.DOWN_C = Crossinglist(UP_filter,DOWN_filter)



