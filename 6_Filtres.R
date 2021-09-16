options(stringsAsFactors = FALSE)
library(stringr) 

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"

#### Initialisation des filtres ####
infogenes = read.table(paste0(path,"Resumer_DEgenes.tab"), header = T)
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")


Tous_les_genes = unique(na.omit(infogenes$ID))
Genes_codant = unique(na.omit(infogenes$ID[infogenes$TYPE == "Coding Gene"]))
UP_PGM = unique(na.omit(infogenes$ID[infogenes$PGM_LATE_REGULATION == "Up-regulated"]))
UP_KU80c = unique(na.omit(infogenes$ID[infogenes$KU80c_LATE_REGULATION == "Up-regulated"]))
UP_XRCC4 = unique(na.omit(infogenes$ID[infogenes$XRCC4_LATE_REGULATION == "Up-regulated"]))
UP_EZL1 = unique(na.omit(infogenes$ID[infogenes$EZL1_LATE_REGULATION == "Up-regulated"]))
DOWN_CTIP = unique(na.omit(infogenes$ID[infogenes$CTIP_EARLY_REGULATION == "Down-regulated" | infogenes$CTIP_INTER_REGULATION == "Down-regulated"]))
UP_PGM_KU80c = unique(intersect(UP_PGM, UP_KU80c))
UP_PKX = unique(intersect(UP_PGM_KU80c, UP_XRCC4))
UP_PKXE= unique(intersect(UP_PKX, UP_EZL1))

# selected_genes = merge(annotation,infogenes[,c(1,grep("REGULATION", colnames(infogenes)))], by = "ID")
# selected_genes = selected_genes[which(is.element(selected_genes$ID, UP_PKXE)),]
# write.table(selected_genes,paste0(path, "Genes_UP_PKXE.tab") ,row.names = F, sep = "\t")


UP_PKX_DOWN_CTIP = unique(intersect(UP_PKX, DOWN_CTIP))
UP_DOWN = unique(intersect(UP_PKXE, DOWN_CTIP))
DE_autogamie = unique(na.omit(infogenes$ID[infogenes$EXPRESSION_PROFIL != "none"]))
Pic_intermediaire = unique(na.omit(infogenes$ID[infogenes$EXPRESSION_PROFIL == "Intermediate peak"]))
UP_PKX_inter = unique(intersect(Pic_intermediaire,UP_PKX))
UP_PKXE_inter = unique(intersect(Pic_intermediaire,UP_PKXE))
UP_PKX_DOWN_inter = unique(intersect(Pic_intermediaire,UP_PKX_DOWN_CTIP))
UP_DOWN_inter = unique(intersect(Pic_intermediaire,UP_DOWN))
Pic_early = unique(na.omit(infogenes$ID[infogenes$EXPRESSION_PROFIL == "Early peak"]))
Early_rep = unique(na.omit(infogenes$ID[infogenes$EXPRESSION_PROFIL == "Early repression"]))

DOWN_PGM = unique(na.omit(infogenes$ID[infogenes$PGM_LATE_REGULATION == "Down-regulated"]))
DOWN_KU80c = unique(na.omit(infogenes$ID[infogenes$KU80c_LATE_REGULATION == "Down-regulated"]))
DOWN_XRCC4 = unique(na.omit(infogenes$ID[infogenes$XRCC4_LATE_REGULATION == "Down-regulated"]))
DOWN_EZL1 = unique(na.omit(infogenes$ID[infogenes$EZL1_LATE_REGULATION == "Down-regulated"]))
UP_CTIP = unique(na.omit(infogenes$ID[infogenes$CTIP_EARLY_REGULATION == "Up-regulated" | infogenes$CTIP_INTER_REGULATION == "Up-regulated"]))
DOWN_PGM_KU80c = unique(intersect(DOWN_PGM, DOWN_KU80c))
DOWN_PKX = unique(intersect(DOWN_PGM_KU80c, DOWN_XRCC4))
DOWN_PKXE= unique(intersect(DOWN_PKX, DOWN_EZL1))
DOWN_PKX_UP_CTIP = unique(intersect(DOWN_PKX, UP_CTIP))
DOWN_UP = unique(intersect(DOWN_PKXE, UP_CTIP))

DOWN_PKX_inter = unique(intersect(Pic_intermediaire,DOWN_PKX))
DOWN_PKXE_inter = unique(intersect(Pic_intermediaire,DOWN_PKXE))
DOWN_PKX_UP_inter = unique(intersect(Pic_intermediaire,DOWN_PKX_UP_CTIP))
DOWN_UP_inter = unique(intersect(Pic_intermediaire,DOWN_UP))


FILTRES = list(Tous_les_genes,Genes_codant,
               UP_PGM, UP_KU80c, UP_XRCC4, UP_EZL1, DOWN_CTIP,
               UP_PGM_KU80c, UP_PKX, UP_PKXE, 
               UP_PKX_DOWN_CTIP, UP_DOWN ,
               DE_autogamie, Pic_intermediaire,
               UP_PKX_inter, UP_PKXE_inter,
               UP_PKX_DOWN_inter, UP_DOWN_inter)

names(FILTRES) = c("Tous_les_genes","Genes_codant",
                   "UP_PGM", "UP_KU80c", "UP_XRCC4", "UP_EZL1", "DOWN_CTIP",
                   "UP_PGM_KU80c", "UP_PKX", "UP_PKXE", 
                   "UP_PKX_DOWN_CTIP", "UP_DOWN" ,
                   "DE_autogamie", "Pic_intermediaire",
                   "UP_PKX_inter", "UP_PKXE_inter",
                   "UP_PKX_DOWN_inter", "UP_DOWN_inter")

ANTIFILTRES = list(Tous_les_genes,Genes_codant,
                   DOWN_PGM, DOWN_KU80c, DOWN_XRCC4, DOWN_EZL1, UP_CTIP,
                   DOWN_PGM_KU80c, DOWN_PKX, DOWN_PKXE, 
                   DOWN_PKX_UP_CTIP, DOWN_UP ,
                   DE_autogamie, Pic_intermediaire,
                   DOWN_PKX_inter, DOWN_PKXE_inter,
                   DOWN_PKX_UP_inter, DOWN_UP_inter)

names(ANTIFILTRES) = c("Tous_les_genes","Genes_codant",
                       "DOWN_PGM", "DOWN_KU80c", "DOWN_XRCC4", "DOWN_EZL1", "UP_CTIP",
                       "DOWN_PGM_KU80c", "DOWN_PKX", "DOWN_PKXE", 
                       "DOWN_PKX_UP_CTIP", "DOWN_UP" ,
                       "DE_autogamie", "Pic_intermediaire",
                       "DOWN_PKX_inter", "DOWN_PKXE_inter",
                       "DOWN_PKX_UP_inter", "DOWN_UP_inter")