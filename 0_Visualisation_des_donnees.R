source("0_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")
path = "./Graph/Normalisations/Color"

#######################################################
# Partie 1 : Comparaison des méthodes de normalisaion #
#######################################################
# Comptage

Type = c("EXPRESSION","RPM", "RPKM")
for (i in Type){
  png(paste0(path,i,"_Boxplot.png"))
    tab = ConcatTab(i)
    rownames(tab)=tab$ID
    tab = tab[,-1]
    tab = tab[,c(12:ncol(tab),1:11)]
    if (i == "EXPRESSION"){type = "raw data"}else{type = i}
    
  CountBoxplot(tab, type, color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
  dev.off()
}

# Récupération des données de comptage normalisées DESeq2
tab = list.files("./DATA/DESeq2/")[grep("tout",list.files("./DATA/DESeq2/"))]
tab = tab[grep("normalisation", tab)]
tab = read.table(paste0("./DATA/DESeq2/", tab), header = T, sep = "\t")
tab = tab[,c(12:ncol(tab),1:11)]

png(paste0(path,"DESeq2_Boxplot.png"))
  CountBoxplot(tab, "DESeq2", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
dev.off()

########################################################
# Partie 2 : Visualisation des données pour clustering #
########################################################
Type = c("RPKM", "RPM","DESeq2")
# RNAi à analyser ensemble
tout = sub(paste0("_expression_table_RPKM.tab"),"",list.files("./DATA/RPKM/"))
rnai_list = list(
  sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
  sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
  controles = tout[which(is.element(tout,c("ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))],
  XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
  XRCC4ctrl2020 = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","CTIP_CTRL")))],
  XRCC4tousctrl = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","ND7", "ICL7", "CTIP_CTRL")))],
  XRCC4xseq2014 = tout[which(is.element(tout,c("XRCC4","XRCC4_CTRL","ICL7","KU80c","ND7","PGM" )))],
  CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
  CTIPctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL","XRCC4_CTRL")))],
  CTIPtousctrl = tout[which(is.element(tout,c("CTIP","ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))],
  tout = tout
)


# type = Type[1]
# i = names(rnai_list)[1]
# distance_methode = "Pearson"
# method_utilisee = "HCL"
# 
# Nom_de_fichier = data_tab
# # rm(data_tab)
# nb_cluster = 4
# graph_type = c("profils")


pdf(paste0(path,"4Cluster/Tout.pdf"))
for (type in Type){
  for (i in names(rnai_list)){
    # Création du tableau de donnée à analyser ensemble
    if (type == "DESeq2"){
      data_tab = read.table(paste0("./DATA/DESeq2/",i,"_",paste(rnai_list[[i]],collapse = "-"),"_normalisation_DESeq2.tab"), header = T, sep="\t")
    }else{
      data_tab = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])
    }
    
    # Passage de la colonne des ID en rowname
    if (colnames(data_tab)[1]=="ID"){
      row.names(data_tab)=data_tab$ID
      data_tab = data_tab[,-1]
    }
    
    print(paste(type,i, "- Commencé"))
    
    # Analyse en composante principale
    # print(paste(type,i, "-----> début ACP"))
    # resExp=PCA_plot_generator(data_tab,colors = NULL, save_path = paste0(path,"/ACP/",type,"_",i,"_"), main = paste0("ACP ",i," (",type,")"))
    # 
    # resHCPC = HCPC(resExp, graph = FALSE)
    # 
    # png(paste0(path,"/ACP/",type,"_",i,"_PCA_ClusteringHCPC_dendog.png"))
    # p = fviz_dend(resHCPC, 
    #               cex = 0.5,                     # Taille du text
    #               palette = "jco",               # Palette de couleur ?ggpubr::ggpar
    #               rect = TRUE, rect_fill = TRUE, # Rectangle autour des groupes
    #               rect_border = "jco",           # Couleur du rectangle
    #               main = paste0("ACP Dendrogram ",i," (",type,")"), 
    #               labels_track_height = 2000)   # Augment l'espace pour le texte
    # print(p)
    # dev.off()
    # 
    # png(paste0(path,"/ACP/",type,"_",i,"_PCA_ClusteringHCPC_clust.png"))
    # p=fviz_cluster(resHCPC,
    #                repel = TRUE,            # Evite le chevauchement des textes
    #                show.clust.cent = F,     # Montre le centre des clusters
    #                palette = "jco",         # Palette de couleurs, voir ?ggpubr::ggpar
    #                ggtheme = theme_bw(),
    #                main = paste0("ACP Factor map ",i," (",type,")"))
    # print(p)
    # dev.off()
    # print(paste(type,i, "-----> fin ACP"))
    
    
    # Analyse de clusering
    for (distance in c("Pearson", "Spearman")){
      for (method in c("kmeans", "HCL")){
        print(paste(type,i, "----->",distance, method))
        
        Clustering(data_tab = data_tab,
                   distance = distance,
                   nb_cluster = 4,
                   method = method,
                   graph_type = c("profils"),
                   titre = paste(type, i))
       
      }
    }
    
    print(paste(type,i," - Fait"))
  }
}
dev.off()
print("-FIN-")



