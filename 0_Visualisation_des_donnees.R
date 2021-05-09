#install.packages("stringr") 
# install.packages("pheatmap")


source("0_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")
source("0_Cluster.R")
library("stringr")  
library("pheatmap")
path = "./Graph/Normalisations/Color/"

#######################################################
# Partie 1 : Comparaison des méthodes de normalisaion #
#######################################################
# Comptage

# Type = c("EXPRESSION","RPM", "RPKM")
# for (i in Type){
#   png(paste0(path,i,"_Boxplot.png"), width = 500, height = 500)
#     tab = ConcatTab(i)
#     rownames(tab)=tab$ID
#     tab = tab[,-1]
#     tab = tab[,c(12:ncol(tab),1:11)]
#     if (i == "EXPRESSION"){type = "raw data"}else{type = i}
#     
#   CountBoxplot(tab, type, color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
#   dev.off()
# }
# 
# # Récupération des données de comptage normalisées DESeq2
# tab = list.files("./DATA/DESeq2/")[grep("tout",list.files("./DATA/DESeq2/"))]
# tab = tab[grep("normalisation", tab)]
# tab = read.table(paste0("./DATA/DESeq2/", tab), header = T, sep = "\t")
# tab = tab[,c(12:ncol(tab),1:11)]
# 
# png(paste0(path,"DESeq2_Boxplot.png"),  width = 500, height = 500)
#   CountBoxplot(tab, "DESeq2", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
# dev.off()

########################################################
# Partie 2 : Visualisation des données pour clustering #
########################################################
# Type = c("RPKM", "RPM","DESeq2")
# # RNAi à analyser ensemble
# tout = sub(paste0("_expression_table_RPKM.tab"),"",list.files("./DATA/RPKM/"))
# rnai_list = list(
#   sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
#   sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
#   controles = tout[which(is.element(tout,c("ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))],
#   XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
#   XRCC4ctrl2020 = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","CTIP_CTRL")))],
#   XRCC4tousctrl = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","ND7", "ICL7", "CTIP_CTRL")))],
#   XRCC4xseq2014 = tout[which(is.element(tout,c("XRCC4","XRCC4_CTRL","ICL7","KU80c","ND7","PGM" )))],
#   CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
#   CTIPctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL","XRCC4_CTRL")))],
#   CTIPtousctrl = tout[which(is.element(tout,c("CTIP","ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))],
#   tout = tout
# )


# pdf(paste0(path,"4Cluster/Tout.pdf"))
# for (type in Type){
#   for (i in names(rnai_list)){
#     # Création du tableau de donnée à analyser ensemble
#     if (type == "DESeq2"){
#       data_tab = read.table(paste0("./DATA/DESeq2/",i,"_",paste(rnai_list[[i]],collapse = "-"),"_normalisation_DESeq2.tab"), header = T, sep="\t")
#     }else{
#       data_tab = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])
#     }
#     
#     # Passage de la colonne des ID en rowname
#     if (colnames(data_tab)[1]=="ID"){
#       row.names(data_tab)=data_tab$ID
#       data_tab = data_tab[,-1]
#     }
#     
#     print(paste(type,i, "- Commencé"))
#     
#     # Analyse en composante principale
#     # print(paste(type,i, "-----> début ACP"))
#     # resExp=PCA_plot_generator(data_tab,colors = NULL, save_path = paste0(path,"/ACP/",type,"_",i,"_"), main = paste0("ACP ",i," (",type,")"))
#     # 
#     # resHCPC = HCPC(resExp, graph = FALSE)
#     # 
#     # png(paste0(path,"/ACP/",type,"_",i,"_PCA_ClusteringHCPC_dendog.png"))
#     # p = fviz_dend(resHCPC, 
#     #               cex = 0.5,                     # Taille du text
#     #               palette = "jco",               # Palette de couleur ?ggpubr::ggpar
#     #               rect = TRUE, rect_fill = TRUE, # Rectangle autour des groupes
#     #               rect_border = "jco",           # Couleur du rectangle
#     #               main = paste0("ACP Dendrogram ",i," (",type,")"), 
#     #               labels_track_height = 2000)   # Augment l'espace pour le texte
#     # print(p)
#     # dev.off()
#     # 
#     # png(paste0(path,"/ACP/",type,"_",i,"_PCA_ClusteringHCPC_clust.png"))
#     # p=fviz_cluster(resHCPC,
#     #                repel = TRUE,            # Evite le chevauchement des textes
#     #                show.clust.cent = F,     # Montre le centre des clusters
#     #                palette = "jco",         # Palette de couleurs, voir ?ggpubr::ggpar
#     #                ggtheme = theme_bw(),
#     #                main = paste0("ACP Factor map ",i," (",type,")"))
#     # print(p)
#     # dev.off()
#     # print(paste(type,i, "-----> fin ACP"))
#     
#     
#     # Analyse de clusering
#     # Choisir le mode de calcule des distances
        # if (distance == "Pearson"){
        #   matDist = as.matrix(cor(data_tab))
        #   pheatmap(matDist, main = paste("Pheatmap Pearson", type, i))
        #   matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
        #   
        # }else if (distance == "Spearman"){
        #   matDist = as.matrix(cor(data_tab,method="spearman"))
        #   pheatmap(matDist, main = paste("Pheatmap Spearman", type, i))
        #   matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
        # }
        # 
        # 
        # 
        # for (method in c("kmeans", "HCL")){
        #   print(paste(type,i, "----->",distance, method))
        #   Clustering(matDist = matDist,
        #              nb_cluster = 5,
        #              method = method,
        #              titre = paste(type,i),
        #              colors = color)
#        
#       }
#     }
#     
#     print(paste(type,i," - Fait"))
#   }
# }
# dev.off()
# print("-FIN-")

#######################################################
# Partie 3 : Analyse ciblé et coloration de graphique #
#######################################################
# Utilisation uniquement des normalisations DESeq2
type = "DESeq2"
# RNAi à analyser ensemble
tout = sub("_expression_table_RPKM.tab","",list.files("./DATA/RPKM/"))
rnai_list = list(
  tout = tout,
  sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
  controles_2014 = tout[which(is.element(tout,c("ND7", "ICL7")))],
  sequencage_2014bis = tout[which(is.element(tout,c("KU80c","ND7","PGM" )))],
  sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
  #controles_2020 = tout[which(is.element(tout,c( "CTIP_CTRL","XRCC4_CTRL")))],
  XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
  CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
  CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL", "XRCC4_CTRL")))]
)

i = "CTIPseulctrl2020"
for (i in names(rnai_list)){
  
  # Création du tableau de donnée à analyser ensemble
  data_tab = read.table(paste0("./DATA/DESeq2/",i,"_normalisation_DESeq2.tab"), header = T, sep="\t")

  # Passage de la colonne des ID en rowname
  if (colnames(data_tab)[1]=="ID"){
    row.names(data_tab)=data_tab$ID
    data_tab = data_tab[,-1]
  }

  # Créaction du vecteur de couleur par cluster
  color = colnames(data_tab)
  for (j in rnai_list[[i]]){
    color[grep(j, color)]=cluster_color[[j]]
  }
  
  print(paste(i, "- Commencé"))
  
  # Changer le nom des colonnes controles
  colnames(data_tab) = str_replace_all(colnames(data_tab),"ND7","ND7_K")
  colnames(data_tab) = str_replace_all(colnames(data_tab),"CTIP_CTRL","ND7_C")
  colnames(data_tab) = str_replace_all(colnames(data_tab),"XRCC4_CTRL","ND7_X")
  
  # Analyse en composante principale
  # print(paste(i, "-----> début ACP"))
  # PCA_plot_generator(data_tab,colors = color, 
  #                           save_path = paste0(path,"/ACP/DESeq2_",i,"_"),
  #                           main = paste0("ACP ",i," (",type,")"))


  # Analyse de clusering
  for (distance in c("Pearson", "Spearman")){
    png(paste0(path,"4Cluster/",i,"_Matrice_",distance,".png"))
    # Choisir le mode de calcule des distances
    if (distance == "Pearson"){
      matDist = as.matrix(cor(data_tab))
      pheatmap(matDist, main = paste("Pheatmap Pearson", type, i))
      matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))

    }else if (distance == "Spearman"){
      matDist = as.matrix(cor(data_tab,method="spearman"))
      pheatmap(matDist, main = paste("Pheatmap Spearman", type, i))
      matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
    }
    dev.off()


    for (method in c("kmeans", "HCL")){
      print(paste(type,i, "----->",distance, method))
      png(paste0(path,"4Cluster/",i,"_Cluster_",method,"_",distance,".png"))
      Clustering(matDist = matDist,
                 nb_cluster = 5,
                 method = method,
                 titre = paste(type,i))
      dev.off()
    }
  }

}

