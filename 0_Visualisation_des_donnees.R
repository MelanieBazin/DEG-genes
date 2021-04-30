source("0_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")
path = "./Graph/Normalisations/"

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
    if (i == "EXPRESSION"){type = "raw data"}else{type = i}
  CountBoxplot(tab, type)
  dev.off()
}
  
# Récupération des données de comptage normalisées DESeq2
tab = list.files("./DATA/DESeq2/")[grep("tout",list.files("./DATA/DESeq2/"))]
tab = tab[grep("normalisation", tab)]
tab = read.table(paste0("./DATA/DESeq2/", tab), header = T, sep = "\t")
  
png(paste0(path,"DESeq2_Boxplot.png"))
  CountBoxplot(tab, "DESeq2")
dev.off()

########################################################
# Partie 2 : Visualisation des données pour clustering #
########################################################
Type = c("RPM", "RPKM","DESeq2")
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
# i = names(rnai_list)[4]
# distance_methode = "Euclidean"
# method_utilisee = "kmeans"



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
    
    # Analyse en composante principale
    PCA_plot_generator(data_tab,colors = NULL, save_path = paste0(path,type,"_",i,"_"))
    
    # Nom_de_fichier = data_tab
    # rm(data_tab)
    # distance = distance_methode
    # nb_cluster = 4
    # method = method_utilisee
    # graph_type = c("heatmap","profils")
    
    
    # Analyse de clusering
    # for (distance_methode in c("Euclidean", "Correlation")){
    #   for (method_utilisee in c("kmeans", "HCL")){
    #     pdf(paste0(path, paste(type, i, distance_methode,method_utilisee,"Clustering",sep="_"),".pdf"))
    #     Clustering(Nom_de_fichier = data_tab,
    #                distance = distance_methode,
    #                nb_cluster = 4,
    #                method = method_utilisee,
    #                graph_type = c("heatmap","profils"))
    #     dev.off()
    #   }
    # }
  }
}





