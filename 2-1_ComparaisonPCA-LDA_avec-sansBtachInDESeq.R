source("0_Cluster.R")
source("2_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")



path = "./Analyse/DESeq2_test02/"
# Utilisation uniquement des normalisations DESeq2
type = "DESeq2"
# RNAi à analyser ensemble
tout = sub("_expression_table_RPKM.tab","",list.files("./DATA/RPKM/"))
rnai_list = list(
  tout = tout,
  controles = tout[which(is.element(tout,c("ND7","ICL7","CTIP_CTRL","XRCC4_CTRL" )))]
  # ,
  # sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
  # controles_2014 = tout[which(is.element(tout,c("ND7", "ICL7")))],
  # sequencage_2014bis = tout[which(is.element(tout,c("KU80c","ND7","PGM" )))],
  # sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
  # #controles_2020 = tout[which(is.element(tout,c( "CTIP_CTRL","XRCC4_CTRL")))],
  # XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
  # CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
  # CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL", "XRCC4_CTRL")))]
)

for (i in names(rnai_list)){

  ##### Création du tableau de donnée à analyser ensemble ####
  #Ouverture des fichiers et création de l'objet countdata
  countdata = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])
  row.names(countdata)=countdata$ID
  countdata=countdata[,-1]
  
  # Boxplot des comptages non-normalisés / longeur de gènes
  png(paste0(path,i,"_Row_Boxplot.png"))
    CountBoxplot(apply(countdata,1,), "DESeq2") ### Ajouter une divison par la longeur des gènes pour le box plot
  dev.off()
  
  # Création du tableau avec les info des colonnes
  infodata=CreatInfoData3(conditions = rnai_list[[i]])  # La verison 3 ajoute une colonne batch avec les année de séquancage
  
  # Mise en forme des données
  countdata =  as.matrix(countdata)
  library(DESeq2)
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ RNAi + Timing + Batch) ### Ajout de la variable batch
  
  
  # Analyse DESeq2
  deseq = DESeq(deseq)
  
  # Graphique du paramètre de dispersion
  png(paste0(path,i,"_dipression_DESeq2.png"))
    plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()
  
  
  # Récupération des données de comptage normalisées
  tab=counts(deseq,normalized=T)
  
  ####
  

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
  print(paste(i, "-----> début ACP"))
  PCA_plot_generator(data_tab,colors = color,
                            save_path = paste0(path,"ACP/DESeq2_",i,"_"),
                            main = paste0("ACP ",i," (",type,")"))
  
  # Analyse de discrimination linéaire (LDA)
  print(paste(i, "-----> début LDA"))


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