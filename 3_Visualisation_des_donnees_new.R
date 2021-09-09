
#### A lancer depuis 4_Analyse.R #####
# Permet de générler des visulaisation 
#  - des donnée avant et après clustering (heatmp et plot)
#  - ACP, LDA, heatmap, clustering herarchique

#  condition correspond aux names(rnai_list) et donc aux groupe de RNA-seq analysés ensembles


##### Boxplot des comptages normalisés divisé par la taille des gènes #####
print(paste( condition, "-----> Creation BoxPlot normalise"))
png(paste0(path, "Comptage_bolxplot_DESeq.png"))
CountBoxplot(data_tab, "DESeq2_seize", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
dev.off()

if (is.element(paste0( condition,"_DESeq2-seize.tab"),list.files(paste0("./DATA/DESeq2-seize/")))){
  data_tab_seize = read.table(paste0("./DATA/DESeq2-seize/", condition,"_DESeq2-seize.tab"), header = T, sep = "\t")
}else {
  print("Calcule de la table DESeq/taille")
  data_tab_seize = DivideByGeneSeize(data_tab)
  write.table(data_tab_seize,paste0("./DATA/DESeq2-seize/", condition,"_DESeq2-seize.tab"), sep="\t",row.names=F,quote=F)
}

png(paste0(path, "Comptage_bolxplot_row_DESeq-seize.png"))
CountBoxplot(data_tab_seize, "DESeq2_seize", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
dev.off()

print("Boxplot fini")

##### Dessiner les profils d'une selection de gènes ####

ExpressionProfils(type = "DESeq2", 
                  condition = "tout", 
                  file = path_dir,
                  select_ID = select_ID)

##### Analyse multi-variée des données pour clustering  #####
# Créaction du vecteur de couleur par anné de séquancage
color = colnames(data_tab)
for (j in rnai_list[[ condition]]){
  if(sum(grepl(j, colnames(data_tab)))>0){
    color[grep(j, color)]=seq_color[[j]]
  }
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"Visualisation/ACP/color2/"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  if(sum(grepl(j, colnames(data_tab)))>0){
    print(color[grep(j, color)])
    print(cluster_color[[j]])
    color[grep(j, color)]=cluster_color[[j]]

  }
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"Visualisation/ACP/color4/"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")


##### Clustering hierarchique  #####
print(paste( condition, "-----> Clustering en cours"))
dir.create(paste0(path,"4Cluster/"),recursive=T,showWarnings=F)
for (distance in c("Pearson", "Spearman")){
  png(paste0(path,"4Cluster/", condition,"_Matrice_",distance,".png"),  width = 600, height = 600)
  # Choisir le mode de calcule des distances
  if (distance == "Pearson"){
    matDist = as.matrix(cor(data_tab))
    p= pheatmap(matDist, main = paste("Pheatmap Pearson DESeq2",  condition))
    print(p)
    matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
    
  }else if (distance == "Spearman"){
    matDist = as.matrix(cor(data_tab,method="spearman"))
    p= pheatmap(matDist, main = paste("Pheatmap Spearman DESeq2",  condition))
    print(p)
    matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
  }
  dev.off()
  
  
  for (method in c("kmeans", "HCL")){
    print(paste(distance, method))
    png(paste0(path,"4Cluster/", condition,"_Cluster_",method,"_",distance,".png"),  width = 800, height = 600)
    Clustering(matDist = matDist,
               nb_cluster = 5,
               method = method,
               titre = paste("DESeq2", condition),
               colors = color)
    dev.off()
  }
}

##### Heatmap et profils  ####
print(paste( condition, "-----> Conception des heatmap"))

# Avant moyenne par cluster
ProfilsPNG(save_path = paste0(path,"Visualisation/profils/"), data_tab, condition =  condition)
ProfilsPDF(save_path = paste0(path,"Visualisation/profils/"), data_tab, condition =  condition)

# problème avec EZL1 a cause des EZL1bis
# MyHeatmaps(path = paste0(path,"Visualisation/Heatmap/"),data_tab, condition =  condition, sortie = "png")
# MyHeatmaps(paste0(path,"Visualisation/HeatmapNoLog/"),data_tab, condition =  condition, Log = F, sortie = "png")

# Avec calcul des moyennes sur les clusters
mean_data_tab = MeanTabCalculation(data_tab, rnai_list, cluster, condition)

ProfilsPNG(save_path = paste0(path,"Visualisation/profils/"), mean_data_tab, moyenne = T, condition =  condition)
ProfilsPDF(save_path = paste0(path,"Visualisation/profils/"), mean_data_tab, moyenne = T, condition =  condition)

# MyHeatmaps(paste0(path,"Visualisation/Heatmap/"),mean_data_tab, moyenne = T, condition =  condition, sortie = "png")
# MyHeatmaps(paste0(path,"Visualisation/HeatmapNoLog/"),mean_data_tab, moyenne = T, condition =  condition, Log = F, sortie = "png")


