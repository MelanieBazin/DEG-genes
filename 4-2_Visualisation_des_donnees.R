
#### A lancer depuis 4_Analyse.R #####
# Permet de générer des visualisations 
#  - des données avant et après clustering (heatmp et plot)
#  - ACP, heatmap, clustering herarchique

#  condition correspond aux names(rnai_list) et donc aux groupe de RNA-seq analysés ensembles


##### Boxplot des comptages normalisés divisé par la taille des gènes #####
print(paste( condition, "-----> Creation BoxPlot normalise"))
png(paste0(path, "Comptage_bolxplot_DESeq.png"))
CountBoxplot(data_tab, "DESeq2_seize", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
dev.off()

dir.create("./DATA/DESeq2-seize/", recursive=T,showWarnings=F)
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
ExpressionProfils(type = "vst",
                  condition = condition,
                  file = path_dir,
                  select_ID = select_ID)

for (color_type in c("methods","replicates")){
  # Choix de la couleur utilisé
  print(paste( condition, "-----> Setting",color_type,"colors"))
  if (color_type == "methods"){
    # Créaction du vecteur de couleur par méthode de séquencage
    color = Batch_color(data_tab, cluster_list = cluster)
  }else if (color_type == "replicates"){
    # Créaction du vecteur de couleur par groupe de pseudo_réplicat
    color = Culster_color_info(data_tab,infodata)
  }
  
  #### Analyse multi-variée des données pour clustering  ####
  print(paste( condition, "-----> PCA analysis :", color_type))
  PCA_plot_generator(data_tab,
                     colors = color,
                     police_seize = 2,
                     save_path = paste0(path,"PCA_",color_type,"/"),
                     main = paste0("PCA ", condition," (DESeq2)"),
                     sortie = "png")
  
  
  #### Matrice de distance et clusterng hiérarchique  ####
  print(paste( condition, "-----> Hierarchical clustering :", color_type))
  
  matDist = as.matrix(cor(data_tab))
  png(paste0(path, condition,"_Matrice_pearson.png"),  width = 600, height = 600)
  pheatmap(matDist, main = paste("Pheatmap Pearson",  condition), cluster_rows = F, cluster_cols = F)
  dev.off()
  
  matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
  res = hclust(matDist)
  res = as.dendrogram(res)
  labels_colors(res)= as.character(color)[order.dendrogram(res)]
  
  png(paste0(path, condition,"_HCL_",color_type,".png"),  width = 800, height = 200)
  plot(res, main = "pearson_vst")
  dev.off()
}

##### Heatmap ####
print(paste( condition, "-----> Conception des heatmap"))

data_tab_ord = OrderColumn(data_tab, infodata_collapse)

# Avant moyenne par cluster
MyHeatmaps(path = paste0(path,"Visualisation/"),data_tab_ord,infodata_collapse, condition)

print ("Heatmap 1 done")
# Avec calcul des moyennes sur les clusters
mean_data_tab = MeanTabCalculation(data_tab, infodata_collapse)
write.table(mean_data_tab,paste0(path,condition ,"_MEANexpression_table_vst.tab"), sep="\t",row.names=T,quote=F)

MyHeatmaps(paste0(path,"Visualisation/"),mean_data_tab,infodata_collapse, condition,  moyenne = T)

print ("Heatmap mean done")
