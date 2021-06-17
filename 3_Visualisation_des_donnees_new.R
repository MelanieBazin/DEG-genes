
#### A lancer depuis 4_Analyse.R #####
# Permet de générler des visulaisation 
#  - des donnée avant et après clustering (heatmp et plot)
#  - ACP, LDA, heatmap, clustering herarchique

# i correspond aux names(rnai_list) et donc aux groupe de RNA-seq analysés ensembles


##### Boxplot des comptages normalisés divisé par la taille des gènes #####
print(paste(i, "-----> Creation BoxPlot normalise"))
# pdf(paste0(path,i,"_DESeq_Boxplot.pdf"))
png(paste0(path,i,"_DESeq_Boxplot.png"))
CountBoxplot(data_tab, "DESeq2_seize", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
dev.off()

if (is.element(paste0(i,"_DESeq2-seize.tab"),list.files(paste0("./DATA/DESeq2-seize/")))){
  data_tab_seize = read.table(paste0("./DATA/DESeq2-seize/",i,"_DESeq2-seize.tab"), header = T, sep = "\t")
}else {
  print("Calcule de la table DESeq/taille")
  data_tab_seize = DivideByGeneSeize(data_tab)
  write.table(data_tab_seize,paste0("./DATA/DESeq2-seize/",i,"_DESeq2-seize.tab"), sep="\t",row.names=F,quote=F)
}

# pdf(paste0(path,i,"_DESeq-seize_Boxplot.pdf"))
png(paste0(path,i,"_DESeq-seize_Boxplot.png"))
CountBoxplot(data_tab_seize, "DESeq2_seize", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
dev.off()

print("Boxplot fini")


##### Heatmap et profils  ####
print(paste(i, "-----> Conception des heatmap"))

# Avant moyenne par cluster
ProfilsPNG(save_path = paste0(path,"profils/"), data_tab, condition = i)
ProfilsPDF(save_path = paste0(path,"profils/"), data_tab, condition = i)


MyHeatmaps(path = paste0(path,"Heatmap/"),data_tab, condition = i, sortie = "png")
MyHeatmaps(paste0(path,"HeatmapNoLog/"),data_tab, condition = i, Log = F, sortie = "png")

# Avec calcul des moyennes sur les clusters
mean_data_tab = MeanTabCalculation(data_tab, rnai_list, cluster,i)

ProfilsPNG(save_path = paste0(path,"profils/"), mean_data_tab, moyenne = T, condition = i)
ProfilsPDF(save_path = paste0(path,"profils/"), mean_data_tab, moyenne = T, condition = i)

MyHeatmaps(paste0(path,"Heatmap/"),mean_data_tab, moyenne = T, condition = i, sortie = "png")
MyHeatmaps(paste0(path,"HeatmapNoLog/"),mean_data_tab, moyenne = T, condition = i, Log = F, sortie = "png")



##### Analyse multi-variée des données pour clustering  #####
# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[i]]){
  color[grep(j, color)]=cluster_color[[j]]
}

# Analyse en composante principale
print(paste(i, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = path,
                   main = paste0("ACP ",i," (DESeq2)"),
                   sortie = "png")

# Analyse de discrimination linéaire (LDA)
print(paste(i, "-----> Analyse LDA"))

lda_data_tab=scale(t(data_tab)) #s'assurer que la sd est de 1 et la moyenne à 0 (prédicat des lda)
# summary(apply(lda_data_tab,2,mean)) #verification que la moyenne est à 0 ou très proche
# summary(apply(lda_data_tab,2,sd)) #verification que la sd est à 1

keep = c()
for (l in 1:ncol(lda_data_tab)){
  keep = c(keep, !is.element(T, is.na(lda_data_tab[,l])))
}
lda_data_tab = lda_data_tab[,keep]

lda_model = lda(lda_data_tab, grouping = infodata$Cluster)
# lda_model$prior
# summary(lda_model$scaling)
lda_pred = predict(lda_model)

LDA_plot_generator("LDA",lda_data_tab,infodata, lda_model, path, i, color, sortie = "png")

# print("Teste des prédictions")
# EvaluPrediction("LDA", data_tab, infodata, i, path)  # Evaluer la prédiction


##### Clustering hierarchique  #####
print(paste(i, "-----> Clustering en cours"))
dir.create(paste0(path,"4Cluster/"),recursive=T,showWarnings=F)
for (distance in c("Pearson", "Spearman")){
  png(paste0(path,"4Cluster/",i,"_Matrice_",distance,".png"))
  # pdf(paste0(path,"4Cluster/",i,"_Matrice_",distance,".pdf"))
  # Choisir le mode de calcule des distances
  if (distance == "Pearson"){
    matDist = as.matrix(cor(data_tab))
    p= pheatmap(matDist, main = paste("Pheatmap Pearson DESeq2", i))
    print(p)
    matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
    
  }else if (distance == "Spearman"){
    matDist = as.matrix(cor(data_tab,method="spearman"))
    p= pheatmap(matDist, main = paste("Pheatmap Spearman DESeq2", i))
    print(p)
    matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
  }
  dev.off()
  
  
  for (method in c("kmeans", "HCL")){
    print(paste(distance, method))
    png(paste0(path,"4Cluster/",i,"_Cluster_",method,"_",distance,".png"))
    # pdf(paste0(path,"4Cluster/",i,"_Cluster_",method,"_",distance,".pdf"))
    Clustering(matDist = matDist,
               nb_cluster = 5,
               method = method,
               titre = paste("DESeq2",i),
               colors = color)
    dev.off()
  }
}


