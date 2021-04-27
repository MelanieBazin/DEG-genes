options(stringsAsFactors = FALSE)

#############################################
# PCA
#############################################
#Olivier


library(FactoMineR)
library(gtools)

require(ggplot2)

PCA_plot_generator <- function(Expression_Mat, colors,max_dim=3,barplot_max_dim=3,image_prefix="PCA_",show_barplot=T, vline=0, ...) {
  resExp = PCA(t(Expression_Mat), graph = F)
  if(show_barplot) {
    eigenvalues <- resExp$eig
    pdf(paste0(image_prefix,"_PCA_Variance.pdf"))
    barplot(eigenvalues[1:barplot_max_dim, 2], names.arg=1:barplot_max_dim, 
            main = "Variances",
            xlab = "Principal Components",
            ylab = "Percentage of variances",
            col ="gray", ...)
    if(vline!=0) {
      abline(v=vline,lty=2,lwd=2)
    }
    dev.off()
  }    
  
  for (i in 1:dim(combn(1:max_dim,2))[2]) {
    
    gp<-plot.PCA(resExp, axes = combn(1:max_dim,2)[,i], habillage = "ind", col.hab = colors, ...)
    ggsave(paste0(image_prefix,i,".pdf"), plot = gp)
  }
  
  
}

##############
# Dessiner des profils
##############
# Gaëlle

plotGenes <- function(expData, title = "", yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(0, ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}


###################
# Clustering
####################

# Fonction 1 : Lecture des données
F1_lecture_donnée <- function(Nom_de_fichier, Chemin_acces = "./"){
  expMatrix = read.table(paste0(Chemin_acces,Nom_de_fichier), header = T, row.names = 1)
  
  # Si  des gènes ne sont pas exrimer
  if (is.element(T, rowSums(expMatrix) <= ncol(expMatrix))){
    # Mettre à part les genes qui ne sont pas exprimés
    expMatrix_1 = expMatrix[rowSums(expMatrix) <= ncol(expMatrix),]
    
    # Supprimer les genes qui ne sont pas exprimés
    expMatrix = expMatrix[rowSums(expMatrix) > ncol(expMatrix),]
    
    return(list(expMatrix, expMatrix_1)) 
    # L'element retourné par la fonction est une liste contantant les 2 tableau de donnée (celui des gènes non eximer et celui des autres gènes)
  } else {
    return(expMatrix)
    # L'élement retourné est le jeu de donné entier
  }
  
}


# Fonction F2 : Calcul de la matrice de distance
F2_matrice_distance <- function(data, distance){
  # Choisir le mode de calcule des distances
  if (distance == "Euclidean"){
    matDist = dist(data)
  }else if (distance == "Correlation"){
    matDist = as.dist(1 - cor(t(data)))
  } 
}


# Fonction 3 : Application de l’algorithme de regroupement
F3_Algorithme_regroupement <- function(matDist, nb_cluster, method){
  # Choisir le type d'algorithme utilisé pour faire les clusters
  if (method  == "kmeans"){
    res = kmeans(matDist, nb_cluster)
    vecCluster = res$cluster
  }else if(method  == "HCL"){
    res = hclust(matDist)
    vecCluster = cutree(res, nb_cluster)
  }
  
  
  
}


# Fonction 4 : Extraction des profils de gènes pour un cluster donné
F4_Extraction_profil_un_cluster <- function(data, vecCluster, numero_de_cluster){
  geneCluster = names(which(vecCluster == numero_de_cluster))
  cluster = data[geneCluster,]
}



# Fonction 5 : Représentation graphique des résultats
F5_Representation_graphique <- function(data, cluster, graph_type, selected_cluster, distance, method, nb_cluster){
  
  if (is.element(TRUE,graph_type == "profils")) {
    plotGenes(cluster, 
              title = paste("Cluster", selected_cluster, "\n",
                            "Distance :",distance,"-",
                            "Algorithme :", method, "avec ", nb_cluster, "clusters"),
              #yMax = max(data)
              yMax = max(cluster)
    ) 
  }
  
  if (is.element(TRUE,graph_type == "heatmap")){
    
    #Permet d'éviter les erreur génréer par les cluster ne contenat qu'un gènes pour lesquel on ne peux pas produire de heatmap
    if (nrow(cluster)>1){ 
      
      heatmap(as.matrix(cluster),
              Colv = NA, Rowv = NA,
              main  = paste("\n","\n","Cluster", selected_cluster,"\n","Distance :",distance,"-",
                            "Algorithme :", method,  nb_cluster, "clusters"))
    }
    
  }
}


# Fonction finale : fonction permettant de lancer les fonctions précédentes dans l'ordre et qui vas créer les graph pour tous les clusters
Clustering <- function(Nom_de_fichier, Chemin_acces = "./",
                            distance, nb_cluster, method,
                            graph_type){
  
  # Permet de gérer le type de sortie de l'ouverture des fichiers selon qu'1 ou 2 tableau soient générer par la fonction F1
  if(is.list(F1_lecture_donnée(Nom_de_fichier))){
    expMatrix = F1_lecture_donnée(Nom_de_fichier)[[1]]
    expMatrix_1 = F1_lecture_donnée(Nom_de_fichier)[[2]]
  } else {
    expMatrix = F1_lecture_donnée(Nom_de_fichier)
    expMatrix_1 = NULL
  }
  
  
  matDist = F2_matrice_distance(expMatrix, distance)
  
  vecCluster = F3_Algorithme_regroupement(matDist,nb_cluster, method)
  
  if (is.element(FALSE,graph_type == "profils")){
    par(mfrow = c(1,1))
  }else {
    par(mfrow = c(2,2))
  }
  
  # Si un tableau avec les données des gnène non exrpimé est générer alors on fait les graphiques correspondnats
  if(!is.null(expMatrix_1)){
    vecCluster_0 = rep(0, nrow(expMatrix_1))
    names(vecCluster_0) = rownames(expMatrix_1)
    
    cluster = F4_Extraction_profil_un_cluster(expMatrix_1, vecCluster_0,
                                              0)
    F5_Representation_graphique(expMatrix_1, cluster, graph_type,
                                0,
                                distance, method, nb_cluster)
    
  }
  
  # Faire les graphiques pour chacun des clusters générés la la fonction F3
  for (selected_cluster in 1:nb_cluster){
    cluster = F4_Extraction_profil_un_cluster(expMatrix, vecCluster,
                                              selected_cluster)
    
    F5_Representation_graphique(expMatrix, cluster, graph_type,
                                selected_cluster,
                                distance, method, nb_cluster)
  }
  
  
  
  
}

##########










