options(stringsAsFactors = FALSE)

OpenTabInList <- function(type){
  path = paste0("./DATA/", type)
  count = paste(path,list.files(path), sep = "/")
  
  list_count  = list()
  for (i in count){
    tab = read.table(i, sep = "\t", header = T)
    name = sub(paste0("./DATA/",type,"/"),"",sub(paste0("_expression_table_",type,".tab"),"",i))
    list_count[[name]] = tab
  } 
return(list_count)
}
  
ConcatTab <- function(type, conditions = NULL){
  annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")
  path = paste0("./DATA/", type)
  count = list.files(path)
  if (type == "EXPRESSION"){
    extention = ".tab"
  }else{
    extention = paste0("_expression_table_",type,".tab")
  }
  
  if (!is.null(conditions)){
    count =  count[which(is.element(count, paste0(conditions,extention)))]
  }
  
  count = paste(path,count, sep = "/")
  tab_count = array(annotation$ID, dim = c(nrow(annotation),1))
  colnames(tab_count)="ID"
  for (i in count){
    tab = read.table(i, sep = "\t", header = T)
    name = sub(".tab","",sub(paste0("./DATA/",type,"/"),"",sub(paste0("_expression_table_",type),"",i)))
    colnames(tab)[2:ncol(tab)]=paste(name, colnames(tab)[2:ncol(tab)], sep = "_")
    tab_count = merge(tab_count, tab, by = "ID")
  } 
  return(tab_count)
}

CountBoxplot <- function (tab, type){
  boxplot(log(tab + 1), ylab = "count values (log scale)",
          main = paste0("Count data (",type,")"), xaxt="n", yaxt="n")
  axis(side = 1, labels = FALSE, tick = F)
  axis(side = 2,
       ## Rotate labels perpendicular to y-axis.
       las = 2,
       ## Adjust y-axis label positions.
       mgp = c(3, 0.75, 0))
  text(x = 1:ncol(tab),
       ## Move labels to just below bottom of chart.
       y = par("usr")[3] - 0.1,
       ## Use names from the data list.
       labels = colnames(tab),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 90,
       ## Adjust the labels to almost 100% right-justified.
       adj = 1,
       ## Increase label size.
       cex = 0.5)
}

CreatInfoData2 <- function(conditions=NULL){
  countdata= ConcatTab("EXPRESSION")
  if (colnames(countdata)[1]=="ID"){
    countdata=countdata[,-1]
  }
  infodata=matrix(NA,nrow = ncol(countdata), ncol=3)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Name","RNAi","Timing")
  infodata[,"Name"] = row.names(infodata)

  CTIP = c("T0", "T5.5", "T12.5", "T25", "Veg")
  CTIP_CTRL = c("T0", "T5", "T10", "T20", "T30", "Veg")
  ICL7 = c("T0", "T5", "T10", "T20", "T35", "T50", "Veg")
  KU80c = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
  ND7 = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
  PGM = c( "T2", "T5", "T10", "T20", "T30", "T40", "Veg")
  XRCC4 = c( "T2", "T7", "T22", "T32","Veg")
  XRCC4_CTRL = c( "T2", "T7", "T22", "T32","Veg")
  
  infodata[,"Timing"] = c(CTIP, CTIP_CTRL , ICL7, KU80c , ND7, PGM, XRCC4, XRCC4_CTRL)
  
  condi = sub(".tab","",list.files("./DATA/EXPRESSION"))
  l = list(CTIP, CTIP_CTRL , ICL7, KU80c , ND7, PGM, XRCC4, XRCC4_CTRL)
  rnai = c()
  for (i in 1:length(l)){
    rnai = c(rnai, rep(condi[i],length(l[[i]])))
  }
  infodata[,"RNAi"] = rnai
  infodata = as.data.frame(infodata)
  
  if (!is.null(conditions)){
    infodata =  infodata[which(is.element(infodata$RNAi, conditions)),]
  }
return(infodata)
}

#############################################
# PCA
#############################################
#Olivier


library(FactoMineR)
library("factoextra")
library(ggplot2)
library(gtools)

PCA_plot_generator <- function(Expression_Mat, colors,save_path, main,max_dim=3,barplot_max_dim=3,image_prefix="PCA_",show_barplot=T, vline=0, ...) {
  resExp = PCA(t(Expression_Mat), graph = F)
  
  if(show_barplot) {
    eigenvalues <- resExp$eig
    pdf(paste0(save_path,image_prefix,"_PCA_Variance.pdf"))
      barplot(eigenvalues[1:barplot_max_dim, 2], names.arg=1:barplot_max_dim, 
              main = "Variances",
              xlab = "Principal Components",
              ylab = "Percentage of variances",
              col ="gray")
      if(vline!=0) {
        abline(v=vline,lty=2,lwd=2)
      }
    dev.off()
  }    
  
  for (i in 1:dim(combn(1:max_dim,2))[2]) {
    
    gp<-plot.PCA(resExp, axes = combn(1:max_dim,2)[,i], habillage = "ind", col.hab = colors, title = main,
                 ggoptions = list(size=1.75))
    ggsave(paste0(save_path,image_prefix,i,".png"), device = "png", plot = gp)
    
  }
 return(resExp)
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

#install.packages("pheatmap")
library("pheatmap")
# Fonction F2 : Calcul de la matrice de distance
F2_matrice_distance <- function(data_tab, distance){
  # Choisir le mode de calcule des distances
  if (distance == "Pearson"){
    matDist = as.matrix(cor(data_tab))
    pheatmap(matDist, main = "Pheatmap Pearson")
    matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
    
  }else if (distance == "Spearman"){
    matDist = as.matrix(cor(data_tab,method="spearman"))
    pheatmap(matDist, main = "Pheatmap Spearman")
    matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
  }
}


# Fonction 5 : Représentation graphique des résultats
F5_Representation_graphique <- function(expMatrix, cluster, graph_type, selected_cluster, distance, method, nb_cluster){
  if (is.element(TRUE,graph_type == "profils")) {
    plotGenes(cluster, 
              title = paste("Cluster", selected_cluster, "\n",
                            "Distance :",distance,"-",
                            "Algorithme :", method, "avec ", nb_cluster, "clusters"),
              #yMax = max(expMatrix)
              yMax = max(cluster)
    ) 
  }
  
  if (is.element(TRUE,graph_type == "heatmap")){
    
    #Permet d'éviter les erreur génréer par les cluster ne contenat qu'un gènes pour lesquel on ne peux pas produire de heatmap
    if (nrow(cluster)>1){ 
      
      heatmap(as.matrix(cluster),
              Colv = NA, Rowv = NA, cexCol = 0.2,
              main  = paste("\n","\n","Cluster", selected_cluster,"\n","Distance :",distance,"-",
                            "Algorithme :", method,  nb_cluster, "clusters"))
    }
    
  }
}


# Fonction finale : fonction permettant de lancer les fonctions précédentes dans l'ordre et qui vas créer les graph pour tous les clusters
Clustering <- function(data_tab, Chemin_acces = "./",
                            distance, nb_cluster, method, graph_type){
  

  ## Créaction de la matrice de distance
  matDist = F2_matrice_distance(data_tab, distance)
  
  ## Créaction d'un vecteur contennat le clusering calculé a partir de la matrice de distance
  # Choisir le type d'algorithme utilisé pour faire les clusters
  if (method  == "kmeans"){
    res = kmeans(matDist, nb_cluster)
    vecCluster = res$cluster
    
  }else if(method  == "HCL"){
    res = hclust(matDist)
    vecCluster = cutree(res, nb_cluster)
    
    #Fait un dendrogramme
    p= plot(res, main = paste("Dendrogramme", nb_cluster, "cluster calculé par", method, distance))
    print(p)
  }
  
  ## Sépare la fentêtre de grahique en 4 seulemnt si seul les profils sont demmandés
  if (is.element(FALSE,graph_type == "profils")){
    par(mfrow = c(1,1))
  }else {
    par(mfrow = c(2,2))
  }
  

  ## Faire les graphiques pour chacun des clusters générés la la fonction F3
  for (selected_cluster in 1:nb_cluster){
    geneCluster = names(which(vecCluster == selected_cluster))
    cluster = data_tab[,geneCluster]
    
    F5_Representation_graphique(data_tab, cluster, 
                                graph_type, selected_cluster, 
                                distance, method, nb_cluster)
  }

  
}

##########










