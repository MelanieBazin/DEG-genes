options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)

set.seed(10111)

# Récupérer les paramètre de clustring
source("0_Cluster.R")

# Récupérer les fonction necessaire au représentaion graphique et la mise en forme des données
source("0_Visualisation_fonction.R")


analyseName = paste0("Test_Combatseq")

analyseName = paste0(Sys.Date(),"_", analyseName)
path_dir = paste0("./Analyse/",analyseName,"/")
dir.create(path_dir,recursive=T,showWarnings=F)
dir.create("./DATA/DESeq2/",recursive=T,showWarnings=F)

condition = names(rnai_list["HiSeqvsNextSeq"])

print(paste("On analyse le jeu de donnee :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))

path = paste0(path_dir,condition ,"/")
dir.create(path,recursive=T,showWarnings=F)

##### Analyse DESeq2 ####
for (correction in c("corrected","uncorrected")){
  if (correction == "corrected"){
    c = "AVEC"
  }else if (correction == "uncorrected"){
    c = "SANS"
  }
  dir.create(paste0(path,"/",c,"_ComBat/"),recursive=T,showWarnings=F)
  print(paste("Teste de conditions",c,"correction ComBatSeq"))
  
  countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_",correction,".tab"), sep="\t",row.names=1,header =  T)
  
  # Création du tableau avec les info des colonnes
  infodata = CreatInfoData(countdata, conditions = condition , rnai_list, cluster)
  
  if (condition == names(rnai_list["HiSeqvsNextSeq"])){
    # Se limiter au données présentes dans les 2 séquancages
    time = unique(infodata$Timing[which(infodata$Seq_method == "NextSeq")])
    countdata = countdata[,is.element(infodata$Timing, time)]
    infodata = infodata[is.element(infodata$Timing, time),]
  }
  
  
  # Creataion de l'objet DESeq2
  countdata = as.matrix(countdata)
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Samples)
  
  
  # Analyse DESeq2
  print(paste(condition , "-----> Analyse DESeq2"))
  deseq = DESeq(deseq)
  
  # Récupération des données de comptage normalisées
  data_tab = assay(vst(deseq, blind = T))
  write.table(data_tab,paste0("./DATA/DESeq2/",condition,"_expression_table_",correction,".tab"), sep="\t",row.names=T,quote=F)
  
  ##### Analyse multi-variée des données pour clustering  #####
  print(paste( condition, "-----> PCA Analyse"))
  # Créaction du vecteur de couleur par année de sequancage
  seq_color = Batch_color(condition, collapse = T)
  names(seq_color) = colnames(data_tab)
  
  print(seq_color)
  
  # Analyse en composante principale
  PCA_plot_generator(data_tab,
                     colors = seq_color,
                     save_path = paste0(path,"/",c,"_ComBat/PCA/"),
                     main = paste0("PCA ", condition," (DESeq2)"),
                     sortie = "png")
  
  # Créaction du vecteur de couleur par cluster
  if (condition == names(rnai_list["HiSeqvsNextSeq"])){
    cluster2 = list(
      ICL7 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
      ICL7bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
      EZL1 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
      EZL1bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2))
    )
  }
  
  clust_color = Culster_color(condition, cluster_list = cluster2)
  names(clust_color) = colnames(data_tab)
  
  print(clust_color)
  
  # Analyse en composante principale
  PCA_plot_generator(data_tab,
                     colors = clust_color,
                     save_path = paste0(path,"/",c,"_ComBat_cluster/PCA/"),
                     main = paste0("PCA ", condition," (DESeq2)"),
                     sortie = "png")
  
  print(paste( condition, "-----> Clustering hiérachique en cours"))
  dir.create(paste0(path,c,"_ComBat/Cluster/"),recursive=T,showWarnings=F)
  for (distance in c("Pearson", "Spearman")){
    
    # Choisir le mode de calcule des distances
    if (distance == "Pearson"){
      matDist = as.matrix(cor(data_tab))
      p= pheatmap(matDist, main = paste("Pheatmap Pearson DESeq2",  condition), cluster_rows = F, cluster_cols = F)
      matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
      
    }else if (distance == "Spearman"){
      matDist = as.matrix(cor(data_tab,method="spearman"))
      p= pheatmap(matDist, main = paste("Pheatmap Spearman DESeq2",  condition), cluster_rows = F, cluster_cols = F)
      matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
    }
    
    png(paste0(path,c,"_ComBat/Cluster/", condition,"_Matrice_",distance,".png"),  width = 600, height = 600)
    print(p)
    dev.off()
    
    
    for (method in c("kmeans", "HCL")){
      print(paste(distance, method))
      png(paste0(path,c,"_ComBat/Cluster/", condition,"_Cluster_",method,"_",distance,"_clutercolor.png"),  width = 800, height = 600)
      Clustering(matDist = matDist,
                 nb_cluster = 4,
                 method = method,
                 titre = paste("DESeq2", condition),
                 colors = clust_color)
      dev.off()
      
      png(paste0(path,c,"_ComBat/Cluster/", condition,"_Cluster_",method,"_",distance,"_methodcolor.png"),  width = 800, height = 600)
      Clustering(matDist = matDist,
                 nb_cluster = 4,
                 method = method,
                 titre = paste("DESeq2", condition),
                 colors = seq_color)
      dev.off()
    }
  }
}


if (condition == names(rnai_list["HiSeqvsNextSeq"])){
  print("Correlation matrix")
  #### Martice de corrélation ####
  
  sansBatch = read.table(paste0("./DATA/DESeq2/",condition,"_expression_table_uncorrected.tab"), sep="\t",row.names=1,header =  T)
  
  avecBatch = read.table(paste0("./DATA/DESeq2/",condition,"_expression_table_corrected",".tab"), sep="\t",row.names=1,header =  T)
  colnames(avecBatch) = paste0(colnames(avecBatch),"_corrected")
  
  data_tab = cbind(sansBatch, avecBatch)
  data_tab = data_tab[,c(grep("EZL1", colnames(data_tab)),grep("ICL7", colnames(data_tab)))]
  
  dir.create(paste0(path,"AVEC_ComBat"),recursive=T,showWarnings=F)
  for (distance in c("Pearson", "Spearman")){
    
    # Choisir le mode de calcule des distances
    if (distance == "Pearson"){
      matDist = as.matrix(cor(data_tab))
      p= pheatmap(matDist, main = paste("Pheatmap Pearson DESeq2",  condition), cluster_rows = F, cluster_cols = F)
      matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
      
    }else if (distance == "Spearman"){
      matDist = as.matrix(cor(data_tab,method="spearman"))
      p= pheatmap(matDist, main = paste("Pheatmap Spearman DESeq2",  condition), cluster_rows = F, cluster_cols = F)
      matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
    }
    
    png(paste0(path,"AVEC_ComBat/", condition,"_Matrice_",distance,".png"),  width = 600, height = 600)
    print(p)
    dev.off()
  }
  
  
} 


