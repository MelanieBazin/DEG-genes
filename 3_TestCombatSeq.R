options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)
library(dendextend)

set.seed(10111)

# Récupérer les paramètre de clustring
source("0_Cluster.R")

# Récupérer les fonction necessaire au représentaion graphique et la mise en forme des données
source("0_Visualisation_fonction.R")


analyseName = paste0("Test_Combatseq")

analyseName = paste0(Sys.Date(),"_", analyseName)
path_dir = paste0("./Analyse/",analyseName,"/")
dir.create("./DATA/DESeq2/",recursive=T,showWarnings=F)

condition = names(rnai_list["HiSeqvsNextSeq"])

print(paste("Analysis of dataset :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))

for (correction in c("corrected","uncorrected")){
  print(paste("Test of conditions",correction, "by ComBatSeq"))
  
  path = paste0(path_dir,condition ,"/",correction,"/")
  dir.create(path,recursive=T,showWarnings=F)

  #### Analyse DESeq2 ####
  # Ouerture du tableau de donnée
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
  print(paste(condition , "-----> DESeq2 analysis"))
  deseq = DESeq(deseq)
  
  # Récupération des données de comptage normalisées
  data_tab = assay(vst(deseq, blind = T))
  
  write.table(data_tab,paste0(path,condition,"_expression_table_",correction,".tab"), sep="\t",row.names=T,quote=F)
  write.table(infodata,paste0(path,condition,"_infodata_",correction,".tab"), sep="\t",row.names=T,quote=F)
  
  #### Créer les vecteurs de couleurs  ####
  # Créaction du vecteur de couleur 
  if (condition == names(rnai_list["HiSeqvsNextSeq"])){
    cluster2 = list(
      ICL7 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
      ICL7bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
      EZL1 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
      EZL1bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2))
    )
  }else{
    cluster2 = cluster
  }
  
  for (color_type in c("methods","replicates")){
    # Choix de la couleur utilisé
    print(paste( condition, "-----> Setting",color_type,"colors"))
    if (color_type == "methods"){
      # Créaction du vecteur de couleur par méthode de séquencage
      color = Batch_color(data_tab, cluster_list = cluster2)
    }else if (color_type == "replicates"){
      # Créaction du vecteur de couleur par groupe de pseudo_réplicat
      color = Culster_color(data_tab, cluster_list = cluster2)
    }
    
    names(color) = colnames(data_tab)
    
    #### Analyse multi-variée des données pour clustering  ####
    print(paste( condition, "-----> PCA analysis :", color_type))
    PCA_plot_generator(data_tab,
                       colors = color,
                       save_path = paste0(path,"PCA_",color_type,"/"),
                       main = paste0("PCA ", condition," (DESeq2)"),
                       sortie = "png")
    
    
    #### Matrice de distance et clusterng hiérarchique  ####
    print(paste( condition, "-----> Hierarchical clustering :", color_type))
    
    matDist = as.matrix(cor(data_tab))
    p= pheatmap(matDist, main = paste("Pheatmap Pearson",  condition), cluster_rows = F, cluster_cols = F)
    
    png(paste0(path, condition,"_Matrice_pearson.png"),  width = 600, height = 600)
    print(p)
    dev.off()
    
    matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
    res = hclust(matDist)
    res = as.dendrogram(res)
    labels_colors(res)= as.character(color)[order.dendrogram(res)]
    
    png(paste0(path, condition,"_HCL_",color_type,".png"),  width = 800, height = 600)
    plot(res, main = "pearson_vst")
    dev.off()
  }
  
  
}


if (condition == names(rnai_list["HiSeqvsNextSeq"])){
  print("Correlation matrix")
  #### Martice de corrélation ####
  
  sansBatch = read.table(paste0(path_dir,condition,"/uncorrected/",condition,"_expression_table_uncorrected.tab"), sep="\t",row.names=1,header =  T)
  
  avecBatch = read.table(paste0(path_dir,condition,"/corrected/",condition,"_expression_table_corrected",".tab"), sep="\t",row.names=1,header =  T)
  colnames(avecBatch) = paste0(colnames(avecBatch),"_corrected")
  
  data_tab = cbind(sansBatch, avecBatch)
  data_tab = data_tab[,c(grep("EZL1", colnames(data_tab)),grep("ICL7", colnames(data_tab)))]
  
  
  # Choisir le mode de calcule des distances
  
  matDist = as.matrix(cor(data_tab))
  p= pheatmap(matDist, main = paste("Pheatmap Pearson DESeq2",  condition), cluster_rows = F, cluster_cols = F)
  matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
  
  png(paste0(path_dir,condition,"/" , condition,"_Matrix.png"),  width = 600, height = 600)
  print(p)
  dev.off()
  
} 
