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

condition = names(rnai_list)[2]

analyseName = paste0("Clustering")
analyseName = paste0(Sys.Date(),"_", analyseName)

path_dir = paste0("./Analyse/",analyseName)
dir.create(path_dir,recursive=T,showWarnings=F)

for (corr in c("Corrected", "Uncorrected")){
  
  path = paste0(path_dir,"/", condition ,"/", corr, "/")
  dir.create(path,recursive=T,showWarnings=F)
  dir.create(paste0(path,"/Visualisation/Cluster/"),recursive=T,showWarnings=F)
  
  if (corr == "Uncorrected"){
    # Ouverture des fichiers sans correction de l'effet Batch
    countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)
  }else if (corr == "Corrected"){
    # Ouverture des fichiers avec correction de l'effet Batch sans la condition de groupe
    countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_pour_DESeq_v1.tab"), sep="\t",row.names=1,header =  T)
  }
  
  ##### Analyse DESeq2 ####
  # Création du tableau avec les info des colonnes
  infodata = CreatInfoData(countdata, conditions = condition , rnai_list, cluster)
  
  # Créataion de l'objet DESeq2
  countdata = as.matrix(countdata)
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Condition)
  
  # Definition des réplicats techniques
  deseq = collapseReplicates(deseq,
                             groupby = deseq$Samples,
                             run     = deseq$Names)
  
  color = c()
  for (r in colnames(deseq)){
    clust = infodata[r, "Cluster"]
    if (clust == "VEG"){
      color = c(color, clust_color["veg"])
    } else if(clust == "EARLY" ){
      color = c(color, clust_color["early"])
    } else if(clust == "INTER" ){
      color = c(color, clust_color["inter"])
    } else if(clust == "LATE" ){
      color = c(color, clust_color["late"])
    }
  }
  
  
  # Analyse DESeq2
  deseq = DESeq(deseq)
  
  ### Recupération des donnée normalisée ####
  # Variance stabilizing transformation
  vsd = assay(vst(deseq, blind = T))
  write.table(vsd,paste0(path,condition ,"_expression_table_vst.tab"), sep="\t",row.names=T,quote=F)
  
  PCA_plot_generator(vsd,
                     colors = color,
                     save_path = paste0(path,"Visualisation/ACP/"),
                     main = paste0("ACP ", condition," (vst)"),
                     sortie = "png")
  
  matDist = as.dist(1-cor(log2(vsd+1), method="pearson"))
  res = hclust(matDist)
  res = as.dendrogram(res)
  labels_colors(res)= as.character(color)[order.dendrogram(res)]
  
  png(paste0(path,"/Visualisation/Cluster/", condition,"_Cluster_pearson_vst.png"),  width = 800, height = 600)
  plot(res, main = "pearson_vst")
  dev.off()
  
  res = kmeans(matDist, 4)
  png(paste0(path,"/Visualisation/Cluster/", condition,"_Kmean_pearson_vst.png"),  width = 800, height = 600)
  fviz_cluster(res, data = matDist, geom = c("point",  "text"), labelsize = 10, repel = T, 
               show.clust.cent = F, ellipse = T, ggtheme = theme_bw(),
               main = "pearson_vst_kmean", 
               xlab = "Principal Component 1",
               ylab = "Principal Component 2")
  dev.off()
  
  # MyHeatmaps(path = paste0(path,"/Visualisation/Cluster/"),
  #            vsd, infodata, condition)
}