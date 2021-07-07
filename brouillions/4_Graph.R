options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)

set.seed(10111)

# Récupérer les paramètre de clustring
source("0_Cluster.R")

# Récupérer les fonction necessaire au représentaion graphique et la mise en forme des données
source("3_Visualisation_des_donnees_fonction.R")


analyseName = paste0("Data_DESeq2_toutBatch")

path_dir = paste0("./Analyse/",analyseName,"/")
dir.create(path_dir,recursive=T,showWarnings=F)

condition = names(rnai_list)[1]
# for (condition in names(rnai_list)){
  print(paste("On analyse le jeu de donnee :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))
  
  path = paste0(path_dir,condition ,"/")
  dir.create(path,recursive=T,showWarnings=F)
  
  ##### Analyse DESeq2 ####
  # Ouverture des fichiers countdata sans correction de l'effet Batch
  countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)
  
  # Boxplot des comptages avant normalisation #
  print(paste(condition , "-----> Creation BoxPlot non-normalise"))
  # pdf(paste0(path,condition ,"_row.pdf"))
  png(paste0(path,condition ,"_row.png"))
  CountBoxplot(countdata, "row", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
  dev.off()
  
  # Création du tableau avec les info des colonnes
  infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)
  
  # Ouverture des fichiers countdata avec correction de l'effet Batch pour ICL7
  countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_pour_DESeq.tab"), sep="\t",row.names=1,header =  T)
  
  # Ouverture des fichiers countdata SANS correction de l'effet Batch pour ICL7
  # countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)
  # path = paste0(path, "sans_correction_Batch/")
  # dir.create(path,recursive=T,showWarnings=F)

  # Creataion de l'objet DESeq2
  countdata = as.matrix(countdata)
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Condition
                                 # design   = ~ Feeding + Cluster
  )
  
  
  # Analyse DESeq2
  print(paste(condition , "-----> Analyse DESeq2"))
  deseq = DESeq(deseq)
  
  # Graphique du paramètre de dispersion
  # pdf(paste0(path,condition ,"_dipression_DESeq2.pdf"))
  png(paste0(path,condition ,"_dipression_DESeq2.png"))
  plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()
  
  
  # Récupération des données de comptage normalisées
  data_tab = counts(deseq,normalized=T)
  
  write.table(data_tab,paste0(path,condition ,"_expression_table_normaliserDESeq2.tab"), sep="\t",row.names=F,quote=F)
  
  #### Lancer les visulalisation des données ####
  source("3_Visualisation_des_donnees_new.R")
# }