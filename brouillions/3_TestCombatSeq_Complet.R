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


analyseName = "2022-02-21_Test_Combatseq"

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
  # # Ouerture du tableau de donnée
  # countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_",correction,".tab"), sep="\t",row.names=1,header =  T)
  # 
  # # Création du tableau avec les info des colonnes
  # infodata = CreatInfoData(countdata, conditions = condition , rnai_list, cluster)
  # 
  # if (condition == names(rnai_list["HiSeqvsNextSeq"])){
  #   # Se limiter au données présentes dans les 2 séquancages
  #   time = unique(infodata$Timing[which(infodata$Seq_method == "NextSeq")])
  #   countdata = countdata[,is.element(infodata$Timing, time)]
  #   infodata = infodata[is.element(infodata$Timing, time),]
  # }
  # 
  # 
  # # Creataion de l'objet DESeq2
  # countdata = as.matrix(countdata)
  # deseq = DESeqDataSetFromMatrix(countData = countdata,
  #                                colData  = infodata,
  #                                design   = ~ Samples)
  # 
  # 
  # # Analyse DESeq2
  # print(paste(condition , "-----> DESeq2 analysis"))
  # deseq = DESeq(deseq)
  # 
  # # Récupération des données de comptage normalisées
  # data_tab = assay(vst(deseq, blind = T))
  # 
  # write.table(data_tab,paste0(path,condition,"_expression_table_",correction,".tab"), sep="\t",row.names=T,quote=F)
  # write.table(infodata,paste0(path,condition,"_infodata_",correction,".tab"), sep="\t",row.names=T,quote=F)
  
  #### Ouverture des fichier ####
  data_tab = read.table(paste0(path,condition,"_expression_table_",correction,".tab"), sep="\t", header = T)
  infodata = read.table(paste0(path,condition,"_infodata_",correction,".tab"), sep="\t", header = T)
  
  #### Créer les vecteurs de couleurs  ####
  # Créaction du vecteur de couleur 
  cluster2 = list(
    ICL7 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
    ICL7bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
    EZL1 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
    EZL1bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2))
  )
  
  # Créaction du vecteur de couleur par méthode de séquencage
  color = Batch_color(data_tab, cluster_list = cluster2)
  
  
  # names(color) = colnames(data_tab)
  
  #### Analyse multi-variée des données pour clustering  ####
  print(paste( condition, "-----> PCA analysis : methods"))
  PCA_plot_generator(data_tab,
                     colors = color,
                     police_seize = 5,
                     save_path = paste0(path,"PCA_methods/"),
                     main = paste0("PCA ", condition," (DESeq2)"),
                     sortie = "pdf",
                     rename = F)
  
  PCA_ggplot_generator(data_tab,
                       infodata,
                       police_seize = 5,
                       point_seize = 3,
                       save_path = paste0(path,"PCA_methods/ggplot/"),
                       main = paste0("PCA ", condition," (DESeq2)"),
                       sortie = "pdf",
                       rename = F)
  
  
}


}
