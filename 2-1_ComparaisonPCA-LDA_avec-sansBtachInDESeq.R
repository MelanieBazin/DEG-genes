source("2_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")
source("4_Functions.R")
library(sva)
library(DESeq2)
library(limma)
library(caret)
library(e1071)
library(MASS)
set.seed(10111)

analyseName = paste0("DESeq2_test05_NewCluster-BtachICL7")

path_dir = paste0("./Analyse/",analyseName,"/")
dir.create(path_dir,recursive=T,showWarnings=F)
# Utilisation uniquement des normalisations DESeq2
type = "DESeq2"

# Definiton des variables DESeq2
FC = 1.5 #Mini 1.5 -> XRCC4 = 2
pvalue = 0.05 #Maxi 0.05 -> XRCC4 = 0.01
RNAi = "PGM"

### Création  de liste de gènes filtrés (retirés de l'analyse) ###
Filtering= list()
Filtering= NULL
# par exemple : retirer les gènes qui sont DE pendant une manip de silencing, 
# ou entre plusieurs manip control ==> des faux positifs


### Vecteur de couleur pour les heatmap
hmcol = colorRampPalette(brewer.pal(10,"RdBu"))(255)
#hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
hmcol = rev(hmcol)

### Limiter le fichier annotation aux gènes avec synonyme ###
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
annotation_synonyms = annotation[annotation$SYNONYMS != "",]
rownames(annotation)=annotation$ID


i = "CTIPseulctrl2020"  

for (i in names(rnai_list)){
  
  path = paste0(path_dir,i,"/")
  dir.create(path,recursive=T,showWarnings=F)
  
  ##### Création du tableau de donnée à analyser ensemble ####
  #Ouverture des fichiers et création de l'objet countdata
  countdata = ConcatTab(type = "EXPRESSION", conditions = rnai_list[[i]])

  # Création du tableau avec les info des colonnes
  infodata = CreatInfoData3(countdata, conditions = i, rnai_list, cluster)
  
  # Correction de l'effet batch avec ComBat
  if (length(grep("ICL7", colnames(countdata)))>0){
    print(paste(i, "-----Correction de l'effet Batch"))
    countdata = ComBat_seq(countdata, batch = infodata$Labo, group = infodata$Cluster)
  }

  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Feeding + Cluster)
  
  
  # Analyse DESeq2
  print(paste(i, "-----> Analyse DESeq2"))
  deseq = DESeq(deseq)

  
  # Graphique du paramètre de dispersion
  png(paste0(path,i,"_dipression_DESeq2.png"))
    plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()
  
  
  # Récupération des données de comptage normalisées
  data_tab=counts(deseq,normalized=T)
  
  write.table(data_tab,paste0("./DATA/DESeq2/",i,"_expression_table_DESeq2.tab"), sep="\t",row.names=F,quote=F)
  
  
  # Boxplot des comptages normalisés divisé par la taille des gènes
  print(paste(i, "-----> Création BoxPlot normalisé"))
  png(paste0(path,i,"_DESeq_Boxplot.png"))
    CountBoxplot(data_tab, "DESeq2_seize", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
  dev.off()
  
  data_tab_seize = DivideByGeneSeize(data_tab)
  write.table(data_tab_seize,paste0("./DATA/DESeq2-seize/",i,"_expression_table_DESEQsurseize.tab"), sep="\t",row.names=F,quote=F)
  
  data_tab_seize = read.table(paste0("./DATA/DESeq2-seize/",i,"_expression_table_DESEQsurseize.tab"), header = T, sep = "\t")
  png(paste0(path,i,"_DESeq-seize_Boxplot.png"))
    CountBoxplot(data_tab_seize, "DESeq2_seize", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
  dev.off()
  
  print("Boxplot terminé")
  
  
  # Créaction du vecteur de couleur par cluster
  color = colnames(data_tab)
  for (j in rnai_list[[i]]){
    color[grep(j, color)]=cluster_color[[j]]
  }

  # Analyse en composante principale
  print(paste(i, "-----> Analyse ACP"))
  PCA_plot_generator(data_tab,colors = color,
                      save_path = path,
                      main = paste0("ACP ",i," (",type,")"))
  
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

  LDA_plot_generator("LDA",lda_data_tab,infodata, lda_model, path, i, color)

  # print("Teste des prédictions")
  # EvaluPrediction("LDA", data_tab, infodata, i, path)  # Evaluer la prédiction

  #### SVM
  # print(paste(i, "-----> Analyse SVM"))
  # svm_model = svm(t(data_tab), y = as.factor(infodata$Cluster), scale = F, kernel = "radial", cost = 5)
  # 
  # svm_tab = cbind(t(data_tab), infodata$Cluster)
  # colnames(svm_tab)[ncol(svm_tab)]="Cluster"
  # 
  # plot(svm_model,as.datat.frame(t(data_tab)), y ~. )
  # 
  # 
  # 
  # trctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3)
  # svm_model = train(t(data_tab), y = as.factor(infodata$Cluster), method = "svmRadialWeights",
  #                     trControl=trctrl)
  # svm_pred = predict(svm_model)
  # 
  # plot
  # 
  # LDA_plot_generator("SVM",t(data_tab),infodata, svm_model, path, i, color)
  # 
  # print("Teste des prédictions")
  # EvaluPrediction("SVM", data_tab)  # Evaluer la prédiction
  
  
  
  # Analyse de clusering
  print(paste(i, "-----> Clustering en cours"))
  dir.create(paste0(path,"4Cluster/"),recursive=T,showWarnings=F)
  for (distance in c("Pearson", "Spearman")){
    png(paste0(path,"4Cluster/",i,"_Matrice_",distance,".png"))
    # Choisir le mode de calcule des distances
    if (distance == "Pearson"){
      matDist = as.matrix(cor(data_tab))
      pheatmap(matDist, main = paste("Pheatmap Pearson", type, i))
      matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))

    }else if (distance == "Spearman"){
      matDist = as.matrix(cor(data_tab,method="spearman"))
      pheatmap(matDist, main = paste("Pheatmap Spearman", type, i))
      matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
    }
    dev.off()


    for (method in c("kmeans", "HCL")){
      print(paste(distance, method))
      png(paste0(path,"4Cluster/",i,"_Cluster_",method,"_",distance,".png"))
      Clustering(matDist = matDist,
                 nb_cluster = 5,
                 method = method,
                 titre = paste(type,i),
                 colors = color)
      dev.off()
    }
  }
  
  ##### Heatmap avant moyenne ####

  data_tab = as.matrix(data_tab)
  
  MyHeatmaps(path = paste0(path,"/Heatmap/"),data_tab, condition = i)
  ProfilsPNG(save_path = paste0(path,"/profils/"), data_tab, condition = i)
  ProfilsPDF(save_path = paste0(path,"/profils/"), data_tab, condition = i)
  
  # Heatmap avec calcul des moyennes

  mean_data_tab = MeanTabCalculation(data_tab, rnai_list, cluster,i)

  MyHeatmaps(paste0(path,"/Heatmap/"),mean_data_tab, moyenne = T, condition = i)
  ProfilsPNG(save_path = paste0(path,"/profils/"), mean_data_tab, moyenne = T, condition = i)
  ProfilsPDF(save_path = paste0(path,"/profils/"), mean_data_tab, moyenne = T, condition = i)
  
  # Heatmap sans log
  MyHeatmaps(paste0(path,"/HeatmapNoLog/"),data_tab, condition = i, Log = F)
  MyHeatmaps(paste0(path,"/HeatmapNoLog/"),mean_data_tab, moyenne = T, condition = i, Log = F)

}
  